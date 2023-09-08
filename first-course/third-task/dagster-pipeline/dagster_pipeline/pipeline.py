from dagster import job, op, Output, Out
import subprocess


@op
def retrieve_files() -> dict:
    fna_file = "third_stage.fna"
    fastq_file = 'SRR24502286.fastq'
    ref_file = 'ref.mmi'
    return {
        'fna': fna_file,
        'fastq': fastq_file,
        'ref': ref_file
    }


@op
def rename_output(context, data) -> None:
    subprocess.run(['mv', '*.html', 'qcreport.html'])


@op
def fastqc(context, data: dict) -> dict:
    context.log.info("running fastqc on fastq file")
    subprocess.run(['fastqc', data['fastq']])
    return data


@op
def make_reference(context, data: dict) -> dict:
    context.log.info("making reference for fna file")
    subprocess.run(['minimap2', '-d', data['ref'], data['fna']])
    return data


@op
def index_genome(context, data: dict, output_file: str = 'alignment.sam') -> str:
    context.log.info("indexing the genome")
    subprocess.run(['minimap2', '-a', data['ref'], data['fastq']], stdout=open(output_file, 'w'))
    return output_file


@op
def samtools_view(context, alignment_file: str = 'alignment.sam', output_file: str = 'alignment.bam') -> str:
    context.log.info("making samtools view on alignment file")
    subprocess.run(['samtools', 'view', '-bS', alignment_file], stdout=open(output_file, 'w'))
    return output_file


@op
def make_flagstat(context, alignment_file: str = 'alignment.sam', output_file = 'flagstat.txt') -> str:
    context.log.info("making file for samtools flagstat")
    subprocess.run(['samtools', 'flagstat', alignment_file], stdout=open(output_file, 'w'))
    return output_file


@op(out={
    'ok_result': Out(is_required=False),
    'bad_result': Out(is_required=False)
}
)
def check_percent(context, filename: str = 'flagstat.txt') -> Output:
    context.log.info("checking percent")
    with open(filename, 'r') as file:
        contents = file.read()
        percent = float(contents.split()[0].replace('%', ''))
    if percent > 90.0:
        result = "OK ✅✅✅"
        context.log.info(f"percent = {percent}% - {result}")
        yield Output(1, 'ok_result')
    else:
        result = "BAD ❌❌❌"
        context.log.info(f"percent = {percent}% - {result}")
        yield Output(2, 'bad_result')


@op
def print_ok(context, ok):
    context.log.info(ok)
    return ok


@op
def print_not_ok(context, not_ok):
    context.log.info(not_ok)
    return not_ok


@op
def make_samtools_sort(context, ok: Output, alignment_file: str="alignment.sam") -> str:
    context.log.info("sorting alignment file")
    subprocess.run(['samtools', 'sort', '-o', 'sorted_alignment.sam', alignment_file])
    return "sorted_alignment.sam"


@op
def make_faidx(context, alignment_file, data=None) -> None:
    if data is None:
        data = {'fna': 'third_stage.fna'}

    context.log.info("making fai file from fna_file")
    subprocess.run(['samtools', 'faidx', data['fna']])


@op
# tmp is operand to show in graph
def convert_sam_to_bam(context, sorted_alignment_file="sorted_alignment.sam", tmp = None) -> str:
    context.log.info("converting sam to bam format")
    subprocess.run(['samtools', 'view', '-bS', '-o', 'alignment.bam', sorted_alignment_file])
    return 'alignment.bam'


@op
def freebayes_run(context, bam_alignment='alignment.bam', data=None, output_file='result.vcf') -> str:
    if data is None:
        data = {'fna': 'third_stage.fna'}
    context.log.info("running freebayes")
    subprocess.run(['freebayes', '-f', data['fna'], bam_alignment], stdout=open(output_file, 'w'))
    return output_file


@job
def run_job():
    data = retrieve_files()
    data1 = fastqc(data)
    rename_output(data1)
    data2 = make_reference(data1)
    alignment_sam = index_genome(data2)
    bam_alignment = samtools_view(alignment_sam)
    flagstat_file = make_flagstat(alignment_sam)
    ok, not_ok = check_percent(flagstat_file)
    ok1 = print_ok(ok)
    print_not_ok(not_ok)

    sorted_alignment = make_samtools_sort(ok1, alignment_sam)
    # res = make_faidx(sorted_alignment)
    bam_alignment = convert_sam_to_bam(sorted_alignment)
    result = freebayes_run(bam_alignment)
    print(result)
