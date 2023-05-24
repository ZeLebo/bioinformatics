# Result of work

# Table of contents
1. [Working in command line](#pipeline-with-arms--ruchki)
2. [Writing bash script for list of commands](#writing-the-bash-script-to-make-the-commands-above)
3. [Working with Dagster as pipeline creator tool]
4. [How to make a new project and run it]
5. [How to visualize the work]

### References
   E.coli link - [link to run](https://www.ncbi.nlm.nih.gov/sra/SRX20287202[accn])

### Mnogo files
![alt-text](mnogo_fails.jpg)

### Pipeline with arms (ruchki)

1. Download all the data
2. Run the fastqc command on fastqc file
    ```bash
    fastqc SRR24502286.fastq
    mv *.html qcreport.html
    ```
3. Index the genome
    ```bash
    minimap2 -d <prefix.mmi> <reference.fasta>
    minimap2 -d ref.mmi third_stage.fna

    # indexing with reference
    minimap2 -a ref.mmi SRR24502286.fastq > alignment.sam
    ```
4. Make samtools view
    ```bash
    samtools view alignment.sam
    ```
5. Make a txt file for flagstat.txt
    ```bash
    samtools flagstat alignment.sam > flagstat.txt
    ```
6. Write a bash script to extract percent from the file
    ```bash

    #!/bin/bash

    # Get the filename from the command line argument
    filename=$1

    # Use grep and sed to extract the percent value
    percent=$(grep -o -P '\d+\.\d+%' $filename | sed 's/%//')

    # trim the first one
    first=$(echo $percent | cut -f1 -d' ')

    echo "Percent is $first%"

    result="BAD"

    if (( $(awk 'BEGIN {print ("'$first'" > "90.0")}') )); then 
        result="OK ✅✅✅"
    else
        result="BAD ❌❌❌"
    fi

    echo $result

    ```
7. Make samtools sort
    ```bash
    samtools sort -o sorted_alignment.sam alignment.sam
    ```
8. Make samtoos faidx
    ```bash
    samtools faidx third_stage.fna
    ```
9. Convert sam to bam
   ```bash
   samtools view -bS -o alignment.bam alignment.sam
   # sort again, but for bam
   samtools sort -o sorted_alignments.bam alighment.bam
   ```
10. Run freebayes
    ```bash
    freebayes -f third_stage.fna sorted_alignments.filtered.bam > output.vcf
    ```
    
#### Result of flagstat
   ```text
    15644 + 0 in total (QC-passed reads + QC-failed reads)
    13760 + 0 primary
    1269 + 0 secondary
    615 + 0 supplementary
    0 + 0 duplicates
    0 + 0 primary duplicates
    14907 + 0 mapped (95.29% : N/A)
    13023 + 0 primary mapped (94.64% : N/A)
    0 + 0 paired in sequencing
    0 + 0 read1
    0 + 0 read2
    0 + 0 properly paired (N/A : N/A)
    0 + 0 with itself and mate mapped
    0 + 0 singletons (N/A : N/A)
    0 + 0 with mate mapped to a different chr
    0 + 0 with mate mapped to a different chr (mapQ>=5)
   ```
### Writing the bash script to make the commands above

I've writen a simple bash script to run the commands above -> [bash script](./bash_script/manual.sh).

Need to mention: before starting the script need to install packages to linux system
```bash
    sudo apt-get install samtools, minimap2, fastqc
```

### Creating new project with DugsterIO
In order to start a new project we need to install the dagster on your system:
```bash
    pip3 install dagster dagit
```
After this we can create new project with Dagster
```bash
    dagster project scaffold --name <name-of-your-project>
```
Then we can navigate to the new-created folder and find setup.py file:
```python
from setuptools import find_packages, setup

setup(
    name="dagster_pipeline",
    packages=find_packages(exclude=["dagster_pipeline_tests"]),
    install_requires=[
        "dagster",
        "dagster-cloud"
    ],
    extras_require={"dev": ["dagit", "pytest"]},
)
```
in section `install_requires=[]` we can add our needed packages for further working.

To install the dependencies run the following command:
```bash
    pip3 install -e ".[dev]"
```
### First program to print `hello world`
Now let's create a new file `say_hello_pipeline.py` in <name-of-your-project> folder and put this as content:
```python
from dagster import job

@job
def say_hello():
    print(f"Hello world from Dagster!")

if __name__ == '__main__':
    say_hello()
```
And run this file as follows:
```bash
  dagster job execute -f <name-of-your-project>/say_hello_pipeline.py > output.txt
```
And here's the output file: [output.txt](dagster-pipeline/output.txt)

Need to mention that you cannot pass `None` as default param in your function, so instead of
```python
def foo(name = None): # it will fail at start due to NoneClass Exception
    if name is None:
      name = "aaa"
    ...
```
write this
```python
def foo(name="aaa"):
    ...
```
### The needed program
You can see the code of `pipeline` here -> [code](dagster-pipeline/dagster_pipeline/pipeline.py)

The result of work is here -> [result in vcf format](dagster-pipeline/dagster_pipeline/output.vcf)

The logs can be found here -> [logs of running command](dagster-pipeline/dagster_pipeline/logs.txt)

### Visualisation of pipeline
You can visualize the pipeline when running command:
```bash
   dagit -f <name-of-your-file>.py
```
I have this result:
![photo](dagster_photo.jpg)
(yeap, dugster creating new description for job itself)

As always, dagster will build up the graph based on the sequence of parameters you're passing in your code.
So, if you write code similar to this:
```python
    retrieve_files()
    fastqc()
    rename_output()
    make_reference()
    index_genome()
    samtools_view()
    make_flagstat()
```
You will end up with only nodes of graph that are not connected to each other in any way,
in order to get the links you need to specify the output data and input data, so the code would be smth like this:
```python
    data = retrieve_files()
    data1 = fastqc(data)
    rename_output(data1)
    data2 = make_reference(data1)
    alignment_sam = index_genome(data2)
    bam_alignment = samtools_view(alignment_sam)
    flagstat_file = make_flagstat(bam_alignment)
    ok, not_ok = check_percent(flagstat_file)
    ok1 = print_ok(ok)
    print_not_ok(not_ok)
```

So the visualisation of pipeline is basically a creativity task (by default all python function return None, so you can pass this `None` to the next function and get the beautiful graph)