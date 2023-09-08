#!/bin/bash

fna_file="third_stage.fna"
fastqc_file="SRR24502286.fastq"
alignment_file="alignment.sam"
sorted_alignment_file="sorted_alignment.sam"

# run fastqc on the file
echo ""
echo "Running fastqc on $fastqc_file"
echo ""
fastqc $fastqc_file
wait
mv *.html qcreport.html

# make reference for genome
echo ""
echo "Making reference for $fna_file"
echo ""
minimap2 -d ref.mmi $fna_file
wait

# index genome
echo ""
echo "Indexing genome for $fna_file"
echo ""
minimap2 -a ref.mmi $fastqc_file > $alignment_file
wait

# run samtools view
echo ""
echo "Running samtools view on $alignment_file"
echo ""
samtools view -bS $alignment_file > alignment.bam
wait

# make flagstat and put into flagstat.txt
echo ""
echo "Making flagstat for $alignment_file"
echo ""
samtools flagstat $alignment_file > flagstat.txt
wait

filename="flagstat.txt"
# extract percent
percent=$(grep -o -P '\d+\.\d+%' $filename | sed 's/%//')
# trim the first one
first=$(echo $percent | cut -f1 -d' ')
result=""

if (( $(awk 'BEGIN {print ("'"$first"'" > "90.0")}') )); then
 result="OK ✅✅✅"
else
 result="BAD ❌❌❌"
 echo "$result"
 echo "percent is too low"
 echo "Exiting..."
 exit 1
fi
wait

echo ""
echo "Percent is $first%"
echo "$result"
echo ""

# make samtools sort
echo ""
echo "Making samtools sort for $alignment_file"
echo ""
samtools sort -o $sorted_alignment_file $alignment_file
wait

# make faidx
echo ""
echo "Making faidx for $alignment_file"
echo ""
samtools faidx $fna_file
wait

# convert sam to bam
echo ""
echo "Converting sam to bam for $sorted_alignment_file"
echo ""
samtools view -bS -o alignment.bam $sorted_alignment_file
wait

# run freebayes
echo ""
echo "Running freebayes on alignment.bam"
echo ""
freebayes -f $fna_file alignment.bam > variants.vcf