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
9. Make samtoos filtering (in bam format)
    ```bash
    bamtools filter -in sorted_alignments.bam -out sorted_alighments.filtered.bam -mapQuality '>30'

    ```
10. Convert sam to bam
    ```bash
    samtools view -bS -o alignment.bam alignment.sam

    # sort again, but for bam
    samtools sort -o sorted_alignments.bam alighment.bam
    # filtering for bam
    bamtools filter -in sorted_alighments.bam -out sorted_alighments.filtered.bam -mapQuality '>30'
    ```
11. Run freebayes
    ```bash
    freebayes -f third_stage.fna sorted_alignments.filtered.bam > output.vcf
    ```
