#!/bin/bash

# Assuming the script receives five arguments:
# 1. Genome index
# 2. FASTQ file 1
# 3. FASTQ file 2
# 4. Output BAM file
# 5. Output log file

genome_index=$1
fastq1=$2
fastq2=$3
output_bam=$4
output_log=$5

# Your alignment command, redirecting stderr to the log file
bowtie2 -x "$genome_index" -1 "$fastq1" -2 "$fastq2" | samtools view -bS -> "$output_bam" 2> "$output_log"
