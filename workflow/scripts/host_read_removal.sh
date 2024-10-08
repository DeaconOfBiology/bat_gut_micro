#!/usr/bin/env bash

################################################################
# author: TJ Rogers                                            #
# input: Name of the trimmed files (without the path)          #
# output: places cleaned reads into data/processed/clean_reads #
#                                                              #
# Removing host sequences.                                     #
################################################################

#### Create bowtie2 database of the host genome
bowtie2-build host_genome_sequence.fasta host_DB

#### Map reads against host sequence database, keeping both mapped and unmapped reads
bowtie2 -p 8 -x host_DB \
-1 SAMPLE_R1.fastq.gz \
-2 SAMPLE_R2.fastq.gz \
-S SAMPLE_mapped_and_unmapped.sam

#### Convert sam files to bam files to save space
samtools view -bS SAMPLE_mapped_and_unmapped.sam > SAMPLE_mapped_and_unmapped.bam

#### Collect reads not mapped to host sequences
samtools view -b -f 12 -F 256 \
    SAMPLE_mapped_and_unmapped.bam \
    > SAMPLE_bothReadsUnmapped.bam 

#### Sort bam files of unmapped reads
samtools sort -n -m 5G -@ 2 SAMPLE_bothReadsUnmapped.bam -o SAMPLE_bothReadsUnmapped_sorted.bam

#### Separate the paired-end reads into separate fastq files
samtools fastq -@ 8 SAMPLE_bothReadsUnmapped_sorted.bam \
  -1 SAMPLE_host_removed_R1.fastq.gz \
  -2 SAMPLE_host_removed_R2.fastq.gz \
  -0 /dev/null -s /dev/null -n