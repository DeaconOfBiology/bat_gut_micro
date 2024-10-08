#!/usr/bin/env bash

################################################################
# author: TJ Rogers                                            #
# input: Name of the trimmed files (without the path)          #
# output: places cleaned reads into data/processed/clean_reads #
#                                                              #
# Removing host sequences.                                     #
################################################################

#Assigning arguments:
sample=$(basename ${$1} -sg_1.fastq.gz) #Grabing the sample name
r1=$1                                   #Forward read input
r2=$2                                   #Reverse read input
fr1=$3                                  #Forward read output
fr2=$4                                  #Reverse read output
host_genome=$5                          #Host genome to index

#### Create bowtie2 database of the host genome
bowtie2-build "$host_genome" host_DB

#### Map reads against host sequence database, keeping both mapped and unmapped reads
bowtie2 -p 8 -x host_DB \
-1 "$r1" \
-2 "$r2" \
-S "$sample"_mapped_and_unmapped.sam

#### Convert sam files to bam files to save space
samtools view -bS "$sample"_mapped_and_unmapped.sam > "$sample"_mapped_and_unmapped.bam

#### Collect reads not mapped to host sequences
samtools view -b -f 12 -F 256 \
    "$sample"_mapped_and_unmapped.bam \
    > "$sample"_bothReadsUnmapped.bam 

#### Sort bam files of unmapped reads
samtools sort -n -m 5G -@ 2 "$sample"_bothReadsUnmapped.bam -o "$sample"_bothReadsUnmapped_sorted.bam

#### Separate the paired-end reads into separate fastq files
samtools fastq -@ 8 "$sample"_bothReadsUnmapped_sorted.bam \
  -1 "$fr1" \
  -2 "$fr2" \
  -0 /dev/null -s /dev/null -n

#### Cleaning up tmp files
