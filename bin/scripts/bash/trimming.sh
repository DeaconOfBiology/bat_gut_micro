#!/usr/bin/env bash

# author: TJ Rogers
# input: Name of the raw files (without the path)
# output: places cleaned reads into data/processed/clean_reads
#
# Trimming the raw reads.



# Cleaning up sample names for input
raw_reads=$(basename ${$1} _R1_001.fastq.gz)
reads=${raw_reads#*-}

# Directories
tmp_dir="data/tmp"
clean_dir="data/processed/clean_reads"
mkdir "$clean_dir"


"$bbduk"/bbduk.sh -Xmx1g in1="$raw_reads"_R1_001.fastq.gz \
                in2="$raw_reads"_R2_001.fastq.gz \
                out1="$tmp_dir"/tmp_"$reads"_1.fastq.gz \
                out2="$tmp_dir"/tmp_"$reads"_2.fastq.gz \
                ref="$bbduk"/adapters.fa ktrim=r k=21 \
                qin=auto mink=11 hdist=2 tbo tpe



"$bbduk"/bbduk.sh  -Xmx1g in1="$tmp_dir"/tmp_"$reads"_1.fastq.gz \
                in2="$tmp_dir"/tmp_"$reads"_2.fastq.gz \
                out1="$clean_dir"/clean_"$reads"_1.fastq.gz \
                out2="$clean_dir"/clean_"$reads"_2.fastq.gz \
                ref={params.ref} k=27 hdist=1 qtrim=rl \
                trimq=17 cardinality=t mingc=0.05 maxgc=0.95