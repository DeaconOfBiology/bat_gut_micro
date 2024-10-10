#!/usr/bin/env bash

################################################################
# author: TJ Rogers                                            #
# input: Name of the clean files (without the path)            #
# output: places filterd cleaned reads into                    #
# data/processed/clean_reads                                   #
#                                                              #
# Removing host reads from clean reads.                        #
################################################################


###########################################################
# This code block dowloads the genomes for representative #
# sequences using the datasets program.                   #
###########################################################
#### Calling in input and output from the rule doc:
accession=$1
ncbi_zip=$2
ncbi_unziped=$3
refdir=$4
#### Running datasets to retrieve sequences
datasets download genome accession --inputfile "$accession" \
    --include gff3,rna,cds,protein,genome,seq-report \
    --filename "$ncbi_zip"

unzip "$ncbi_zip" -d "$refdir"

mv "$refdir"ncbi_dataset/data/GCA_036768555.1/GCA_036768555.1_ASM3676855v1_genomic.fna "$refdir"ncbi_dataset/data/GCA_036768555.1/GCA_036768555.1.fna
mv "$refdir"ncbi_dataset/data/GCA_028538775.1/GCA_028538775.1_mMyoYum1.0.hap1_genomic.fna "$refdir"ncbi_dataset/data/GCA_028538775.1/GCA_028538775.1.fna
mv "$refdir"/ncbi_dataset/data/GCA_004126475.3/GCF_004126475.2_mPhyDis1.pri.v3_genomic.fna "$refdir"ncbi_dataset/data/GCA_004126475.3/GCF_004126475.3.fna
mv "$refdir"/ncbi_dataset/data/GCA_019186645.2/GCF_019186645.2_TTU_PhHast_1.1_genomic.fna "$refdir"ncbi_dataset/data/GCA_019186645.2/GCF_019186645.2.fna
##########################################################
# This code chunk was originally used to test the script #
##########################################################
# datasets download genome accession --inputfile "$refdir"/test_bat_host_accessions.txt \
#     --include gff3,rna,cds,protein,genome,seq-report --filename /projects/raw_lab/projects/Bats/bat_gut_micro/test/bat_references.zip

# unzip /projects/raw_lab/projects/Bats/bat_gut_micro/test/bat_references.zip /projects/raw_lab/projects/Bats/bat_gut_micro/test/
