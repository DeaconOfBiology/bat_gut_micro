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
input=(${snakemake_input[accessions]})
output=(${snakemake_output[ncbi_zip]})

#### Running datasets to retrieve sequences
datasets download genome accession --inputfile "$input" \
    --include gff3,rna,cds,protein,genome,seq-report \
    --filename "$output"




##########################################################
# This code chunk was originally used to test the script #
##########################################################
# datasets download genome accession --inputfile /projects/raw_lab/projects/Bats/bat_gut_micro/data/references/test_bat_host_accessions.txt \
#     --include gff3,rna,cds,protein,genome,seq-report --filename /projects/raw_lab/projects/Bats/bat_gut_micro/test/bat_references.zip

# unzip /projects/raw_lab/projects/Bats/bat_gut_micro/test/bat_references.zip /projects/raw_lab/projects/Bats/bat_gut_micro/test/
