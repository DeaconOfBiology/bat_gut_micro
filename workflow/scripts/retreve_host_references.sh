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
# datasets download genome accession --inputfile snakemake_input[accessions] \
#     --include gff3,rna,cds,protein,genome,seq-report

datasets download genome accession --inputfile /projects/raw_lab/projects/Bats/bat_gut_micro/data/references/bat_host_accessions.txt \
    --include gff3,rna,cds,protein,genome,seq-report