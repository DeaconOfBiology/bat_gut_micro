#!/usr/bin/env bash

################################################################
# author: TJ Rogers                                            #
# input: Name of the raw files (without the path)              #
# output: places trimed reads into data/temp                   #
#                                                              #
# Trimming the raw reads.                                      #
################################################################

###################################################
# Send all stderr from this script to the log file #
###################################################
# exec 2> "${snakemake_log[0]}" 

#####################################
# Call read files from the rule doc #
#####################################
# Raw
# raw_reads=(${snakemake_input[reads]})
# rr1="${raw_reads[0]}"
# rr2="${raw_reads[1]}"


# # Temp
# temp_reads=(${snakemake_output[tem_reads]})
# tr1="${temp_reads[0]}"
# tr2="${temp_reads[1]}"

# # Clean
# clean_reads=(${snakemake_output[final_reads]})
# cr1="${clean_reads[0]}"
# cr2="${clean_reads[1]}"
#####################
# Removing adapters #
#####################
bbduk.sh -Xmx1g in1=$1 in2=$2 out1=$3 out2=$4 ref=/users/troger50/.conda/envs/bbmap/share/bbmap/resources/adapters.fa ktrim=r k=21 qin=auto mink=11 hdist=2 tbo tpe

#################
# Removing phiX #
#################
bbduk.sh  -Xmx1g in1=$3 in2=$4 out1=$5 out2=$6 ref=/users/troger50/.conda/envs/bbmap/share/bbmap/resources/phix174_ill.ref.fa.gz k=27 hdist=1 qtrim=rl trimq=17 cardinality=t mingc=0.05 maxgc=0.95