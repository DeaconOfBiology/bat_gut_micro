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

#####################################
# Call read files from the rule doc #
#####################################

#####################
# Removing adapters #
#####################
# Check if final output files exist
if [[ -f "$5" && -f "$6" ]]; then
  echo "Final files $5 and $6 already exist. Skipping this sample..."
  exit 0
fi

# Check if intermediate files exist
if [[ -f "$3" && -f "$4" ]]; then
    echo "Intermediate files $3 and $4 found. Running step 2 to remove phiX sequences..."
    # Removing phiX
    bbduk.sh -Xmx1g in1=$1 in2=$2 out1=$3 out2=$4 ref=/users/troger50/.conda/envs/bbmap/share/bbmap/resources/adapters.fa ktrim=r k=21 qin=auto mink=11 hdist=2 tbo tpe
else
    echo "Intermediate files not found. Running both steps to remove adapter and phiX sequences..."
    # Removing adapters
    bbduk.sh -Xmx1g in1=$1 in2=$2 out1=$3 out2=$4 ref=/users/troger50/.conda/envs/bbmap/share/bbmap/resources/adapters.fa ktrim=r k=21 qin=auto mink=11 hdist=2 tbo tpe
    # Removing phiX
    bbduk.sh -Xmx1g in1=$3 in2=$4 out1=$5 out2=$6 ref=/users/troger50/.conda/envs/bbmap/share/bbmap/resources/phix174_ill.ref.fa.gz k=27 hdist=1 qtrim=rl trimq=17 cardinality=t mingc=0.05 maxgc=0.95
fi