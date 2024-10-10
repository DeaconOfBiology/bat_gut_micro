import os
import subprocess
import argparse
from glob import glob

###########################################
# Make sure everything is set up properly #
###########################################
# The default length (-l) setting is 1000bp.
l=1000

#########################################
# Aligning reads to make coverage files #
#########################################
# Coverage file(s) for the binning programs
# I need to create coverage files for each binning program. This will be coverage
# files of each assembly accross all metagenomes. This will allow for better bin-
# ing.
def align_reads():
    print()

def make_cov_files():
    print()


#################################################################################
# Setting up functions to run each binning, refining, and dereplication program #
#################################################################################
# Function to run concoct
def run_concoct(assembly_path, output_dir):
    concoct_output_dir = os.path.join(output_dir, "concoct") #define concoct directory inside of ouput_dir
    os.makedirs(concoct_output_dir, exist_ok=True)           #create the concoct output directory

    cmd = f"""
    mamba run -n concoct concoct \
    --composition_file {assembly_path} \
    --coverage_file {assembly_path}.cov \
    -b {concoct_output_dir}
    """
    subprocess.run(cmd, shell=True, check=True)
    print(f"CONCOCT binning completed for {assembly_path}")

# Function to run metabat2
def run_metabat2(assembly_path, output_dir):
    metabat2_output_dir = os.path.join(output_dir, "metabat2")
    os.makedirs(metabat2_output_dir, exist_ok=True)

    cmd = f"""
    mamba run -n metabat2 metabat2 \
    -i {assembly_path} \
    -o {metabat2_output_dir}
    """
    subprocess.run(cmd, shell=True, check=True)
    print(f"Metabat2 binning completed for {assembly_path}")

# Function to run maxbin2
def run_maxbin2(assembly_path, output_dir):
    maxbin2_output_dir = os.path.join(output_dir, "maxbin2")
    os.makedirs(maxbin2_output_dir,exist_ok=True)

    cmd = f"""
    mamaba run -n maxbin2 maxbin2 \
    -contig {assembly_path} \
    -out {maxbin2_output_dir}/maxbin2_output \
    -abund {assembly_path}.abund
    """
    subprocess.run(cmd, shell=True, check=True)
    print(f"MaxBin2 binning completed for {assembly_path}")

def run_bin_refinner():
    print()

def run_drep():
    print()