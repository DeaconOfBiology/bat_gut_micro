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
def align_reads(assembly_path, reads, output_dir):
    os.makedirs(ouput_dir, exist_ok=True)

    outR1 = f"map_{os.path.basename(reads[0])}"
    outR2 = f"map_{os.path.basename(reads[1])}"
    if not outR1.endswith(".gz"):
        outR1 += ".gz"
    if not outR2.endswith(".gz"):
        outR2 += ".gz"
    outfile = [os.path.join(path, outR1), os.path.join(path, outR2)] #Need to figure out what this is doing

    # Iterate over each sample and map the reads
        for sample in samples:
            sample_name = sample['sample_name']
            forward_reads = sample['forward']
            reverse_reads = sample['reverse']

            # Define output paths
            sam_file = os.path.join(assembly_output_dir, f"{sample_name}.sam")
            sorted_bam_file = os.path.join(assembly_output_dir, f"{sample_name}.sorted.bam")

            # Step 1: Run BWA-MEM to map reads and generate SAM file
            bwa_mem_cmd = f"""
            mamba run -n bwa bwa mem {assembly} {forward_reads} {reverse_reads} > {sam_file}
            """
            subprocess.run(bwa_mem_cmd, shell=True, check=True)

            # Step 2: Convert SAM to BAM and sort it
            sam_to_bam_cmd = f"""
            mamba run -n samtools samtools view -bS {sam_file} | mamba run -n samtools samtools sort -o {sorted_bam_file}
            """
            subprocess.run(sam_to_bam_cmd, shell=True, check=True)

            # Step 3: Index the sorted BAM file
            index_bam_cmd = f"mamba run -n samtools samtools index {sorted_bam_file}"
            subprocess.run(index_bam_cmd, shell=True, check=True)

            # Step 4: Remove intermediate SAM file
            os.remove(sam_file)
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