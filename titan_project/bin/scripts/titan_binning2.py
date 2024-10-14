#!/usr/bin/env python3

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
def align_reads(read_path, assembly, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Find forward and reverse reads
    R1 = sorted(glob(os.path.join(read_path, '**', "*_R1_001.fastq.gz"), recursive=True))
    R2 = sorted(glob(os.path.join(read_path, '**', "*_R2_001.fastq.gz"), recursive=True))

    # Check that forward and reverse basenames match
    if sorted(os.path.basename(f).replace('_R1_001.fastq.gz', '') for f in R1) != \
       sorted(os.path.basename(r).replace('_R2_001.fastq.gz', '') for r in R2):
        raise ValueError("Mismatch between forward and reverse reads. Ensure that each forward read has a matching reverse read.")
    
    reads = list(zip(R1, R2))
    print(f"Indexing assembly: {assembly}")
    assembly_name = os.path.basename(os.path.dirname(os.path.dirname(assembly)))
    assembly_output_dir = os.path.join(output_dir, assembly_name)
    os.makedirs(assembly_output_dir, exist_ok=True)
    ref_index= os.path.splitext(os.path.basename(assembly))[0]
    ref_index_path=os.path.join(assembly_output_dir, f"{ref_index}")
    subprocess.run(f"mamba run -n bwa bwa index -p {ref_index_path} {assembly}", shell=True, check=True)


    for i, (forward, reverse) in enumerate(reads):
        sample_name = "_".join(os.path.basename(forward).split("_")[:-2])
        sam_output = os.path.join(assembly_output_dir, f"{sample_name}.sam")
        unbam_output = os.path.join(assembly_output_dir, f"{sample_name}.unsorted.bam")
        bam_output = os.path.join(assembly_output_dir, f"{sample_name}.sorted.bam")
        # Align reads with BWA and convert to BAM using samtools
        print(f"Aligning {forward} and {reverse} to {ref_index_path}")
        bwa_cmd = f"mamba run -n bwa bwa mem -t 15 {ref_index_path} {forward} {reverse} > {sam_output}"
        samtools_cmd = f"mamba run -n samtools samtools view -F 2316 -b -o {unbam_output} {sam_output}"

        # Combine commands
        print(f"Sorting BAM file: {unbam_output}")
        sort_cmd = f"mamba run -n samtools samtools sort {unbam_output} -o {bam_output}"

        print(f"Indexing BAM file: {bam_output}")
        index_cmd = f"mamba run -n samtools samtools index {bam_output}"
        
        # Execute the commands
        subprocess.run(bwa_cmd, shell=True, check=True)
        subprocess.run(samtools_cmd, shell=True, check=True)
        subprocess.run(sort_cmd, shell=True, check=True)
        subprocess.run(index_cmd, shell=True, check=True)

        # Remove unsorted BAM file to save space
        if os.path.exists(f"{unbam_output}"):
            os.remove(f"{unbam_output}")
            os.remove(sam_output)



# #################################################################################
# # Setting up functions to run each binning, refining, and dereplication program #
# #################################################################################
# # Function to run concoct
# def run_concoct(assembly_path, output_dir):
#     concoct_output_dir = os.path.join(output_dir, "concoct") #define concoct directory inside of output_dir
#     os.makedirs(concoct_output_dir, exist_ok=True)           #create the concoct output directory

#     # Find all assembly files


#     # Step 1: Cut contigs into smaller fragments:
#     cut_cmd=f"""mamba run -n concoct cut_up_fasta.py original_contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa"""
#     # Step 2: Generate coverage depth table of assembly across samples
#     cov_cmd=f"""mamba run -n concoct concoct_coverage_table.py contigs_10K.bed mapping/Sample*.sorted.bam > coverage_table.tsv"""
#     # Step 3: Run Concoct
#     con_cmd = f"""
#     mamba run -n concoct concoct \
#         --composition_file {assembly_path} \
#         --coverage_file {assembly_path}.cov \
#         -b {concoct_output_dir}
#     """
#     # Step 4: Merge subcontig clustering into original contig clustering
#     merge_cmd=f"""mamba run -n concoct merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv"""
#     # Step 5: Extract bins as individual FASTA files:
#     mkdir concoct_output/fasta_bins
#     bin_cmd=f"""mamba run -n concoct extract_fasta_bins.py original_contigs.fa concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins"""
    
#     # Execute the commands
#     subprocess.run(cut_cmd, shell=True, check=True)
#     subprocess.run(cov_cmd, shell=True, check=True)
#     subprocess.run(con_cmd, shell=True, check=True)
#     subprocess.run(merge_cmd, shell=True, check=True)
#     subprocess.run(bin_cmd, shell=True, check=True)

#     cmd = f"""
#     mamba run -n concoct concoct \
#     --composition_file {assembly_path} \
#     --coverage_file {assembly_path}.cov \
#     -b {concoct_output_dir}
#     """
#     subprocess.run(cmd, shell=True, check=True)
#     print(f"CONCOCT binning completed for {assembly_path}")

# # Function to run metabat2
# def run_metabat2(assembly_path, output_dir):
#     metabat2_output_dir = os.path.join(output_dir, "metabat2")
#     os.makedirs(metabat2_output_dir, exist_ok=True)

#     cmd = f"""
#     mamba run -n metabat2 metabat2 \
#     -i {assembly_path} \
#     -o {metabat2_output_dir}
#     """
#     subprocess.run(cmd, shell=True, check=True)
#     print(f"Metabat2 binning completed for {assembly_path}")

# # Function to run maxbin2
# def run_maxbin2(assembly_path, output_dir):
#     maxbin2_output_dir = os.path.join(output_dir, "maxbin2")
#     os.makedirs(maxbin2_output_dir,exist_ok=True)

#     cmd = f"""
#     mamaba run -n maxbin2 maxbin2 \
#     -contig {assembly_path} \
#     -out {maxbin2_output_dir}/maxbin2_output \
#     -abund {assembly_path}.abund
#     """
#     subprocess.run(cmd, shell=True, check=True)
#     print(f"MaxBin2 binning completed for {assembly_path}")

# def run_bin_refinner():
#     print()

# def run_drep():
#     print()

def main():
    # Command-line argument parsing
    parser = argparse.ArgumentParser(description="Test read alignment")

    # Add flags for fastq and assembly directories
    parser.add_argument("-f", "--fastq", required=True, type=str, help="Path to the directory with fastq files")
    parser.add_argument("-a", "--assembly", required=True, type=str, help="Path to the directory with assembly files")
    parser.add_argument("-o", "--output", required=True, type=str, help="Directory to save BAM outputs")

    args = parser.parse_args()


    # Ensure output directory exists
    os.makedirs(args.output, exist_ok=True)

    # Get all assembly files matching the specified pattern
    assembly_files = glob(os.path.join(args.assembly, '**', "final.contigs.fa"), recursive=True)
    # Run the selected binning programs for each assembly file
    for assembly in assembly_files:
        # Extract the assembly name (without extension) to create output subdirectories
        assembly_name = os.path.splitext(os.path.basename(assembly))[0]
        assembly_output_dir = os.path.join(args.output, assembly_name)
        align_reads(args.fastq, assembly, assembly_output_dir)


# def main():
#     # Command-line argument parsing
#     parser = argparse.ArgumentParser(description="Run metagenomic binning tools")
    
#     # Add flags for assembly and output directories
#     parser.add_argument("-a", "--assembly", required=True, type=str, help="Path to the assembly files (glob pattern)")
#     parser.add_argument("-o", "--output", required=True, type=str, help="Directory to save binning outputs")
    
#     # Add flags for selecting binning tools
#     parser.add_argument("--concoct", action="store_true", help="Run CONCOCT binning")
#     parser.add_argument("--maxbin2", action="store_true", help="Run MaxBin2 binning (default if no other specified)")
#     parser.add_argument("--metabat2", action="store_true", help="Run MetaBAT2 binning")
    
#     args = parser.parse_args()

#     # Ensure output directory exists
#     os.makedirs(args.output, exist_ok=True)
    
#     # Default to MaxBin2 if no binning tool is selected
#     if not (args.concoct or args.metabat2 or args.maxbin2):
#         args.maxbin2 = True  # Set MaxBin2 as the default option
    
#     # Get all assembly files matching the specified pattern
#     assembly_files = glob(args.assembly)
    
#     # Run the selected binning programs for each assembly file
#     for assembly in assembly_files:
#         # Extract the assembly name (without extension) to create output subdirectories
#         assembly_name = os.path.splitext(os.path.basename(assembly))[0]
#         assembly_output_dir = os.path.join(args.output, assembly_name)
        
#         if args.concoct:
#             run_concoct(assembly, assembly_output_dir)
#         if args.maxbin2:
#             run_maxbin2(assembly, assembly_output_dir)
#         if args.metabat2:
#             run_metabat2(assembly, assembly_output_dir)

if __name__ == "__main__":
    main()