configfile:"/projects/raw_lab/projects/Bats/bat_gut_micro/config/configure.yml"

rule all:
    input:
        final_reads=expand([config["temp_files"] + "{sample}-sg_1.fastq.gz", config["temp_files"] + "{sample}-sg_1.fastq.gz"], sample=config["Samples"]),
        ncbi_zip=config["references"] + "bat_references.zip"

rule clean_raw_reads:
    input:
        reads=[config["raw_reads"] + "raw-{sample}-sg_R1_001.fastq.gz", config["raw_reads"] + "raw-{sample}-sg_R2_001.fastq.gz"]
    output:
        final_reads=[config["temp_files"] + "{sample}-sg_1.fastq.gz", config["temp_files"] + "{sample}-sg_1.fastq.gz"]
    threads: 5
    log: 
        config["logs"] + "bbduk_{sample}-sg"
    conda:
        config["envs"] + "bbmap.yaml"
    script:
        config["scr"] + "trimming.sh"

rule get_host_references:
    input:
        accessions=[config["references"] + "bat_host_accessions.txt"]
    output:
        ncbi_zip=config["references"] + "bat_references.zip"
    threads: 5
    log:
        config["logs"]
    conda:
        config["envs"] + "ncbi_datasets.yaml"
    script:
        config["scr"] + "retreve_host_references.sh"