rule all:


rule clean_raw_reads:
    input:
        reads=["raw-{sample}_R1_001.fastq.gz","raw-{sample}_R2_001.fastq.gz"]
        script = "bin/scripts/bash/trimming.sh"
    output:
        tem_reads=[temp(config["temp_files"] + "/temp_{sample}-sg_1.fastq.gz"),temp(config["temp_files"] + "/temp_{sample}-sg_1.fastq.gz")],
        final_reads=["{sample}_1.fastq.gz","{sample}_1.fastq.gz"]  
    threads:
    log: 
        config["logs"] + "/bbduk_{sample}"
    conda:
        "envs/bbmap.yaml"
    script:
        ""
    
