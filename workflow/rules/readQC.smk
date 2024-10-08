configfile:"/projects/raw_lab/projects/Bats/bat_gut_micro/config/configure.yml"

rule all:
    input:
        final_reads=expand([config["clean_reads"] + "{sample}-sg_1.fastq.gz", config["clean_reads"] + "{sample}-sg_1.fastq.gz"], sample=config["Samples"]),
        ncbi_zip=config["references"] + "bat_references.zip"


rule clean_raw_reads:
    input:
        r1=config["raw_reads"] + "raw-{sample}-sg_R1_001.fastq.gz",
        r2=config["raw_reads"] + "raw-{sample}-sg_R2_001.fastq.gz"
    output:
        tr1=temp(config["temp_files"] + "temp_{sample}-sg_1.fastq.gz"),
        tr2=temp(config["temp_files"] + "temp_{sample}-sg_2.fastq.gz"),
        fr1=config["clean_reads"] + "{sample}-sg_1.fastq.gz", 
        fr2=config["clean_reads"] + "{sample}-sg_2.fastq.gz"
    threads: 5
    log: 
        config["logs"] + "bbduk_{sample}-sg.log"
    params:
        script=config["scr"] + "trimming.sh"
    conda:
        config["envs"] + "bbmap.yaml"
    shell:
        """
            (bash {params.script} {input.r1} {input.r2} {output.tr1} {output.tr2} {output.fr1} {output.fr2} ) 2> {log}
        """


rule get_host_references:
    input:
        accessions=[config["references"] + "bat_host_accessions.txt"]
    output:
        ncbi_zip=config["references"] + "bat_references.zip"
    threads: 5
    log:
        config["logs"] + "bat_host_accessions.log"
    params:
        script=config["scr"] + "retreve_host_references.sh"
    conda:
        config["envs"] + "ncbi_datasets.yaml"
    shell:
        """
            (bash {params.script} {input.accessions} {output.ncbi_zip}) 2> {log}
        """
        

# rule remove_host_reads:
#     input:
#         ncbi_zip=config["references"] + "bat_references.zip",
#         r1=config["clean_reads"] + "{sample}-sg_1.fastq.gz", 
#         r2=config["clean_reads"] + "{sample}-sg_2.fastq.gz"
#     output:
#         r1=config["clean_reads"] + "{sample}-sg_host_removed_2.fastq.gz", 
#         r2=config["clean_reads"] + "{sample}-sg_host_removed_2.fastq.gz"
#     threads: 15
#     log:
#      config["logs"] + ""
#     params:
#         script=config[""] + "host_read_removal.sh",
#         host=cofig["references"] + "ncbi_dataset/data/{host}/{host}_mPhyDis1.pri.v3_genomic.fna" #May need to change this ending.
#     conda:
#         config["envs"] + ""
#     shell:
#         """
#           (bash {params.script} {input.r1} {input.r2} {output.r1} {output.r2} {params.host})2>{log}
#         """