#
# rule bwa_index:
#     input:
#         index=config['alignment']['genome'] + "GRCh38.fna",
#     output:
#         "{genome}.0123",
#         "{genome}.amb",
#         "{genome}.ann",
#         "{genome}.bwt.2bit.64",
#         "{genome}.pac",
#     log:
#         "logs/bwa-mem2_index/{genome}.log",
#     conda:
#         "../envs/bwa2.yaml",
#     shell:
#         """bwa-mem2 index {input.index}"""

rule align_bowtie2_original:
    input:
        fastq1=lambda wildcards: config["replicate"]["input_folder"] + "{sample}_1.fastq",
        fastq2=lambda wildcards: config["replicate"]["input_folder"] + "{sample}_2.fastq",
    output:
        config["alignment"]["output_folder"] + "bowtie2/{sample}_original.bam"
    log:
        config["alignment"]["output_folder"] + "bowtie2/log/{sample}_original.log"
    conda:
        "../envs/bowtie2.yaml"
    params:
        index=config['alignment']['genome'] + "GRCh38",
    shell:
        """
        {config[alignment][commands][bowtie2]} {params.index} {input.fastq1} {input.fastq2} {output} {log}
        """


rule align_bowtie2_replicates:
    input:
        fastq1=lambda wildcards: gather_checkpoint_outputs(wildcards)[0],
        fastq2=lambda wildcards: gather_checkpoint_outputs(wildcards)[1],
    output:
        config["alignment"]["output_folder"] + "bowtie2/{sample}_sh{n}.bam"\
        if config["replicate"]["replicate_type"] == "sh"\
        else config["alignment"]["output_folder"] + "bowtie2/{sample}_rc.bam"
    log:
        config["alignment"]["output_folder"] + "bowtie2/{sample}_sh{n}.log"\
        if config["replicate"]["replicate_type"] == "sh"\
        else config["alignment"]["output_folder"] + "bowtie2/log/{sample}_rc.log"
    conda:
        "../envs/bowtie2.yaml"
    params:
        index=config['alignment']['genome'] + "GRCh38",
    shell:
        """
        # bowtie2 -x {params.index} -1 {input.fastq1} -2 {input.fastq2} | samtools view -bS -> {output}
        {config[alignment][commands][bowtie2]} {params.index} {input.fastq1} {input.fastq2} {output} {log}
        """
