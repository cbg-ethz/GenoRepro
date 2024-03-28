
rule bwa_index:
    input:
        index=config['alignment']['genome'] + "GRCh38.fna"
    output:
        "{genome}.0123",
        "{genome}.amb",
        "{genome}.ann",
        "{genome}.bwt.2bit.64",
        "{genome}.pac"
    log:
        "logs/bwa-mem2_index/{genome}.log"
    conda:
        "../envs/bwa2.yaml"
    shell:
        """
        bwa-mem2 index {input.index}
        """

rule align_bwa2_original:
    input:
        index=config['alignment']['genome'] + "bwa2/GRCh38.fna",
        fastq1=(lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_1.fastq" if config['replicate']['pair_type'] == 'paired' else lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}.fastq"),
        fastq2=(lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_2.fastq" if config['replicate']['pair_type'] == 'paired' else ""),
    output:
        config["alignment"]["output_folder"] + "bwa2/{sample}_original.bam"
    log:
        config["alignment"]["output_folder"] + "bwa2/log/{sample}_original.log"
    conda:
        "../envs/bwa2.yaml"
    shell:
        """
        if [[ config['replicate']['pair_type'] == "paired" ]]; then
            # Single-end alignment command
            {config[alignment][command][bwa]} {input.index} {input.fastq1} > {output} 2> {log}
        else
            # Paired-end alignment command
            {config[alignment][command][bwa]} {input.index} {input.fastq1} {input.fastq2} > {output} 2> {log}
        fi
        """


rule align_bwa2_replicates:
    input:
        index = config['alignment']['genome'] + "bwa2/GRCh38.fna",
        fastq1=lambda wildcards: gather_checkpoint_outputs(wildcards)[0],
        fastq2=(lambda wildcards: gather_checkpoint_outputs(wildcards)[1] if config['replicate']['pair_type'] == 'paired' else "")
    output:
        config["alignment"]["output_folder"] + "bwa2/{sample}_sh{n}.bam"\
        if config["replicate"]["replicate_type"] == "sh"\
        else config["alignment"]["output_folder"] + "bwa2/{sample}_rc.bam"
    log:
        config["alignment"]["output_folder"] + "bwa2/{sample}_sh{n}.log"\
        if config["replicate"]["replicate_type"] == "sh"\
        else config["alignment"]["output_folder"] + "bwa2/log/{sample}_rc.log"
    conda:
        "../envs/bwa2.yaml"
    shell:
        """
        if [[ config['replicate']['pair_type'] == "paired" ]]; then
            # Single-end alignment command
            {config[alignment][command][bwa]} {input.index} {input.fastq1} > {output} 2> {log}
        else
            # Paired-end alignment command
            {config[alignment][command][bwa]} {input.index} {input.fastq1} {input.fastq2} > {output} 2> {log}
        fi
        """

