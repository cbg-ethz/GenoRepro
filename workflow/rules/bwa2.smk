rule bwa_index:
    input:
        index=config['alignment']['genome'] + "GRCh38.fna",
    output:
        "{genome}.0123",
        "{genome}.amb",
        "{genome}.ann",
        "{genome}.bwt.2bit.64",
        "{genome}.pac",
    log:
        "logs/bwa-mem2_index/{genome}.log",
    conda:
        "../envs/bwa2.yaml",
    shell:
        """bwa-mem2 index {input.index"""



rule align_bwa2_original:
    input:
        index=config['alignment']['genome'] + "bwa2/GRCh38.fna",
        fastq1=(lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_{wildcards.pair}1.fastq"
                if config['replicate']['pair_type'] == 'paired'
                else config["replicate"]["input_folder"] + f"{wildcards.sample}_{wildcards.pair}.fastq"),
        fastq2=lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_{wildcards.pair}2.fastq"
                if config['replicate']['pair_type'] == 'paired' else [],
    output:
        config["alignment"]["output_folder"] + "bwa2/{sample}_{pair}.bam"
    log:
        config["alignment"]["output_folder"] + "bwa2/log/{sample}_{pair}.log"
    params:
        pair_type=config["replicate"]["pair_type"]
    conda:
        "../envs/bwa2.yaml"
    wildcard_constraints:
        sample="[A-Za-z0-9]+",
        pair="S|R"
    shell:
        """
        if [[ {params.pair_type} == "single" ]]; then
            echo "Single-end alignment for {wildcards.sample}"
            bwa-mem2 mem {input.index} {input.fastq1} | samtools view -bS -> {output} 2> {log}
        else
            echo "Paired-end alignment for {wildcards.sample}"
            bwa-mem2 mem {input.index} {input.fastq1} {input.fastq2} | samtools view -bS -> {output} 2> {log}
        fi
        """

# check if alignment is really run with single or paired-data, split single and paired results to different folders if needed

rule align_bwa2_replicates:
    input:
        index = config['alignment']['genome'] + "bwa2/GRCh38.fna",
        fastq1=lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[0] \
            if config['replicate']['pair_type'] == 'paired' else \
            gather_checkpoint_outputs_single(wildcards),
        fastq2=lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[1] \
            if config['replicate']['pair_type'] == 'paired' else [],
    output:
        config["alignment"]["output_folder"] + "bwa2/{sample}_{Rtype}{n}_{pair}.bam"\
        if config["replicate"]["replicate_type"] == "sh"\
        else config["alignment"]["output_folder"] + "bwa2/{sample}_{Rtype}_{pair}.bam"
    log:
        config["alignment"]["output_folder"] + "bwa2/{sample}_{Rtype}{n}_{pair}.log"\
        if config["replicate"]["replicate_type"] == "sh"\
        else config["alignment"]["output_folder"] + "bwa2/log/{sample}_{Rtype}_{pair}.log"
    params:
        pair_type=config["replicate"]["pair_type"]
    conda:
        "../envs/bwa2.yaml"
    wildcard_constraints:
        sample="[A-Za-z0-9]+",
        pair="S|R",
        Rtype="sh|rc"
    shell:
         """
        if [[ {params.pair_type} == "single" ]]; then
            echo "Single-end alignment for {wildcards.sample}"
            bwa-mem2 mem {input.index} {input.fastq1} | samtools view -bS -> {output} 2> {log}

        else
            # Paired-end alignment command
            echo "Paired-end alignment for {wildcards.sample}"
            bwa-mem2 mem {input.index} {input.fastq1} {input.fastq2} | samtools view -bS -> {output} 2> {log}
        fi
        """
