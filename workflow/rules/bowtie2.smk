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
#         "../envs/bowtie2.yaml",
#     shell:
#         """bwa-mem2 index {input.index}"""

rule align_bowtie2_original:
    input:
        fastq1=(lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_{wildcards.pair}1.fastq"
                if config['replicate']['pair_type'] == 'paired'
                else config["replicate"]["input_folder"] + f"{wildcards.sample}_{wildcards.pair}.fastq"),
        fastq2=lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_{wildcards.pair}2.fastq"
                if config['replicate']['pair_type'] == 'paired' else [],
    output:
        config["alignment"]["output_folder"] + "bowtie2/{sample}_{pair}.bam"
    log:
        config["alignment"]["output_folder"] + "bowtie2/log/{sample}_{pair}.log"
    params:
        pair_type=config["replicate"]["pair_type"],
        index=config['alignment']['genome'] + "bowtie2/GRCh38",
    conda:
        "../envs/bowtie2.yaml"
    wildcard_constraints:
        sample="[A-Za-z0-9]+",
        pair="S|R"
    shell:
        """
        if [[ {params.pair_type} == "single" ]]; then
            echo "Single-end alignment for {wildcards.sample}"
            bowtie2 -x {params.index} -U {input.fastq1} | samtools view -bS -> {output} 2> {log}
        else
            echo "Paired-end alignment for {wildcards.sample}"
            bowtie2 -x {params.index} -1 {input.fastq1} -2 {input.fastq2} | samtools view -bS -> {output} 2> {log}
        fi
        """


rule align_bowtie2_replicates:
    input:
        fastq1=lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[0] \
            if config['replicate']['pair_type'] == 'paired' else \
            gather_checkpoint_outputs_single(wildcards),
        fastq2=lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[1] \
            if config['replicate']['pair_type'] == 'paired' else [],
    output:
        config["alignment"]["output_folder"] + "bowtie2/{sample}_{Rtype}{n}_{pair}.bam"\
        if config["replicate"]["replicate_type"] == "sh"\
        else config["alignment"]["output_folder"] + "bowtie2/{sample}_{Rtype}_{pair}.bam"
    log:
        config["alignment"]["output_folder"] + "bowtie2/{sample}_{Rtype}{n}_{pair}.log"\
        if config["replicate"]["replicate_type"] == "sh"\
        else config["alignment"]["output_folder"] + "bowtie2/log/{sample}_{Rtype}_{pair}.log"
    params:
        pair_type=config["replicate"]["pair_type"],
        index=config['alignment']['genome'] + "bowtie2/GRCh38",
    wildcard_constraints:
        sample="[A-Za-z0-9]+",
        pair="S|R",
        Rtype="sh|rc"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        if [[ {params.pair_type} == "single" ]]; then
            echo "Single-end alignment for {wildcards.sample}"
            bowtie2 -x {params.index} -U {input.fastq1} | samtools view -bS -> {output} 2> {log}
        else
            echo "Paired-end alignment for {wildcards.sample}"
            bowtie2 -x {params.index} -1 {input.fastq1} -2 {input.fastq2} | samtools view -bS -> {output} 2> {log}
        fi
        """
