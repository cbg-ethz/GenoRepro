rule index_subread:
    input:
        index="{genome}",
    output:
        multiext(
            "{genome}",
            ".00.b.array",
            ".00.b.tab",
            ".files",
            ".log",
            ".reads"
        ),
    conda:
        "../envs/subread.yaml",
    params:
        base_name=config['alignment']['genome_base']
    shell:
        """
        subread-buildindex -o {input.index} {input.index}
        """

rule align_subread_original:
    input:
        index = multiext(config['alignment']['genome'],".00.b.array",".00.b.tab",".files",".log",".reads"),
        fastq1=( lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_1.fastq"
            if config['replicate']['pair_type'] == 'paired'
            else config["replicate"]["input_folder"] + f"{wildcards.sample}.fastq"),
        fastq2=lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_2.fastq"
            if config['replicate']['pair_type'] == 'paired' else [],
    output:
        config["alignment"]["output_folder"] + "subread/seed_{seed}/" + "bam/{sample}_{ending}.bam"
    log:
        config["alignment"]["output_folder"] + "subread/seed_{seed}/" + "log/{sample}_{ending}.log",
    params:
        idx = config['alignment']['genome']
    threads: 1
    conda:
        "../envs/subread.yaml"
    wildcard_constraints:
        ending="o",
    shell:
        """
        subread-align -t 1 -i {params.idx} -r {input.fastq1} -R {input.fastq2} -o {output}
        """

rule align_subread_replicates:
    input:
        index=multiext(config['alignment']['genome'],".00.b.array",".00.b.tab",".files",".log",".reads"),
        fastq1= lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[0] \
            if config['replicate']['pair_type'] == 'paired' else \
            gather_checkpoint_outputs_single(wildcards),
        fastq2=lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[1] \
            if config['replicate']['pair_type'] == 'paired' else [],
    output:
        config["alignment"]["output_folder"] + "subread/seed_{seed}/" + "bam/{sample}_{ending}.bam"
    log:
        config["alignment"]["output_folder"] + "subread/seed_{seed}/" + "log/{sample}_{ending}.log",
    params:
        idx = config['alignment']['genome']
    threads: 1
    conda:
        "../envs/subread.yaml"
    wildcard_constraints:
        ending="(sh\\d+|both\\d+|rc)",
    shell:
        """
        subread-align -t 1 -i {params.idx} -r {input.fastq1} -R {input.fastq2} -o {output}
        """
