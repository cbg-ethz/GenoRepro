rule index_hisat2:
    input:
        index="{genome}"
    output:
        multiext(
            "{genome}",
            ".1.ht2",
            ".2.ht2",
            ".3.ht2",
            ".4.ht2",
            ".5.ht2",
            ".6.ht2",
            ".7.ht2",
            ".8.ht2",
        ),
    params:
        base_name=config['alignment']['genome']
    conda:
        "../envs/hisat2.yaml"
    shell:
        """
        hisat2-build {input.index} {params.base_name}
        """

rule align_hisat2_original:
    input:
        index=multiext(config['alignment']['genome'], ".1.ht2",
            ".2.ht2", ".3.ht2", ".4.ht2", ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"),
        fastq1=(lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_1.fastq"
            if config['replicate']['pair_type'] == 'paired'
            else config["replicate"]["input_folder"] + f"{wildcards.sample}.fastq"),
            fastq2=lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_2.fastq"
        if  config['replicate']['pair_type'] == 'paired' else [],
    output:
        config["alignment"]["output_folder"] + "hisat2/seed_{seed}/" + "bam/{sample}_{ending}.bam"
    log:
        config["alignment"]["output_folder"] + "hisat2/seed_{seed}/" + "log/{sample}_{ending}.log",
    params:
        idx=config['alignment']['genome'],
        base_name=config['alignment']['genome']
    threads: 1
    conda:
        "../envs/hisat2.yaml"
    wildcard_constraints:
        ending="o",
    shell:
        """
        hisat2 --no-spliced-alignment -q -x {params.base_name} -1 {input.fastq1} -2 {input.fastq2} | samtools view -bS -> {output} 2> {log}
        """

rule align_hisat2_replicates:
    input:
        index=multiext(config['alignment']['genome'],".1.ht2",
            ".2.ht2",".3.ht2",".4.ht2",".5.ht2",".6.ht2",".7.ht2",".8.ht2"),
        fastq1= lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[0] \
        if config['replicate']['pair_type'] == 'paired' else \
        gather_checkpoint_outputs_single(wildcards),
        fastq2=lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[1] \
            if config['replicate']['pair_type'] == 'paired' else [],

    output:
        config["alignment"]["output_folder"] + "hisat2/seed_{seed}/" + "bam/{sample}_{ending}.bam"
    log:
        config["alignment"]["output_folder"] + "hisat2/seed_{seed}/" + "log/{sample}_{ending}.log",
    params:
        idx=config['alignment']['genome'],
        base_name=config['alignment']['genome']
    threads: 1
    conda:
        "../envs/hisat2.yaml"
    wildcard_constraints:
        ending="(sh\\d+|both\\d+|rc)"
    shell:
        """
        hisat2 --no-spliced-alignment -q -x {params.base_name} -1 {input.fastq1} -2 {input.fastq2} | samtools view -bS -> {output} 2> {log}
        """