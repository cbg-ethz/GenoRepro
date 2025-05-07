rule index_smalt:
    input:
        index="{genome}",
    output:
        "{genome}.sma",
        "{genome}.smi",
    conda:
        config["alignment"]["envs"]["smalt"],

    shell:
        """
        echo "input index: {input.index}"

        smalt index {input.index} {input.index}
        """

rule align_smalt_original:
    input:
        index = multiext(config['alignment']['genome'],".sma",".smi"),
        fastq1=( lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_1.fastq"
            if config['replicate']['pair_type'] == 'paired'
            else config["replicate"]["input_folder"] + f"{wildcards.sample}.fastq"),
        fastq2=lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_2.fastq"
            if config['replicate']['pair_type'] == 'paired' else [],
    output:
        config["alignment"]["output_folder"] + "smalt/seed_{seed}/" + "bam/{sample}_{ending}.bam"
    log:
        config["alignment"]["output_folder"] + "smalt/seed_{seed}/" + "log/{sample}_{ending}.log",
    params:
        idx = config['alignment']['genome']
    threads: 1
    conda:
        config["alignment"]["envs"]["smalt"]
    wildcard_constraints:
        ending="o",
    shell:
        """
        smalt map -r 1 {params.idx} {input.fastq1} {input.fastq2} | samtools view -bS -> {output} 2> {log}
        """

rule align_smalt_replicates:
    input:
        index = multiext(config['alignment']['genome'],".sma",".smi"),
        fastq1= lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[0] \
            if config['replicate']['pair_type'] == 'paired' else \
            gather_checkpoint_outputs_single(wildcards),
        fastq2=lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[1] \
            if config['replicate']['pair_type'] == 'paired' else [],
    output:
        config["alignment"]["output_folder"] + "smalt/seed_{seed}/" + "bam/{sample}_{ending}.bam"
    log:
        config["alignment"]["output_folder"] + "smalt/seed_{seed}/" + "log/{sample}_{ending}.log",
    params:
        idx = config['alignment']['genome']
    threads: 1
    conda:
        "../envs/smalt.yaml"
    wildcard_constraints:
        ending="(sh\\d+|both\\d+|rc)",
    shell:
        """
        smalt map -r 1 {params.idx} {input.fastq1} {input.fastq2} | samtools view -bS -> {output} 2> {log}
        """
