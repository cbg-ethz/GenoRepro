rule index_minimap2:
    input:
        index="{genome}"
    output:
        "{genome}.mmi"
    conda:
        "../envs/minimap2.yaml",
    threads: 1
    shell:
        """
        minimap2 -x sr -d {output} {input}
        """

rule align_minimap2_original:
    input:
        index=config['alignment']['genome'] + ".mmi",
        fastq1=( lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_1.fastq"
        if config['replicate']['pair_type'] == 'paired'
        else config["replicate"]["input_folder"] + f"{wildcards.sample}.fastq"),
        fastq2=lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_2.fastq"
        if config['replicate']['pair_type'] == 'paired' else [],
    output:
        config["alignment"]["output_folder"] + "minimap2/seed_{seed}/" + "bam/{sample}_{ending}.bam"
    log:
        config["alignment"]["output_folder"] + "bwa2/seed_{seed}/" + "log/{sample}_{ending}.log"
    params:
        idx=config['alignment']['genome']
    conda:
        "../envs/minimap2.yaml"
    wildcard_constraints:
        ending="o",
    threads: 1
    shell:
        """
        minimap2 -ax sr {params.idx} {input.fastq1} {input.fastq2} | samtools view -bS -> {output} 2> {log}
        """


rule align_minimap2_replicates:
    input:
        index = config['alignment']['genome'] + ".mmi",
        fastq1= lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[0] \
            if config['replicate']['pair_type'] == 'paired' else \
            gather_checkpoint_outputs_single(wildcards),
        fastq2=lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[1] \
            if config['replicate']['pair_type'] == 'paired' else [],
    output:
        config["alignment"]["output_folder"] + "minimap2/seed_{seed}/" + "bam/{sample}_{ending}.bam",
    log:
        config["alignment"]["output_folder"] + "minimap2/seed_{seed}/" + "log/{sample}_{ending}.log"
    params:
        idx=config['alignment']['genome']
    threads: 1
    conda:
        "../envs/minimap2.yaml"
    wildcard_constraints:
        ending="(sh\\d+|both\\d+|rc)"
    shell:
        """
        minimap2 -ax sr {params.idx} {input.fastq1} {input.fastq2} | samtools view -bS -> {output} 2> {log}
        """




