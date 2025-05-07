rule index_nextgenmap:
    input:
        index="{genome}"
    output:
        "{genome}-enc.2.ngm",
        "{genome}-ht-13-2.3.ngm",
    conda:
        config["alignment"]["envs"]["nextgenmap"],
    shell:
        """
        ngm -r {input.index}
        """


rule align_nextgenmap_original:
    input:
        index=multiext(config['alignment']['genome'], "-enc.2.ngm", "-ht-13-2.3.ngm"),
        fastq1=(lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_1.fastq"
        if config['replicate']['pair_type'] == 'paired'
        else config["replicate"]["input_folder"] + f"{wildcards.sample}.fastq"),
        fastq2=lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_2.fastq"
        if config['replicate']['pair_type'] == 'paired' else [],
    output:
        config["alignment"]["output_folder"] + "nextgenmap/seed_{seed}/" + "bam/{sample}_{ending}.bam",
    log:
        config["alignment"]["output_folder"] + "nextgenmap/seed_{seed}/" + "log/{sample}_{ending}.log",
    params:
        idx=config['alignment']['genome']
    threads: 1
    conda:
        config["alignment"]["envs"]["nextgenmap"]
    wildcard_constraints:
        ending="o",
    shell:
        """(
        ngm -r {params.idx} -1 {input.fastq1} -2 {input.fastq2} --bam -o {output}
        ) >{log} 2>&1"""

rule align_nextgenmap_replicates:
    input:
        index=multiext(config['alignment']['genome'],"-enc.2.ngm","-ht-13-2.3.ngm"),
        fastq1= lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[0] \
        if config['replicate']['pair_type'] == 'paired' else \
        gather_checkpoint_outputs_single(wildcards),
        fastq2=lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[1] \
            if config['replicate']['pair_type'] == 'paired' else [],
    output: config["alignment"]["output_folder"] + "nextgenmap/seed_{seed}/" + "bam/{sample}_{ending}.bam",
    log:
        config["alignment"]["output_folder"] + "nextgenmap/seed_{seed}/" + "log/{sample}_{ending}.log"
    params:
        idx=config['alignment']['genome']
    threads: 1
    conda:
        "../envs/nextgenmap.yaml"
    wildcard_constraints:
        ending="(sh\\d+|both\\d+|rc)"
    shell:
        """(
        ngm -r {params.idx} -1 {input.fastq1} -2 {input.fastq2} --bam -o {output}
        ) >{log} 2>&1"""



