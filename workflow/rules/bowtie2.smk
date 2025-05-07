rule index_bowtie2:
    input:
        index="{genome}",
    output:
        "{genome}.1.bt2",
        "{genome}.2.bt2",
        "{genome}.3.bt2",
        "{genome}.4.bt2",
        "{genome}.rev.1.bt2",
        "{genome}.rev.2.bt2",

    conda:
        config["alignment"]["envs"]["bowtie2"],
    shell:
        """
        echo "bowtie2 indexing {input.index}"

        bowtie2-build {input.index} {wildcards.genome}
        """


rule align_bowtie2_original:
    input:
        index = multiext(config['alignment']['genome'],".1.bt2",".2.bt2",".3.bt2",".4.bt2",".rev.1.bt2",".rev.2.bt2"),
        fastq1=(lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_1.fastq"
        if config['replicate']['pair_type'] == 'paired'
        else config["replicate"]["input_folder"] + f"{wildcards.sample}.fastq"),
        fastq2=lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_2.fastq"
        if config['replicate']['pair_type'] == 'paired' else [],
    output:
        config["alignment"]["output_folder"] + "bowtie2/seed_{seed}/" + "bam/{sample}_{ending}.bam"
    log:
        config["alignment"]["output_folder"] + "bowtie2/seed_{seed}/" + "log/{sample}_{ending}.log"
    params:
        pair_type=config["replicate"]["pair_type"],
        idx=config['alignment']['genome']
    threads: 6
    conda:
        "../envs/bowtie2.yaml"
    wildcard_constraints:
        ending="o",
    shell:
        """
        if [[ {params.pair_type} == "single" ]]; then
            echo "Single-end alignment for {wildcards.sample}"
            bowtie2 -x {input.index} -U {input.fastq1} | samtools view -bS -> {output} 2> {log}
        else
            echo "Paired-end alignment for {wildcards.sample}"
            bowtie2 -x {params.idx} -1 {input.fastq1} -2 {input.fastq2} | samtools view -bS -> {output} 2> {log}
        fi
        """

rule align_bowtie2_replicates:
    input:
        index = multiext(config['alignment']['genome'],".1.bt2",".2.bt2",".3.bt2",".4.bt2",".rev.1.bt2",".rev.2.bt2"),
        fastq1=lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[0] \
            if config['replicate']['pair_type'] == 'paired' else \
            gather_checkpoint_outputs_single(wildcards),
        fastq2=lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[1] \
            if config['replicate']['pair_type'] == 'paired' else [],

    output:
        config["alignment"]["output_folder"] + "bowtie2/seed_{seed}/" + "bam/{sample}_{ending}.bam",
    log:
        config["alignment"]["output_folder"] + "bowtie2/seed_{seed}/" + "log/{sample}_{ending}.log"
    params:
        pair_type=config["replicate"]["pair_type"],
        idx=config['alignment']['genome']
    threads: 6
    conda:
        config["alignment"]["envs"]["bowtie2"],
    wildcard_constraints:
        ending="(sh\\d+|both\\d+|rc)"
    shell:
        """
        if [[ {params.pair_type} == "single" ]]; then
            echo "Single-end alignment for {wildcards.sample}"
            bowtie2 -x {input.index} -U {input.fastq1} | samtools view -bS -> {output} 2> {log}
        else
            echo "Paired-end alignment for {wildcards.sample}"
            bowtie2 -x {params.idx} -1 {input.fastq1} -2 {input.fastq2} | samtools view -bS -> {output} 2> {log}
        fi
        """



