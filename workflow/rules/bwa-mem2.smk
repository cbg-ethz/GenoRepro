rule index_bwa:
    input:
        index="{genome}",
    output:
        "{genome}.0123",
        "{genome}.amb",
        "{genome}.ann",
        "{genome}.bwt.2bit.64",
        "{genome}.pac",
    # log:
    #     "logs/bwa-mem2_index/{genome}.log",sna
    conda:
        "../envs/bwa-mem2.yaml",

    shell:
        """
        echo "input index: {input.index}"

        bwa-mem2 index {input.index}
        """

rule align_bwa2_original:
    input:
        index=multiext(config['alignment']['genome'], ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123"),
        fastq1=(lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_1.fastq"
            if config['replicate']['pair_type'] == 'paired'
            else config["replicate"]["input_folder"] + f"{wildcards.sample}.fastq"),
        fastq2=lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_2.fastq"
            if config['replicate']['pair_type'] == 'paired' else [],
    output:
        config["alignment"]["output_folder"] + "bwa2/seed_{seed}/" + "bam/{sample}_{ending}.bam"
    log:
        config["alignment"]["output_folder"] + "bwa2/seed_{seed}/" + "log/{sample}_{ending}.log"
    params:
        pair_type=config["replicate"]["pair_type"],
        idx=config['alignment']['genome']
    threads: 1
    conda:
        "../envs/bwa-mem2.yaml"
    wildcard_constraints:
        ending="o",
    shell:
        """
        if [[ {params.pair_type} == "single" ]]; then
            echo "Single-end alignment for {wildcards.sample}"
            bwa-mem2 mem {input.index} {input.fastq1} | samtools view -bS -> {output} 2> {log}
        else
            echo "Paired-end alignment for {wildcards.sample} input: {input} output: {output}"
            bwa-mem2 mem {params.idx} {input.fastq1} {input.fastq2} | samtools view -bS -> {output} 2> {log}
            # bwa-mem2.avx mem {params.idx} {input.fastq1} {input.fastq2} > {output}
        fi
        """

# check if alignment is really run with single or paired-data, split single and paired results to different folders if needed

rule align_bwa2_replicates:
    input:
        index=multiext(config['alignment']['genome'], ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123"),
        fastq1=lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[0] \
            if config['replicate']['pair_type'] == 'paired' else \
            gather_checkpoint_outputs_single(wildcards),
        fastq2=lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[1] \
            if config['replicate']['pair_type'] == 'paired' else [],
    output:
        config["alignment"]["output_folder"] + "bwa2/seed_{seed}/" + "bam/{sample}_{ending}.bam",
    log:
        config["alignment"]["output_folder"] + "bwa2/seed_{seed}/" + "log/{sample}_{ending}.log"
    params:
        pair_type=config["replicate"]["pair_type"],
        idx=config['alignment']['genome']
    threads: 1
    conda:
        "../envs/bwa-mem2.yaml"
    wildcard_constraints:
        ending="(sh\\d+|both\\d+|rc)"
    shell:
        """
       if [[ {params.pair_type} == "single" ]]; then
           echo "Single-end alignment for {wildcards.sample}"
           bwa-mem2 mem {input.index} {input.fastq1} | samtools view -bS -> {output} 2> {log}

       else
           # Paired-end alignment command
           echo "Paired-end alignment for {wildcards.sample} {input.index}"
           bwa-mem2 mem {params.idx} {input.fastq1} {input.fastq2} | samtools view -bS -> {output} 2> {log}
           # bwa-mem2.avx mem {params.idx} {input.fastq1} {input.fastq2} > {output}

       fi
       """