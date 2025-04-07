rule index_snap:
    input:
        index=config['alignment']['genome_path']
    output:
        Genome=config['alignment']['genome_path'] + 'Genome',
        GenomeIndex=config['alignment']['genome_path'] + 'GenomeIndex',
        GenomeIndexHash=config['alignment']['genome_path'] + 'GenomeIndexHash',
        OverflowTable=config['alignment']['genome_path'] + 'OverflowTable'
    conda:
        "../envs/snap.yaml"
    params:
        out_folder=config['alignment']['genome_path']
    shell:
        """
        snap-aligner index {input} {params.out_folder}
        """

rule align_snap_original:
    input:
        index=[config['alignment']['genome_path'] + 'Genome',
               config['alignment']['genome_path'] + 'GenomeIndex',
               config['alignment']['genome_path'] + 'GenomeIndexHash',
               config['alignment']['genome_path'] + 'OverflowTable',],
        fastq1=( lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_1.fastq"
            if config['replicate']['pair_type'] == 'paired'
            else config["replicate"]["input_folder"] + f"{wildcards.sample}.fastq"),
        fastq2=lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_2.fastq"
            if config['replicate']['pair_type'] == 'paired' else [],
    output:
        config["alignment"]["output_folder"] + "snap/seed_{seed}/" + "bam/{sample}_{ending}.bam"
    log:
        config["alignment"]["output_folder"] + "snap/seed_{seed}/" + "log/{sample}_{ending}.log",
    params:
        idx=config['alignment']['genome']
    threads: 1
    conda:
        "../envs/snap.yaml"
    wildcard_constraints:
        ending="o",
    shell:
        """
        snap-aligner {params.idx} {input.fastq1} {input.fastq2} -o -bam {output}
        """