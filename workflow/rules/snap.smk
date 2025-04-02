rule index_snap:
    input:
        index=config['alignment']['genome_path']
    output:
        Genome = config['alignment']['genome_path'] + 'Genome',
        GenomeIndex = config['alignment']['genome_path'] + 'GenomeIndex',
        GenomeIndexHash = config['alignment']['genome_path'] + 'GenomeIndexHash',
        OverflowTable = config['alignment']['genome_path'] + 'OverflowTable'
    conda:
        "../envs/snap.yaml"
    params:
        out_folder=config['alignment']['genome_path']
    shell:
        """
        snap-aligner index {input} {params.out_folder}
        """