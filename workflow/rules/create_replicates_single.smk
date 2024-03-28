checkpoint create_replicates_single:
    input:
        fastq=config["replicate"]["input_folder"] + "{sample}_{pair}.fastq"
    output:
        config["replicate"]["output_folder"] + "{sample}_{pair}_sh{n}.fastq"
        if config["replicate"]["replicate_type"] == "sh"
        else config["replicate"]["output_folder"] + "{sample}_{pair}_rc.fastq",
    params:
        out_folder=config["replicate"]["output_folder"],
        rep_type=config["replicate"]["replicate_type"],
        rep_num=config["replicate"]["replicate_number"],
        pair_type=config["replicate"]["pair_type"],
    conda:
        "../envs/replicate.yaml"
    wildcard_constraints:
        sample="[A-Za-z0-9]+",
        pair="S"
    shell:
        """
        echo "single fastq replicates"
        python {REPROFLOW_BASEDIR}/scripts/create_replicates.py \
            -f1 {input.fastq} \
            -r {params.rep_type} \
            -n {params.rep_num} \
            -p {params.pair_type} \
            -o {params.out_folder}
        """
