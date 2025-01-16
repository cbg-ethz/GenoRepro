checkpoint create_replicates_single:
    input:
        fastq=config["replicate"]["input_folder"] + "{sample}.fastq",
    output:
        out=config["replicate"]["output_folder"] + "seed_{seed}/" + "{sample}_{Rtype}{n}.fastq"
        if config["replicate"]["replicate_type"] == "sh" or config["replicate"]["replicate_type"] == "both"
        else config["replicate"]["output_folder"] + "seed_{seed}/" + "{sample}_{Rtype}.fastq",
    params:
        out_folder=config["replicate"]["output_folder"],
        rep_type=config["replicate"]["replicate_type"],
        rep_num=config["replicate"]["replicate_number"],
        pair_type=config["replicate"]["pair_type"],
        seed=str(config["replicate"]["seed"]),
    conda:
        "../envs/replicate.yaml"
    wildcard_constraints:
        sample = "[A-Za-z0-9]+",
        Rtype = "sh|rc|both"
    shell:
        """
        echo "output1 fastq1: {output.out}"
        
        python {REPROFLOW_BASEDIR}/scripts/create_replicates.py \
            -f1 {input.fastq} \
            -r {params.rep_type} \
            -n {params.rep_num} \
            -p {params.pair_type} \
            -o {params.out_folder} \
            -s {params.seed}
        """