checkpoint create_replicates_single:
    input:
        fastq=config["replicate"]["input_folder"] + "{sample}.fastq",
    output:
        out=config["replicate"]["output_folder"] + "seed_{seed}/" + "{sample}_{ending}.fastq"
    params:
        out_folder=config["replicate"]["output_folder"],
        rep_types=config["replicate"]["replicate_types"],
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
            -r {params.rep_types} \
            -n {params.rep_num} \
            --all \
            -p {params.pair_type} \
            -o {params.out_folder} \
            -s {params.seed}
        """