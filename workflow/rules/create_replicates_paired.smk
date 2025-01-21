checkpoint create_replicates_paired:
    input:
        fastq1=config["replicate"]["input_folder"] + "{sample}_1.fastq",
        fastq2=config["replicate"]["input_folder"] + "{sample}_2.fastq"
    output:
        out1=config["replicate"]["output_folder"] + "seed_{seed}/" + "{sample}_{ending}_1.fastq",
        out2=config["replicate"]["output_folder"] + "seed_{seed}/" + "{sample}_{ending}_2.fastq"
    params:
        out_folder=config["replicate"]["output_folder"],
        rep_types=config["replicate"]["replicate_types"],
        rep_num=config["replicate"]["replicate_number"],
        pair_type=config["replicate"]["pair_type"],
        seed=str(config["replicate"]["seed"]),
    conda:
        "../envs/replicate.yaml"
    wildcard_constraints:
        sample="[A-Za-z0-9]+",
        Rtype="sh|rc|both"
    shell:
        """
        echo "output1 fastq1: {output.out1}"

        python {REPROFLOW_BASEDIR}/scripts/create_replicates.py \
            -f1 {input.fastq1} \
            -f2 {input.fastq2} \
            -r {params.rep_types} \
            -n {params.rep_num} \
            --all \
            -p {params.pair_type} \
            -o {params.out_folder} \
            -s {params.seed}
        """

