rule compare_bam:
    input:
        first_file=config["assessment"]["parser_out"] +
                   "{tool}/seed_{seed}/" +
                   "csv/{sample}_o.csv",
        second_file=config["assessment"]["parser_out"] +
                    "{tool}/seed_{seed}/" +
                    "csv/{sample}_{ending}.csv",
        fastq_file = config["replicate"]["input_folder"] + "{sample}_1.fastq"
    output:
        csv_file = config["assessment"]["comparer_out"] +
        "{tool}/seed_{seed}" +
        "/csv/{sample}_o_{ending}.csv",
        done_file=config["assessment"]["comparer_out"] +
        "{tool}/seed_{seed}/" + "done/{sample}_o_{ending}.txt"
    conda:
        "../envs/comparer.yaml"
    wildcard_constraints:
        ending="(sh\\d+|both\\d+|rc)"
    shell:
        """
        echo "done_file: {output.done_file}"

        python {SNAKEFILE_DIR}/scripts/comparer.py \
            {input.first_file} {input.second_file} {input.fastq_file} \
            -x all -f pos -f edit_dist -f quality \
            -o {output.csv_file}

        echo "done" > {output.done_file}
        """
