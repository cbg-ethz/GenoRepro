rule compare_bam:
    input:
        first_file=config["assessment"]["parser_out"]
        + "{tool}/{sample}_original_{pair}.csv",
        second_file=config["assessment"]["parser_out"]
        + "{tool}/{sample}_{type}_{pair}.csv",
        info_file=config["replicate"]["output_folder"] + "{sample}_{pair}_info.csv",
    output:
        config["assessment"]["comparer_out"]
        + "{tool}/{sample}_{type}-original_{pair}.csv",
    wildcard_constraints:
        sample="[A-Za-z0-9]+",
        pair="S|R",
    conda:
        "../envs/comparer.yaml"
    shell:
        """
        python {REPROFLOW_BASEDIR}/scripts/comparer.py \
            {input.first_file} {input.second_file} {input.info_file} \
            -o {output}
        #&& touch "../work/comparer_out/done.txt" 
        """
