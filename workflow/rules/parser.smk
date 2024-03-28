rule parse_original_bam:
    input:
        config["alignment"]["output_folder"] + "{tool}/{sample}_{pair}.bam",
    output:
        config["assessment"]["parser_out"] + "{tool}/{sample}_original_{pair}.csv",
    wildcard_constraints:
        sample="[A-Za-z0-9]+",
        pair="S|R"
    conda:
         "../envs/parser.yaml"
    shell:
         """
        python {REPROFLOW_BASEDIR}/scripts/parser.py {input} -o {output}
        """

rule parse_replicates_bam:
    input:
        config["alignment"]["output_folder"] + "{tool}/{sample}_{type}_{pair}.bam",
    output:
        config["assessment"]["parser_out"] + "{tool}/{sample}_{type}_{pair}.csv",
    wildcard_constraints:
        sample="[A-Za-z0-9]+",
        pair="S|R",
    conda:
         "../envs/parser.yaml"
    shell:
         """
        python {REPROFLOW_BASEDIR}/scripts/parser.py {input} -o {output}
        """
