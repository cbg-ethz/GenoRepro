rule parse_bam:
    input:
        config["alignment"]["output_folder"] + "{tool}/seed_{seed}/" + "bam/{sample}_{ending}.bam"
    output:
        config["assessment"]["parser_out"] + "{tool}/seed_{seed}/" + "csv/{sample}_{ending}.csv"
    conda:
        "../envs/parser.yaml"
    wildcard_constraints:
        ending="(sh\\d+|both\\d+|rc+|o)"
    shell:
        """
        python {REPROFLOW_BASEDIR}/scripts/parser_bam.py {input} -o {output}
        """






