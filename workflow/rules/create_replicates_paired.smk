checkpoint create_replicates_paired:
    input:
        fastq1=config["replicate"]["input_folder"] + "{sample}_{pair}1.fastq",
        fastq2=config["replicate"]["input_folder"] + "{sample}_{pair}2.fastq"
    output:
        out1=config["replicate"]["output_folder"] + "{sample}_{Rtype}{n}_{pair}1.fastq"
                if config["replicate"]["replicate_type"] == "sh"
                else config["replicate"]["output_folder"] + "{sample}_{Rtype}_{pair}1.fastq",
        out2=config["replicate"]["output_folder"] + "{sample}_{Rtype}{n}_{pair}2.fastq"
                if config["replicate"]["replicate_type"] == "sh"
                else config["replicate"]["output_folder"] + "{sample}_{Rtype}_{pair}2.fastq",
    params:
        out_folder=config["replicate"]["output_folder"],
        rep_type=config["replicate"]["replicate_type"],
        rep_num=config["replicate"]["replicate_number"],
        pair_type=config["replicate"]["pair_type"],
    conda:
        "../envs/replicate.yaml"
    wildcard_constraints:
        sample="[A-Za-z0-9]+",
        pair="R",
        Rtype="sh|rc"
    shell:
        """
        python {REPROFLOW_BASEDIR}/scripts/create_replicates.py \
            -f1 {input.fastq1} \
            -f2 {input.fastq2} \
            -r {params.rep_type} \
            -n {params.rep_num} \
            -p {params.pair_type} \
            -o {params.out_folder}
        """
