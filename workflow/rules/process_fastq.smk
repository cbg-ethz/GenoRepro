rule generate_csv:
    input:
        input1=(lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_{wildcards.pair}1.fastq"
                if config['replicate']['pair_type'] == 'paired'
                else config["replicate"]["input_folder"] + f"{wildcards.sample}_{wildcards.pair}.fastq"),
        input2=lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_{wildcards.pair}2.fastq"
                if config['replicate']['pair_type'] == 'paired' else [],
    output:
        config["replicate"]["output_folder"] + "{sample}_{pair}_info.csv"
    params:
        pair_type=config["replicate"]["pair_type"]
    conda:
        "../envs/process_fastq.yaml"
    shell:
        """
        echo "pair type is {params.pair_type}"
        if [[ {params.pair_type} == "single" ]]; then
            python {REPROFLOW_BASEDIR}/scripts/process_fastq.py \
            --fastq_file1 {input.input1} {params.pair_type} {output}
        else
            python {REPROFLOW_BASEDIR}/scripts/process_fastq.py \
            --fastq_file1 {input.input1} --fastq_file2 {input.input2} {params.pair_type} {output}
        fi
        """
