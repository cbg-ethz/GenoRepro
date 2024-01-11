rule create_replicates:
    input:
        config["replicates"]["input"],
    output:
        directory(config["replicates"]["output"]),
        #dynamic(directory(config["replicates"]["output"])/"replicates_{sample}.fastq")),

        #"{directory}/reiplicates",
    params:
        SEED=2,
        CORES=1,
    conda:
        "../envs/replicate.yaml",
    shell:
        "python {REPROFLOW_BASEDIR}/scripts/create_replicates.py  --csv_file {input} --output {output} --cores {params.CORES} --seed {params.SEED}"
        
