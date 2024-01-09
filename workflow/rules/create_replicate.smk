rule create_replicate:
    input:
        files="../../work/input.csv",
        script="../scripts/create_replicates.py",
    params:
        SEED=2,
        CORES=1,
    shell:
        "python {input.script} --csv_file {input.files}"
