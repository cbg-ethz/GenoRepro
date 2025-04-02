# Rule to start Lima and set up Conda with SNAP-aligner
rule start_lima_and_setup_conda:
    output:
        "lima_ready.txt"
    params:
        conda_env="snap-env"
    shell:
        """
        # Start Lima if not running
        if ! limactl list | grep -q "default.*Running"; then
            limactl start default
        fi

        # Set up Conda and install SNAP-aligner inside Lima
        limactl shell default -- bash -c "
            source ~/.bashrc
            if ! command -v conda &> /dev/null; then
                echo 'Conda not found, please install Miniconda in Lima manually'
                exit 1
            fi
            if ! conda env list | grep -q {params.conda_env}; then
                conda create -y -n {params.conda_env} snap-aligner
            fi
        "

        # Confirm setup is complete
        echo 'Lima and Conda environment ready' > {output}
        """

# Rule to run SNAP-aligner for genome indexing inside Lima
rule index_snap:
    input:
        "lima_ready.txt",
        genome=config['alignment']['genome']
    output:
        [
            f"{config['alignment']['genome_path']}/Genome",
            f"{config['alignment']['genome_path']}/GenomeIndex",
            f"{config['alignment']['genome_path']}/GenomeIndexHash",
            f"{config['alignment']['genome_path']}/OverflowTable"
        ]
    params:
        conda_env="snap-env",
        mount_genome_path=config["alignment"]["genome_path"],
        genome=config["alignment"]["genome"].split("/")[-1]  # Extract genome filename
    shell:
        """
        limactl shell default -- bash -c "
            source ~/.bashrc
            conda activate {params.conda_env}

            # Ensure macOS directory is mounted
            sudo mkdir -p /mnt/lima-host
            sudo mount -t 9p -o trans=virtio host /mnt/lima-host

            # Run SNAP-aligner indexing inside Lima
            snap-aligner index /mnt/lima-host/{params.mount_genome_path}/{params.genome} /mnt/lima-host/{params.mount_genome_path}
        "
        """

# Optional: Stop Lima after Snakemake completes execution
rule stop_lima:
    input:
        "lima_ready.txt"
    shell:
        "limactl stop default"