
import os
rule build_snap_docker_image:
    output:
        "snap_docker_built.txt"
    params:
        dockerfile="workflow/docker/Dockerfile.snap",
        image_name="snap:v2.0.2"
    shell:
        """
        docker build -t {params.image_name} -f {params.dockerfile} . \
        && echo "Docker image {params.image_name} built successfully" > {output}
        """

rule index_snap:
    input:
        index=config['alignment']['genome'],
    output:
        [
            os.path.join(os.path.abspath(os.path.dirname(config["alignment"]["genome_path"])), "Genome"),
            os.path.join(os.path.abspath(os.path.dirname(config["alignment"]["genome_path"])), "GenomeIndex"),
            os.path.join(os.path.abspath(os.path.dirname(config["alignment"]["genome_path"])), "GenomeIndexHash"),
            os.path.join(os.path.abspath(os.path.dirname(config["alignment"]["genome_path"])), "OverflowTable"),
        ],
    params:
        docker_image="snap:v2.0.2",
        mount_genome_path=os.path.abspath(os.path.dirname(config["alignment"]["genome_path"])),# Directory containing the genome file
        mount_output="snap_index/",# Directory where the index files will be stored
        genome=config["alignment"]["genome_base"]
    threads: 1,
    shell:
        """
        docker run -it \
            -v {params.mount_genome_path}:/snap_index \
            {params.docker_image} bash -c "
        ./snap-aligner index ../snap_index/{params.genome} ../snap_index/
        "
        """

rule align_snap_original:
    input:
        index=[
            os.path.join(os.path.abspath(os.path.dirname(config["alignment"]["genome_path"])), "Genome"),
            os.path.join(os.path.abspath(os.path.dirname(config["alignment"]["genome_path"])), "GenomeIndex"),
            os.path.join(os.path.abspath(os.path.dirname(config["alignment"]["genome_path"])), "GenomeIndexHash"),
            os.path.join(os.path.abspath(os.path.dirname(config["alignment"]["genome_path"])), "OverflowTable"),
        ],
        fastq1=lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_1.fastq",
        fastq2=lambda wildcards: config["replicate"]["input_folder"] + f"{wildcards.sample}_2.fastq",
    output:
        bam=config["alignment"]["output_folder"] + "snap/seed_{seed}/bam/{sample}_{ending}.bam",
    log:
        config["alignment"]["output_folder"] + "snap/seed_{seed}/log/{sample}_{ending}.log",
    params:
        docker_image="snap:v2.0.2",
        mount_genome_path=os.path.abspath(os.path.dirname(config["alignment"]["genome_path"])),
        out_folder=config["alignment"]["output_folder"] + "snap/seed_{seed}/bam/",
    threads: 1,
    wildcard_constraints:
        ending="o",
    shell:
        """
        docker run -it --rm \
            -v {params.mount_genome_path}:/snap_index \
            -v $(realpath {input.fastq1}):/samples/sample1.fastq \
            -v $(realpath {input.fastq2}):/samples/sample2.fastq \
            -v $(realpath {params.out_folder}):/output \
            {params.docker_image} bash -c "
            ./snap-aligner paired /snap_index /samples/sample1.fastq /samples/sample2.fastq \
            -o /output/{wildcards.sample}_{wildcards.ending}.bam
            "
        """



rule align_snap_replicates:
    input:
        index=[
            os.path.join(os.path.abspath(os.path.dirname(config["alignment"]["genome_path"])), "Genome"),
            os.path.join(os.path.abspath(os.path.dirname(config["alignment"]["genome_path"])), "GenomeIndex"),
            os.path.join(os.path.abspath(os.path.dirname(config["alignment"]["genome_path"])), "GenomeIndexHash"),
            os.path.join(os.path.abspath(os.path.dirname(config["alignment"]["genome_path"])), "OverflowTable"),
        ],
        fastq1= lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[0] \
            if config['replicate']['pair_type'] == 'paired' else \
            gather_checkpoint_outputs_single(wildcards),
        fastq2=lambda wildcards: gather_checkpoint_outputs_paired(wildcards)[1] \
            if config['replicate']['pair_type'] == 'paired' else [],
    output:
        config["alignment"]["output_folder"] + "snap/seed_{seed}/" + "bam/{sample}_{ending}.bam",
    log:
        config["alignment"]["output_folder"] + "snap/seed_{seed}/" + "log/{sample}_{ending}.log",
    params:
        docker_image = "snap:v2.0.2",
        mount_genome_path = os.path.abspath(os.path.dirname(config["alignment"]["genome_path"])),
        out_folder = config["alignment"]["output_folder"] + "snap/seed_{seed}/bam/"
    threads: 1,
    wildcard_constraints:
        ending="(sh\\d+|both\\d+|rc)"
    shell:
        """
        docker run -it --rm \
            -v {params.mount_genome_path}:/snap_index \
            -v $(realpath {input.fastq1}):/samples/sample1.fastq \
            -v $(realpath {input.fastq2}):/samples/sample2.fastq \
            -v $(realpath {params.out_folder}):/output \
            {params.docker_image} bash -c "
            ./snap-aligner paired /snap_index /samples/sample1.fastq /samples/sample2.fastq \
            -o /output/{wildcards.sample}_{wildcards.ending}.bam
            "
        """