FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="41e6aabe95cbe85ddd5b3475353a1f700d5e6376134f34de5e4b2c1da3a5ff1c"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: envs/bowtie2.yaml
#   prefix: /conda-envs/0ac551e85942d25ba49ec96637298f9f
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - bowtie2=2.5.1
#     - samtools=1.19.2
RUN mkdir -p /conda-envs/0ac551e85942d25ba49ec96637298f9f
COPY envs/bowtie2.yaml /conda-envs/0ac551e85942d25ba49ec96637298f9f/environment.yaml

# Conda environment:
#   source: envs/bwa-mem2.yaml
#   prefix: /conda-envs/1879bc59d779e81a5d242823a80cf3f8
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - bwa-mem2=2.2.1
#     - samtools=1.19.2
RUN mkdir -p /conda-envs/1879bc59d779e81a5d242823a80cf3f8
COPY envs/bwa-mem2.yaml /conda-envs/1879bc59d779e81a5d242823a80cf3f8/environment.yaml

# Conda environment:
#   source: envs/comparer.yaml
#   prefix: /conda-envs/7a104b6018bcc67515ce3276815616ea
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - pandas=1.5.1
#     - pip
#     - pip:
#       - psutil==5.9.3
#       - biopython==1.79
#       - csvsort==1.6.1
#       - typer==0.6.1
RUN mkdir -p /conda-envs/7a104b6018bcc67515ce3276815616ea
COPY envs/comparer.yaml /conda-envs/7a104b6018bcc67515ce3276815616ea/environment.yaml

# Conda environment:
#   source: envs/hisat2.yaml
#   prefix: /conda-envs/dbd436bdbc2da694c9ccf90803de00b7
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - hisat2=2.2.1
#     - samtools=1.19.2
RUN mkdir -p /conda-envs/dbd436bdbc2da694c9ccf90803de00b7
COPY envs/hisat2.yaml /conda-envs/dbd436bdbc2da694c9ccf90803de00b7/environment.yaml

# Conda environment:
#   source: envs/minimap2.yaml
#   prefix: /conda-envs/448dd07561a0ddb91d69e1f01857a4e3
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - minimap2=2.24
#     - samtools=1.19.2
RUN mkdir -p /conda-envs/448dd07561a0ddb91d69e1f01857a4e3
COPY envs/minimap2.yaml /conda-envs/448dd07561a0ddb91d69e1f01857a4e3/environment.yaml

# Conda environment:
#   source: envs/nextgenmap.yaml
#   prefix: /conda-envs/375284d87666c58ad1c45dbbeeddbab0
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - nextgenmap=0.5.5
#     - samtools=1.19.2
RUN mkdir -p /conda-envs/375284d87666c58ad1c45dbbeeddbab0
COPY envs/nextgenmap.yaml /conda-envs/375284d87666c58ad1c45dbbeeddbab0/environment.yaml

# Conda environment:
#   source: envs/parser.yaml
#   prefix: /conda-envs/53dc3cd8cee2f06234f0f26108b5271c
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - prettytable=3.4.1
#     - pysam=0.19.1
#     - pip
#     - pip:
#       - csvsort==1.6.1
#       - psutil==5.9.3
#       - typer==0.6.1
RUN mkdir -p /conda-envs/53dc3cd8cee2f06234f0f26108b5271c
COPY envs/parser.yaml /conda-envs/53dc3cd8cee2f06234f0f26108b5271c/environment.yaml

# Conda environment:
#   source: envs/replicate.yaml
#   prefix: /conda-envs/bd349d6f6f0ad160e2c6f8e171d82ffe
#   channels:
#      - conda-forge
#      - bioconda
#      - defaults
#    dependencies:
#      - chardet=4.0.0
#      - numpy=1.23.4
#      - typer=0.6.1
#      - biopython=1.79
#      - pandas=1.5.1
RUN mkdir -p /conda-envs/bd349d6f6f0ad160e2c6f8e171d82ffe
COPY envs/replicate.yaml /conda-envs/bd349d6f6f0ad160e2c6f8e171d82ffe/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/0ac551e85942d25ba49ec96637298f9f --file /conda-envs/0ac551e85942d25ba49ec96637298f9f/environment.yaml && \
    conda env create --prefix /conda-envs/1879bc59d779e81a5d242823a80cf3f8 --file /conda-envs/1879bc59d779e81a5d242823a80cf3f8/environment.yaml && \
    conda env create --prefix /conda-envs/7a104b6018bcc67515ce3276815616ea --file /conda-envs/7a104b6018bcc67515ce3276815616ea/environment.yaml && \
    conda env create --prefix /conda-envs/dbd436bdbc2da694c9ccf90803de00b7 --file /conda-envs/dbd436bdbc2da694c9ccf90803de00b7/environment.yaml && \
    conda env create --prefix /conda-envs/448dd07561a0ddb91d69e1f01857a4e3 --file /conda-envs/448dd07561a0ddb91d69e1f01857a4e3/environment.yaml && \
    conda env create --prefix /conda-envs/375284d87666c58ad1c45dbbeeddbab0 --file /conda-envs/375284d87666c58ad1c45dbbeeddbab0/environment.yaml && \
    conda env create --prefix /conda-envs/53dc3cd8cee2f06234f0f26108b5271c --file /conda-envs/53dc3cd8cee2f06234f0f26108b5271c/environment.yaml && \
    conda env create --prefix /conda-envs/bd349d6f6f0ad160e2c6f8e171d82ffe --file /conda-envs/bd349d6f6f0ad160e2c6f8e171d82ffe/environment.yaml && \
    conda clean --all -y
