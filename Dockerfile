FROM mambaorg/micromamba:2.0.5

USER root
ENV MAMBA_ROOT_PREFIX=/opt/conda

ARG SNAKEMAKE_VERSION=9.23.1
ARG NICEGUI_VERSION=3.14.0
ARG ANTISMASH_VERSION=8.0.4
ARG GECCO_VERSION=0.10.3
ARG DEEPBGC_VERSION=0.1.31
ARG EGGNOG_MAPPER_VERSION=2.1.13
ARG DBCAN_VERSION=5.2.9
ARG ARTS_COMMIT=8922f296b2a532ba51f4d5daa6a807838c21be24
LABEL org.opencontainers.image.title="BGC-XPLORER" \
      org.opencontainers.image.description="Reproducible biosynthetic gene-cluster discovery workflow" \
      org.opencontainers.image.source="https://github.com/vebaev/bgc-xplorer"

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential default-jdk-headless tar gzip wget curl \
        hmmer ncbi-blast+ diamond-aligner \
    && rm -rf /var/lib/apt/lists/*

# HMMER 2.3.2 is required by antiSMASH for some RiPP modules (hmmpfam2).
RUN cd /tmp \
    && wget http://eddylab.org/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz \
    && tar xzf hmmer-2.3.2.tar.gz \
    && cd hmmer-2.3.2 \
    && ./configure --prefix=/usr/local/hmmer2 \
    && make \
    && make install \
    && ln -sf /usr/local/hmmer2/bin/hmmpfam /usr/local/bin/hmmpfam2 \
    && rm -rf /tmp/hmmer-2.3.2 /tmp/hmmer-2.3.2.tar.gz

# Prodigal is required by DeepBGC; use the upstream static binary since it is
# not available in the base Debian repositories.
RUN cd /tmp \
    && wget https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux -O prodigal \
    && chmod +x prodigal \
    && mv prodigal /usr/local/bin/prodigal

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 snakemake=${SNAKEMAKE_VERSION} pandas pyyaml pip uvicorn \
    && micromamba run -n base pip install nicegui==${NICEGUI_VERSION} \
    && micromamba clean -a -y

RUN micromamba create -y -n antismash -c bioconda -c conda-forge antismash=${ANTISMASH_VERSION} \
    && micromamba clean -a -y

# Ensure antiSMASH uses the externally mounted database directory.
ENV ANTISMASH_DB_DIR=/db/antismash

RUN micromamba create -y -n gecco -c bioconda -c conda-forge gecco=${GECCO_VERSION} \
    && micromamba clean -a -y

RUN micromamba create -y -n deepbgc -c conda-forge -c bioconda \
        python=3.7 hmmer prodigal pip \
    && micromamba run -n deepbgc pip install deepbgc==${DEEPBGC_VERSION} \
    && micromamba clean -a -y

RUN micromamba create -y -n eggnog -c conda-forge -c bioconda \
        python=3.11 eggnog-mapper=${EGGNOG_MAPPER_VERSION} diamond mmseqs2 \
    && ln -sf /opt/conda/envs/eggnog/bin/diamond \
              /opt/conda/envs/eggnog/lib/python3.11/site-packages/eggnogmapper/bin/diamond \
    && ln -sf /opt/conda/envs/eggnog/bin/mmseqs \
              /opt/conda/envs/eggnog/lib/python3.11/site-packages/eggnogmapper/bin/mmseqs \
    && micromamba clean -a -y

RUN micromamba create -y -n dbcan -c conda-forge -c bioconda \
        python=3.10 dbcan=${DBCAN_VERSION} hmmer \
    && micromamba clean -a -y

RUN mkdir -p /opt/arts \
    && wget -qO /tmp/arts.tar.gz "https://github.com/ZiemertLab/ARTS/archive/${ARTS_COMMIT}.tar.gz" \
    && tar xzf /tmp/arts.tar.gz --strip-components=1 -C /opt/arts \
    && rm /tmp/arts.tar.gz

RUN micromamba create -y -n arts -c conda-forge -c bioconda \
        python=3.8 hmmer blast mafft fasttree prodigal pip \
    && micromamba run -n arts pip install -r /opt/arts/requirements.txt \
    && micromamba clean -a -y

RUN cd /opt/arts \
    && tar xzf linux_64bins.tar.gz -C /usr/local/bin/ \
    && chmod +x /usr/local/bin/ranger-dtl-U \
                /usr/local/bin/raxmlHPC-SSE3 \
                /opt/arts/artspipeline1.py

COPY Snakefile /app/Snakefile
COPY rules/ /app/rules/
COPY scripts/ /app/scripts/
RUN /opt/conda/bin/python /app/scripts/startup_checks.py
COPY config/ /app/config/
COPY envs/ /app/envs/
COPY db/manifest.yaml /app/db/manifest.yaml

COPY scripts/entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

WORKDIR /work
EXPOSE 8778 8484

ENTRYPOINT ["/entrypoint.sh"]
