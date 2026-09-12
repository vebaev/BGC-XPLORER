FROM mambaorg/micromamba:2.0.5

RUN micromamba install -y -n base -c conda-forge -c bioconda python=3.7 hmmer prodigal pip \
    && micromamba run -n base pip install deepbgc \
    && micromamba clean -a -y

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["deepbgc", "--help"]
