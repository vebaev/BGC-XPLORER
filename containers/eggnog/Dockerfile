FROM mambaorg/micromamba:2.0.5

RUN micromamba install -y -n base -c conda-forge -c bioconda python=3.11 diamond pip \
    && micromamba run -n base pip install eggnog-mapper \
    && micromamba clean -a -y

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["emapper.py", "--help"]
