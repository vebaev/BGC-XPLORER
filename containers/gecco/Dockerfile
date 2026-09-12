FROM mambaorg/micromamba:2.0.5

RUN micromamba install -y -n base -c conda-forge -c bioconda gecco \
    && micromamba clean -a -y

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["gecco", "--help"]
