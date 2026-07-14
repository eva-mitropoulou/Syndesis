FROM mambaorg/micromamba:1.5.10

WORKDIR /workspace/syndesis

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean -a -y

COPY --chown=$MAMBA_USER:$MAMBA_USER . .
RUN pip install -e .

ENV PYTHONUNBUFFERED=1

CMD ["syndesis", "--help"]

