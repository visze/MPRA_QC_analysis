FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="21b372b6591863d2430349b75ad1e0dfa59e93ce9cc7cda90149fa0b807b1dd3"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/default.yml
#   prefix: /conda-envs/baf473917eb0b6d3a5d8e1b86138cb93
#   ---
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - biopython
#     - matplotlib
#     - numpy
#     - pandas
#     - pysam
#     - regex
#     - scikit-learn
#     - scipy
#     - seaborn
RUN mkdir -p /conda-envs/baf473917eb0b6d3a5d8e1b86138cb93
COPY workflow/envs/default.yml /conda-envs/baf473917eb0b6d3a5d8e1b86138cb93/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/baf473917eb0b6d3a5d8e1b86138cb93 --file /conda-envs/baf473917eb0b6d3a5d8e1b86138cb93/environment.yaml && \
    conda clean --all -y
