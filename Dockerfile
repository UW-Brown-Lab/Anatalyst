FROM ubuntu:22.04
LABEL maintainer="Wesley Blashka <wblashka@wisc.edu>"
LABEL description="Docker image for UW Brown Lab single cell RNA-seq analysis"
LABEL version="1.0"

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

# ---------- Install System Dependencies ---------- #
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    cmake \
    curl \
    g++ \
    gcc \
    git \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libfftw3-dev \
    libigraph-dev \
    libpng-dev \
    libbz2-dev \
    liblzma-dev \
    libboost-all-dev \
    pkg-config \
    python3.11 \
    python3.11-dev \
    python3.11-venv \
    python3-pip \
    software-properties-common \
    wget \
    zlib1g-dev \
    libzmq3-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgeos-dev \
    libgdal-dev \
    libudunits2-dev \
    apt-transport-https \
    gnupg \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ---------- Add R 4.4 ---------- #
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
 && echo "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" > /etc/apt/sources.list.d/r-cran.list \
 && apt-get update && apt-get install -y --no-install-recommends \
    r-base \
    r-base-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ---------- Python Setup ---------- #
RUN python3.11 -m pip install --no-cache-dir --upgrade pip setuptools wheel

RUN python3.11 -m pip install --no-cache-dir \
    scanpy \
    anndata \
    numpy \
    scipy \
    pandas \
    matplotlib \
    seaborn \
    statsmodels \
    scikit-learn \
    python-igraph \
    leidenalg \
    louvain \
    umap-learn \
    pynndescent \
    h5py \
    tables \
    jupyter \
    jupyterlab \
    ipykernel \
    rpy2 \
    markdown \
    scikit-image \
    mkdocs \
    mkdocs-material

# ---------- R Packages ---------- #
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"
RUN R -e "BiocManager::install(c('BiocVersion', 'BiocStyle'))"
RUN R -e "install.packages(c('remotes','tidyverse','devtools','Matrix','ggplot2','dplyr','readr'), repos='https://cloud.r-project.org/', dependencies=TRUE)"
RUN R -e "BiocManager::install(c('SingleCellExperiment','scater','scran','DropletUtils','GenomicFeatures','DelayedArray'))"
RUN R -e "install.packages('Seurat', repos='https://cloud.r-project.org/', dependencies=TRUE)"
RUN R -e "install.packages('SoupX', repos='https://cloud.r-project.org/', dependencies=TRUE)"
RUN R -e "BiocManager::install(c('glmGamPoi','presto'))"
RUN R -e "options(warn = 2); install.packages('IRkernel', repos='https://cloud.r-project.org/'); IRkernel::installspec(user = FALSE)"

# ---------- Defaults ---------- #
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1 \
 && update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1

# Add verification script
RUN echo "#!/bin/bash\nR -e \"installed.packages()[,1]\"" > /usr/local/bin/check-r-packages \
 && chmod +x /usr/local/bin/check-r-packages

# ---------- Workspace configuration (last for caching) ---------- #
WORKDIR /app
COPY config_examples /app/config_examples
COPY sc_pipeline /app/sc_pipeline
COPY scripts /app/scripts
COPY tests /app/tests

ENTRYPOINT ["/bin/bash"]