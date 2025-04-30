FROM ubuntu:22.04
LABEL maintainer="Wesley Blashka <wblashka@wisc.edu>"
LABEL description="Docker image for UW Brown Lab single cell RNA-seq analysis"
LABEL version="1.0"

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

# Set working directory
WORKDIR /app

# Install system dependencies
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
    # Additional dependencies that might help with IRkernel
    libzmq3-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libtiff5-dev \
    libjpeg-dev \
    # Dependencies specifically for Seurat and SoupX
    libgeos-dev \
    libgdal-dev \
    libudunits2-dev \
    apt-transport-https \
    gnupg \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Add R 4.4 repository 
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" > /etc/apt/sources.list.d/r-cran.list

# Install R 4.4
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base \
    r-base-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Verify R version
RUN R --version

# Set up Python environment
RUN python3.11 -m pip install --no-cache-dir --upgrade pip setuptools wheel

# Install Python packages for single cell analysis
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
    markdown

# Install BiocManager first
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"

# Install required Bioconductor dependencies 
RUN R -e "BiocManager::install(c('BiocVersion', 'BiocStyle'))"

# Install CRAN packages
RUN R -e "install.packages(c('remotes', 'tidyverse', 'devtools', 'Matrix', 'ggplot2', 'dplyr', 'readr'), repos='https://cloud.r-project.org/', dependencies=TRUE)"

# Install Bioconductor packages 
RUN R -e "BiocManager::install(c('SingleCellExperiment', 'scater', 'scran', 'DropletUtils', 'GenomicFeatures', 'DelayedArray'))"

# Install Seurat (from CRAN, not Bioconductor)
RUN R -e "install.packages('Seurat', repos='https://cloud.r-project.org/', dependencies=TRUE)"

# Install SoupX
RUN R -e "install.packages('SoupX', repos='https://cloud.r-project.org/', dependencies=TRUE)"

# Install glmGamPoi and presto for speed
RUN R -e "BiocManager::install(c('glmGamPoi','presto'))"

# Make Python 3.11 the default Python
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1

# Install IRkernel as root (before switching users)
RUN R -e "options(warn = 2); install.packages('IRkernel', repos='https://cloud.r-project.org/'); IRkernel::installspec(user = FALSE)"

# Create a non-root user
RUN useradd -m -s /bin/bash researcher

# Create R library directory for the user and fix permissions
RUN mkdir -p /home/researcher/R/library && \
    chown -R researcher:researcher /home/researcher/R && \
    echo "R_LIBS_USER='/home/researcher/R/library'" >> /home/researcher/.Renviron

# Add verification script
RUN echo "#!/bin/bash\nR -e \"installed.packages()[,1]\"" > /usr/local/bin/check-r-packages && \
    chmod +x /usr/local/bin/check-r-packages

# Switch to non-root user
USER researcher

# Set up working directory for the user
WORKDIR /home/researcher/work

# Set up entry point
ENTRYPOINT ["/bin/bash"]