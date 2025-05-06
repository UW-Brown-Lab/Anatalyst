#!/bin/bash
# Script to install R dependencies for documentation build

# Make the script exit on any error
set -e

# Install BiocManager first
Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"

# Install required Bioconductor dependencies 
Rscript -e "BiocManager::install(c('BiocVersion', 'BiocStyle'))"

# Install CRAN packages
Rscript -e "install.packages(c('remotes', 'tidyverse', 'devtools', 'Matrix', 'ggplot2', 'dplyr', 'readr'), repos='https://cloud.r-project.org/', dependencies=TRUE)"

# Install Bioconductor packages 
Rscript -e "BiocManager::install(c('SingleCellExperiment', 'scater', 'scran', 'DropletUtils', 'GenomicFeatures', 'DelayedArray'))"

# Install Seurat (from CRAN, not Bioconductor)
Rscript -e "install.packages('Seurat', repos='https://cloud.r-project.org/', dependencies=TRUE)"

# Install SoupX
Rscript -e "install.packages('SoupX', repos='https://cloud.r-project.org/', dependencies=TRUE)"

# Install glmGamPoi and presto for speed
Rscript -e "BiocManager::install(c('glmGamPoi','presto'))"

echo "R dependencies installation completed successfully"