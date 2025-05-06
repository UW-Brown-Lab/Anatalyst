#!/bin/bash
# Script to install R dependencies for documentation build

# Make the script exit on any error
set -e

# Create user R library directory
USER_R_DIR="$HOME/R/library"
mkdir -p "$USER_R_DIR"

# Set R library path
export R_LIBS_USER="$USER_R_DIR"

echo "Installing R packages to $USER_R_DIR"

# Install BiocManager first
Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org/', lib='$USER_R_DIR')"

# Install required Bioconductor dependencies 
Rscript -e "BiocManager::install(c('BiocVersion', 'BiocStyle'), lib='$USER_R_DIR')"

# Install minimal set of packages needed for documentation
Rscript -e "install.packages(c('Matrix', 'jsonlite'), repos='https://cloud.r-project.org/', lib='$USER_R_DIR')"

# Install bare minimum Bioconductor packages needed for documentation
Rscript -e "BiocManager::install(c('DropletUtils', 'SoupX'), lib='$USER_R_DIR')"

# Install minimal version of Seurat (may still take a while)
Rscript -e "install.packages('Seurat', repos='https://cloud.r-project.org/', lib='$USER_R_DIR', dependencies=FALSE)"

echo "R dependencies installation completed successfully"