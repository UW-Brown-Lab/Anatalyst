# Anatalyst: The Analysis Catalyst - A Modular Pipeline built for Single-cell RNA-seq Analysis

**Anatalyst** is a flexible, modular Python pipeline designed to facilitate boilerplate analytical workflows that leverage existing Python and R programs. The current scope of module support is focused on Single Cell RNA sequencing, built on top of [Scanpy](https://scanpy.readthedocs.io/), providing a customizable workflow for common single-cell analysis tasks.

## Key Features

- **Modular Architecture**: Each analysis step is encapsulated in a module that can be included or excluded as needed
- **Configurable Pipeline**: Simple YAML configuration to customize pipeline parameters
- **Checkpoint System**: Save and resume pipeline execution from checkpoints
- **R Integration**: Seamless integration of R tools (like SoupX) for specialized analyses
- **Reproducible Analysis**: Detailed reports with visualizations and parameter settings
- **Framework Flexibility**: New modules can be custom-built and inserted into the pipeline to allow for any type of sequential analysis

## Workflow Overview

Anatalyst provides a comprehensive workflow for single-cell analysis:

1. **Data Loading**: Import aligned 10X single-cell data using an .h5 file
2. **Quality Control**: Calculate QC metrics and visualize distributions
3. **Ambient RNA Removal**: Remove background RNA contamination using SoupX
4. **Doublet Detection**: Identify and flag potential cell doublets
5. **Cell Filtering**: Filter out low-quality cells and outliers
6. **Pearson Normalization**: Normalize data using Pearson residuals
7. **Dimensionality Reduction**: PCA, UMAP, and t-SNE for visualization and analysis
8. **Report Generation**: Create comprehensive HTML reports with key figures

## Quick Start

TBD - expecting pull Docker container with pipeline already installed
Mount local directory for input/output and extra module insertion?

```bash
# Install the container
docker pull something-or-other

# Run a pipeline with a configuration file
python -m sc_pipeline.scripts.run_pipeline --config my_config.yaml
# This command will likely change and may just be the way the container is launched? 
# Perhaps these args just get passed to docker compose via ENV variables and we make a new entrypoint.sh file?
```

Check out the [Getting Started](getting-started.md) guide for more detailed instructions.

## Pipeline Diagram

```
insert diagram here
```
