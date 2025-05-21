# QCMetrics Module

The QCMetrics module calculates quality control metrics for single-cell RNA-seq data, including gene and UMI counts per cell, as well as the percentage of reads mapping to mitochondrial and ribosomal genes.

## Overview

Quality control (QC) is a critical step in single-cell RNA-seq analysis. This module uses Scanpy's `pp.calculate_qc_metrics` function to compute standard QC metrics and, optionally, generates visualizations of the distributions.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mito_pattern` | string | `^MT-` | Regex pattern to identify mitochondrial genes |
| `ribo_pattern` | string | `^RP[SL]` | Regex pattern to identify ribosomal genes |
| `create_plots` | boolean | `True` | Whether to generate plots of QC metrics |
| `plot_title` | string | None | Title for generated plots |
| `layer_key` | string | None | If provided, calculate metrics on this layer instead of `.X` |

## Example Configuration

```yaml
- name: pre_qc
  type: QCMetrics
  params:
    mito_pattern: "^MT-"
    ribo_pattern: "^RP[SL]"
    create_plots: true
    plot_title: "Pre-filtering QC Metrics"
```

## Input/Output

### Inputs

- `data`: An AnnData object containing single-cell data

### Outputs

- Updates the `data` object with:
  - `.obs['n_genes_by_counts']`: Number of genes with positive counts in each cell
  - `.obs['total_counts']`: Total counts per cell
  - `.obs['pct_counts_mt']`: Percentage of counts in mitochondrial genes
  - `.obs['pct_counts_ribo']`: Percentage of counts in ribosomal genes
  - Visualization of these metrics in the report (if `create_plots` is `True`)

## Functionality

The module performs the following operations:

1. Calculates basic QC metrics (genes and UMIs per cell)
2. Identifies mitochondrial genes based on the provided pattern
3. Calculates the percentage of reads mapping to mitochondrial genes per cell
4. Identifies ribosomal genes based on the provided pattern
5. Calculates the percentage of reads mapping to ribosomal genes per cell
6. Optionally generates violin plots showing the distribution of these metrics

## Visualization

When `create_plots` is set to `True`, the module generates violin plots for:

- Number of genes detected per cell
- Total UMI counts per cell
- Percentage of counts from mitochondrial genes
- Percentage of counts from ribosomal genes

These plots are saved as images and included in the final report.

![QC Metrics Violin Plot](../images/qc_metrics_example.png)

## Usage Notes

- This module may be run twice in a pipeline: once before filtering to identify outliers and again after filtering to confirm the effectiveness of filtering steps.
- Different organisms may require adjusting the `mito_pattern` and `ribo_pattern`. For example:
  - Human: `^MT-` (mitochondrial), `^RP[SL]` (ribosomal)
  - Mouse: `^mt-` (mitochondrial), `^Rp[sl]` (ribosomal)
- High mitochondrial percentage often indicates cell stress or apoptosis
- High ribosomal percentage may indicate cell stress or high protein synthesis activity

## Implementation Details

The module uses Scanpy's `pp.calculate_qc_metrics` function, which computes:

- Number of genes detected per cell
- Total counts per cell
- Percent of counts in feature subsets (e.g., mitochondrial genes)

The module's visualization uses Scanpy's `pl.violin` functionality for creating the multi-panel violin plots.

## See Also

- [Filtering Module](filtering.md): For filtering cells based on QC metrics
- [Scanpy QC Documentation](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html): For more details on the underlying implementation

## API Reference

::: sc_pipeline.modules.qcmetrics.QCMetrics
    options:
      show_root_heading: true
      show_source: true