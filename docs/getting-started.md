# Getting Started

This guide will walk you through running your first analysis with Anatalyst.

## Basic Usage

Anatalyst is designed to be run from a configuration file that defines the modules to use and their parameters. The basic workflow is:

1. Create a configuration file
2. Run the pipeline
3. Review the results

## Creating a Configuration File

Create a YAML file that defines your pipeline. Here's a minimal example:

```yaml
pipeline:
  name: my_first_pipeline
  output_dir: ./output
  r_memory_limit_gb: 8
  figure_defaults:
    width: 8
    height: 6
  
modules:
  - name: data_loading
    type: DataLoading
    params:
      file_path: /path/to/your/filtered_feature_bc_matrix.h5
      
  - name: qc_metrics
    type: QCMetrics
    params:
      mito_pattern: "^MT-"
      ribo_pattern: "^RP[SL]"
      
  - name: report_generator
    type: ReportGenerator
```

Save this as `minimal_config.yaml`.

## Running the Pipeline

Execute the pipeline using the `run_pipeline.py` script:

```bash
python -m /workspace/scripts/run_pipeline.py --config minimal_config.yaml
```

This will:

1. Load your data
2. Calculate QC metrics
3. Generate a report in the output directory

## Example with a Complete Analysis

Here's a more comprehensive example configuration for a full analysis workflow:

```yaml
pipeline:
  name: complete_analysis
  output_dir: ./output
  r_memory_limit_gb: 8
  figure_defaults:
    width: 8
    height: 6
  checkpointing:
    enabled: true
    modules_to_checkpoint: all
    max_checkpoints: 5
  
modules:
  - name: data_loading
    type: DataLoading
    params:
      file_path: /path/to/filtered_feature_bc_matrix.h5
      
  - name: pre_qc
    type: QCMetrics
    params:
      mito_pattern: "^MT-"
      ribo_pattern: "^RP[SL]"
      create_plots: true
      
  - name: ambient_removal
    type: AmbientRNARemoval
    params:
      raw_counts_path: /path/to/raw_feature_bc_matrix.h5
      filtered_counts_path: /path/to/filtered_feature_bc_matrix.h5
      ndims: 30
      resolution: 0.8
      
  - name: doublet_detection
    type: DoubletDetection
    params:
      expected_doublet_rate: 0.05
      
  - name: filtering
    type: Filtering
    params:
      filters:
        n_genes_by_counts:
          type: numeric
        pct_counts_mt:
          type: numeric
        predicted_doublet:
          type: boolean
      
  - name: post_qc
    type: QCMetrics
    params:
      mito_pattern: "^MT-"
      ribo_pattern: "^RP[SL]"
      
  - name: normalization
    type: PearsonNormalization
    params:
      n_top_genes: 2000
      
  - name: dim_reduction
    type: DimensionalityReduction
    params:
      n_pcs: 50
      compute_umap: true
      compute_tsne: true
      
  - name: report_generator
    type: ReportGenerator
    params:
      generate_html: true
```

Save this as `complete_config.yaml` and run it as before:

```bash
python -m /workspace/scripts/run_pipeline.py --config complete_config.yaml
```

## Resuming from a Checkpoint

If your pipeline fails or stops for any reason, you can resume from the last checkpoint:

```bash
python -m /workspace/scripts/run_pipeline.py --config complete_config.yaml --checkpoint doublet_detection
```

This will resume the pipeline from after the "doublet_detection" module.

## Viewing the Results

After running the pipeline, check the output directory:

```
output/
├── checkpoints/          # Pipeline checkpoints
├── images/               # Generated figures
├── analysis_report.md    # Markdown report
└── analysis_report.html  # HTML report 
```

Open the HTML report in a web browser to view the analysis results, including visualizations and parameter settings.

## Next Steps

- Explore the [Configuration](configuration/index.md) documentation to learn about all available options
- Browse the [Modules](modules/index.md) section to understand each analysis step in detail
- Check out the [Examples](examples/index.md) for more use cases and sample analyses