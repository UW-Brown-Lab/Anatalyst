# DataLoading Module

The DataLoading module is the entry point for most single-cell RNA-seq analysis pipelines, responsible for loading aligned data from 10X Genomics Cell Ranger output files into the pipeline's data context. The module as it stands has very little flexibility for file types and formats, so feel free to augment this as needed!

## Overview

This module uses Scanpy's `read_10x_h5` function to load data from 10X HDF5 files and prepares it for downstream analysis. It handles the initial data setup, including making variable names unique and creating the first data layer for tracking data transformations throughout the pipeline.

!!! note "Multimodal Data"
    Right now, this only supports single modal data. Any other layers (i.e. Antibody Capture) are discarded
    Future implementations will shift from `AnnData` to using `muon` and `MuData` to resolve this

## Parameters

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `file_path` | string | Yes | - | Path to the 10X h5 file (typically `sample_filtered_feature_bc_matrix.h5`) |

## Example Configuration

```yaml
- name: data_loading
  type: DataLoading
  params:
    file_path: /path/to/filtered_feature_bc_matrix.h5
```

## Input/Output

### Inputs

- None (this is typically the first module in a pipeline)

### Outputs

- `data`: An AnnData object containing the loaded single-cell data with:
  - `.X`: The main expression matrix
  - `.var_names`: Gene identifiers (made unique)
  - `.obs_names`: Cell barcodes
  - `.layers['loaded_counts']`: A copy of the original counts data

## Functionality

The module performs the following operations:

1. **File Validation**: Checks that the specified file path exists and is accessible
2. **Data Loading**: Uses Scanpy's `read_10x_h5` function to load the HDF5 file
3. **Variable Name Processing**: Ensures all gene names are unique using Scanpy's `var_names_make_unique()` function
4. **Layer Creation**: Creates a `loaded_counts` layer to preserve the original data before any transformations
5. **Data Context Storage**: Stores the loaded AnnData object in the pipeline's data context for use by subsequent modules

## Supported File Formats

Currently, the DataLoading module supports:

- **10X HDF5 files** (`.h5`): The standard output format from Cell Ranger, typically named `sample_filtered_feature_bc_matrix.h5`

The module expects the standard 10X HDF5 structure containing:
- Expression matrix (genes × cells)
- Gene information (features)
- Cell barcode information

## Usage Notes

### File Path Requirements

- The `file_path` parameter must point to a valid 10X HDF5 file
- Relative paths are supported and resolved relative to the working directory
- The file must be readable by the user running the pipeline (which should not be an issue when working inside the container)

### Data Structure

After loading, the AnnData object contains:
- **Observations (`.obs`)**: Cell-level metadata (initially just cell barcodes)
- **Variables (`.var`)**: Gene-level metadata (gene IDs, symbols, etc.)
- **Expression Matrix (`.X`)**: Raw count data (typically sparse matrix format)
- **Layers**: The `loaded_counts` layer preserves the original data

### Memory Considerations

Large datasets may require significant memory. Consider the following:
- 10X filtered matrices are typically manageable for most systems
- Raw matrices (containing all droplets) may require more memory
- This initially do not require large amounts of memory; however, the pipeline creates copies of data in layers, which increases memory usage, especially if some layers are not sparse

## Error Handling

The module has minimal error handling, but will log a clear report for the source of the issue.

## Integration with Other Modules

The DataLoading module is designed to work seamlessly with other pipeline modules, providing the data for most downstream operations:

- **Quality Control**: The loaded data is immediately ready for QC metric calculation
- **Ambient RNA Removal**: Can work with both filtered and raw count matrices
- **Preprocessing**: All subsequent modules expect the AnnData structure created by this module

## Example Usage

### Basic Usage

```yaml
modules:
  - name: load_data
    type: DataLoading
    params:
      file_path: ./data/filtered_feature_bc_matrix.h5
```

### With Absolute Path

```yaml
modules:
  - name: load_data
    type: DataLoading
    params:
      file_path: /home/user/scrnaseq_data/sample1/outs/filtered_feature_bc_matrix.h5
```

### In a Complete Pipeline

```yaml
pipeline:
  name: full_analysis
  output_dir: ./output

modules:
  - name: data_loading
    type: DataLoading
    params:
      file_path: ./data/filtered_feature_bc_matrix.h5
      
  - name: qc_metrics
    type: QCMetrics
    params:
      mito_pattern: "^MT-"
      
  # Additional modules...
```

## Implementation Details

The module leverages several key components:

- **Scanpy Integration**: Uses `sc.read_10x_h5()` for reliable 10X file parsing
- **Layer Management**: Utilizes the pipeline's `save_layer()` utility function to track data versions
- **Error Logging**: Comprehensive logging of all operations and errors
- **Data Validation**: Ensures loaded data meets pipeline requirements

## Output Verification

After successful execution, you can verify the data was loaded correctly:

```python
# The loaded AnnData object should have:
print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")
print(f"Layers available: {list(adata.layers.keys())}")
print(f"Matrix type: {type(adata.X)}")
```

Expected output:
```
Loaded 8381 cells and 36601 genes
Layers available: ['loaded_counts']
Matrix type: <class 'scipy.sparse._matrix.csr_matrix'>
```

## Troubleshooting

### Common Issues

**File Not Found Error**
```
Error loading data: [Errno 2] No such file or directory: 'path/to/file.h5'
```
- Verify the file path is correct
- Check file permissions
- Ensure the file exists and isn't corrupted

**Memory Errors**
```
MemoryError: Unable to allocate array
```
- Monitor system memory usage
- Consider using a machine with more RAM for large datasets, especially if using more memory-hungry modules after this one

**Invalid HDF5 Format**
```
OSError: Unable to open file (file signature not found)
```
- Verify the file is a valid HDF5 file
- Check if the file download completed successfully
- Re-download the file if necessary

## Performance Considerations

- **Loading Time**: Varies depending on the file size, but generally a few seconds at most
- **Storage**: Creates additional layers that increase memory footprint

## See Also

- [QCMetrics Module](qcmetrics.md): For calculating quality control metrics on loaded data
- [AmbientRNARemoval Module](ambientrnaremoval.md): For removing background contamination (requires both raw and filtered files)
- [Scanpy Documentation](https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_10x_h5.html): For more details on the underlying data loading function
