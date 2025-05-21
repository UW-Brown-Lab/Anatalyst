# Configuration Overview

Anatalyst uses YAML configuration files to define the pipeline structure and module parameters. This flexible approach allows you to customize the analysis workflow without changing code.

## Configuration File Structure

A typical configuration file has two main sections:

1. **Pipeline Settings**: Global settings for the entire pipeline
2. **Modules**: List of modules to execute, with their parameters

```yaml
pipeline:
  name: my_analysis_pipeline
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
  - name: module1
    type: ModuleType1
    params:
      param1: value1
      param2: value2
      
  - name: module2
    type: ModuleType2
    params:
      param1: value1
      param2: value2
```

## Pipeline Settings

The `pipeline` section contains global settings that apply to the entire pipeline:

| Setting | Type | Description |
|---------|------|-------------|
| `name` | string | Name of the pipeline (used in reports) |
| `output_dir` | string | Directory where results will be saved |
| `r_memory_limit_gb` | number | Memory limit for R processes (in GB) |
| `figure_defaults` | object | Default settings for generated figures |
| `checkpointing` | object | Configuration for checkpoint system |

### Figure Defaults

The `figure_defaults` section controls the default size and appearance of generated figures:

```yaml
figure_defaults:
  width: 8    # Width in inches
  height: 6   # Height in inches
```

### Checkpointing

The `checkpointing` section configures the checkpoint system, which allows resuming a pipeline from intermediate points:

```yaml
checkpointing:
  enabled: true                  # Whether checkpointing is enabled
  modules_to_checkpoint: all     # Which modules to create checkpoints for
  max_checkpoints: 5             # Maximum number of checkpoints to keep
```

For `modules_to_checkpoint`, you can specify:
- `all`: Checkpoint after every module (default)
- A list of module names: `['data_loading', 'filtering']`

## Module Configuration

Each entry in the `modules` list defines a module to execute in the pipeline:

```yaml
- name: module_name              # A unique name for this module instance
  type: ModuleType               # The type of module to use (class name)
  params:                        # Module-specific parameters
    param1: value1
    param2: value2
```

### Module Execution Order

Modules are executed in the order they appear in the configuration file. Each module:

1. Receives the data context from previous modules
2. Applies its analysis steps
3. Updates the data context for subsequent modules

### Module Parameters

Each module has its own set of parameters, defined in the `params` section. Refer to the individual [module documentation](../modules/index.md) for details on available parameters.

## Validation

The pipeline configuration is validated when loaded:

1. Required sections (`pipeline` and `modules`) must be present
2. The `pipeline` section must have a `name` and `output_dir`
3. Each module must have a `name` and `type`
4. Module parameters are validated against each module's parameter schema

## Example Configuration

Here's a minimal example configuration:

```yaml
pipeline:
  name: minimal_pipeline
  output_dir: ./output
  
modules:
  - name: data_loading
    type: DataLoading
    params:
      file_path: /path/to/data.h5
      
  - name: qc_metrics
    type: QCMetrics
    params:
      mito_pattern: "^MT-"
      ribo_pattern: "^RP[SL]"
      
  - name: report_generator
    type: ReportGenerator
```

For more complex examples, see the [Examples](examples.md) page.

## Command Line Options

When running the pipeline, you can specify additional options:

```bash
python -m /workspace/scripts/run_pipeline.py --config my_config.yaml --checkpoint module_name --log-file pipeline.log
```

Options:
- `--config`: Path to the configuration file (required)
- `--checkpoint`: Resume from a specific checkpoint
- `--log-file`: Path to save log output
- `--log-level`: Logging level (DEBUG, INFO, WARNING, ERROR)

## Next Steps

- Browse [example configurations](examples.md) to see complete pipeline setups
- Check the [module documentation](../modules/index.md) for details on available modules and their parameters
- Learn how to [create custom modules](../customization/new-modules.md) to extend the pipeline