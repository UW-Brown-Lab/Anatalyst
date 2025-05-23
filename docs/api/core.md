# Anatalyst Pipeline Core API

This document provides detailed documentation for the core components of the `Anatalyst` framework. The core modules form the foundation of the pipeline system, handling tasks such as configuration parsing, data management, module execution, and pipeline orchestration.

## Overview

Analyst is designed around a modular architecture that allows for flexible composition of analysis workflows. Each analysis is broken down into a series of modules that each perform a specific task with the data. The core components facilitate this by:

1. **Configuration Management**: Reading and validating YAML configuration files that define pipeline steps
2. **Data Context**: Managing shared data between modules and providing checkpointing capabilities
3. **Module System**: Defining the base interface for all analysis modules
4. **Pipeline Execution**: Orchestrating the execution of modules in sequence

## Core Components

### ConfigParser

The `ConfigParser` class handles the reading and validation of YAML configuration files that define pipeline structure and parameters.

```python
from sc_pipeline.core.config import ConfigParser

parser = ConfigParser()
config = parser.parse("path/to/config.yaml")
```

#### Key Methods

- `parse(config_file)`: Parses a YAML configuration file and returns the validated configuration dictionary.
- `_validate_config(config)`: Internal method that validates the structure of the configuration.
- `_set_defaults(config)`: Internal method that sets default values for optional configuration parameters.

#### Configuration Structure

A basic configuration file should have the following structure:

```yaml
pipeline:
  name: my_pipeline
  output_dir: ./output
  r_memory_limit_gb: 8  # Optional
  figure_defaults:      # Optional
    width: 8
    height: 6
  checkpointing:        # Optional
    enabled: true
    modules_to_checkpoint: all  # or a list of module names
    max_checkpoints: 5
  
modules:
  - name: module1
    type: ModuleType
    params:
      param1: value1
      param2: value2
      
  - name: module2
    type: AnotherModuleType
    params:
      param1: value1
```

#### Scope

In order to provide module level access to the pipeline's configuration, the executor adds the entire parsed configuration to the `DataContext` shortly after its instantiation before any modules are run. This can be accessed from inside any module using `data_context.get('CONFIG')`

### DataContext

The `DataContext` class provides a shared storage space for data passed between modules in the pipeline. It also offers checkpointing capabilities to save and restore pipeline state. Anything handled inside of a module that is not explicitly written to a file space or to the `DataContext` will be unavailable

```python
from sc_pipeline.core.data_context import DataContext

# Create a data context with checkpointing enabled
data_context = DataContext(checkpoint_dir="./checkpoints", max_checkpoints=3)

# Store and retrieve data
context.set("key", value)
value = context.get("key")

# Save and load checkpoints
context.save_checkpoint("after_module1")
context.load_checkpoint("after_module1")
```

#### Figure Generation
To better facilitate record keeping and report generation at the end of the workflow, the `DataContext` class has a built in `add_figure` method. This takes the name of the module generating the figure, and optional title, description, caption, and image_path arguments. This is designed to be compatible with the `ReportGenerator` module, which will essentially string together an html document of these elements. The `save_figure` method of the `AnalysisModule` class returns the path of the saved image, which can then be passed right to the `DataContext.add_figure()` method if you want to include the generated figure.
```python
class MyModule(AnalysisModule):
  def run(self, data_context: DataContext):
    title = "My Figure"
    desc = "A description of my figure"
    caption = "A caption for my figure"
    
    fig = sc.pl.embedding( # Returns matplot figure (from scanpy)
              adata, 
              basis=pca_key, 
              color=None, 
              return_fig=True, 
              title=f"My Figure"
            )
    img_path = self.save_figure(data_context, self.name, fig)
    data_context.add_figure(
        module_name= self.name,
        title= title,
        description= desc,
        caption=caption,
        image_path= img_path
    )
    return
```

#### Key Methods

- `set(key, value)`: Store data by key
- `get(key, default=None)`: Retrieve data by key, with an optional default value
- `__contains__(key)`: Check if a key exists (`key in context`)
- `keys()`: Get all available data keys
- `save_checkpoint(checkpoint_name)`: Save the current state to a checkpoint file (.pkl file of the DataContext itself)
- `load_checkpoint(checkpoint_name)`: Load data from a checkpoint file
- `add_figure(module_name, title=None, description=None, image_path=None, caption=None)`: Add a figure to be included in reports

### AnalysisModule

The `AnalysisModule` class is the base class for all analysis modules in the pipeline. It defines the interface that modules must implement and provides common functionality.

```python
from sc_pipeline.core.module import AnalysisModule

class MyModule(AnalysisModule):
    """Custom module implementation."""
    
    PARAMETER_SCHEMA = {
        'param1': {
            'type': str,
            'required': True,
            'description': 'Description of parameter 1'
        },
        'param2': {
            'type': int,
            'default': 10,
            'description': 'Description of parameter 2'
        }
    }
    
    def __init__(self, name, params):
        super().__init__(name, params)
        self.required_inputs = ["input1", "input2"]
        self.outputs = ["output1"]
        
    def run(self, data_context):
        # Implementation of module functionality
        # Access parameters with self.params.get('param_name', default_value)
        # Access inputs with data_context.get('input_name')
        # Store outputs with data_context.set('output_name', value)
        return True  # Return True if successful, False otherwise
```

#### Key Methods

- `__init__(name, params)`: Initialize the module with a name and parameters
- `run(data_context)`: Execute the module's analysis (must be implemented by subclasses)
- `validate_inputs(data_context)`: Check if all required inputs are available
- `validate_outputs(data_context)`: Check if all expected outputs were created
- `_validate_params(provided_params)`: Validate parameters against schema and apply defaults
- `get_metadata()`: Return metadata about the module
- `save_figure(data_context, module_name, fig, figsize=None, name=None, output_dir=None, dpi=300)`: Save a matplotlib figure

#### Parameter Schema Definition

Module parameters can be defined using the `PARAMETER_SCHEMA` class variable. This schema is used to validate parameters and apply defaults. Each parameter entry should include:

- `type`: The expected Python type (e.g., `str`, `int`, `float`, `bool`, `list`, `dict`)
- `required`: Whether the parameter is required (default: `False`)
- `default`: Default value if not provided (optional)
- `description`: Human-readable description of the parameter (optional)

For list parameters, you can specify the expected element type using `element_type`.

### PipelineExecutor

The `PipelineExecutor` class handles the orchestration of pipeline execution, running modules in sequence according to the configuration.

```python
from sc_pipeline.core.executor import PipelineExecutor

# Create and run a pipeline from a configuration file
executor = PipelineExecutor("path/to/config.yaml")
success = executor.run()

# Optionally, start from a checkpoint
success = executor.run(start_from="after_module1")
```

#### Key Methods

- `__init__(config_file)`: Initialize the executor with a configuration file
- `register_module_type(module_type, module_class)`: Register a module type with its implementing class
- `run(start_from=None)`: Execute the pipeline, optionally starting from a checkpoint
- `_get_module_class(module_type)`: Internal method to get the class for a module type
- `_set_global_settings()`: Internal method to store global pipeline settings

## Utility Modules

### RBridge

The `RBridge` class provides a bridge for calling R functions from Python, using a temporary workspace directory for file exchange.

```python
from sc_pipeline.utils.r_bridge import RBridge

# Create an R bridge with a specified memory limit
r_bridge = RBridge(r_script_dir="./r_scripts", memory_limit_gb=16)

# Run an R script with arguments
success, stdout, stderr = r_bridge.run_r_script("script.R", {"arg1": "value1", "arg2": "value2"})

# Get the path to a file in the workspace
file_path = r_bridge.get_workspace_path("output.txt")

# Clean up when done
r_bridge.cleanup_workspace()
```

### AnnData Utilities

The `adata_utils` module provides utilities for working with AnnData objects, particularly focused on managing layers.

```python
from sc_pipeline.utils.adata_utils import save_layer, set_active_layer

# Save data as a layer in AnnData
adata = save_layer(adata, name="raw_counts", data=None, make_active=False)

# Set a specific layer as the active layer (X matrix)
adata = set_active_layer(adata, layer_name="normalized")
```

## Pipeline Execution Flow

1. The `PipelineExecutor` parses the configuration file using `ConfigParser`
2. Global settings are stored in the `DataContext`
3. Each module is executed in sequence:
   - The module class is dynamically imported based on the module type
   - An instance of the module is created with the specified name and parameters
   - Required inputs are validated
   - The module's `run` method is called with the data context
   - Outputs are validated
   - If checkpointing is enabled, a checkpoint is created
4. If any module fails, execution stops and returns `False`
5. If all modules complete successfully, execution returns `True`

## Checkpointing and Error Recovery

The pipeline supports checkpointing to save and restore state between runs. This allows for recovery from failures without having to restart the entire pipeline.

To enable checkpointing, configure the checkpointing section in your configuration:

```yaml
pipeline:
  # ... other settings ...
  checkpointing:
    enabled: true
    modules_to_checkpoint: all  # or a list of module names
    max_checkpoints: 5 # leaving this at 1 will always save only the last successfully run module's checkpoint
```

To resume from a checkpoint, use the `start_from` parameter when running the pipeline:

```python
executor = PipelineExecutor("path/to/config.yaml")
executor.run(start_from="after_module1")
```

## Module Development Guidelines

When developing new modules for the pipeline, follow these guidelines:

1. Inherit from `AnalysisModule`
2. Define `PARAMETER_SCHEMA` to specify parameters and their validation rules
3. Define `required_inputs` and `outputs` in the `__init__` method
4. Always define a logger instance: `self.logger = logging.getLogger(f"Module.{name}")`
5. Implement the `run(data_context)` method to perform the module's functionality
6. Handle exceptions and log appropriate messages
7. Return `True` if successful, `False` otherwise

Example implementation:

```python
import logging
from sc_pipeline.core.module import AnalysisModule

class MyModule(AnalysisModule):
    """A custom module for the pipeline."""
    
    PARAMETER_SCHEMA = {
        'param1': {
            'type': str,
            'required': True,
            'description': 'Description of parameter 1'
        },
        'param2': {
            'type': int,
            'default': 10,
            'description': 'Description of parameter 2'
        }
    }
    
    def __init__(self, name, params):
        super().__init__(name, params)
        self.logger = logging.getLogger(f"Module.{name}")
        self.required_inputs = ["input1"]
        self.outputs = ["output1"]
        
    def run(self, data_context):
        try:
            # Get inputs
            input_data = data_context.get("input1")
            
            # Get parameters
            param1 = self.params.get('param1')
            param2 = self.params.get('param2', 10)
            
            self.logger.info(f"Processing data with parameters: {param1}, {param2}")
            
            # Process data
            output_data = self._process_data(input_data, param1, param2)
            
            # Store outputs
            data_context.set("output1", output_data)
            
            self.logger.info("Processing completed successfully")
            return True
            
        except Exception as e:
            self.logger.error(f"Error in processing: {e}", exc_info=True)
            return False
            
    def _process_data(self, data, param1, param2):
        # Implementation of data processing
        pass
```

## Troubleshooting and Debugging

### Logging

The pipeline uses Python's logging module to provide detailed logs for troubleshooting. Adjust the log level to see more or less detail:

```python
import logging

# Set log level for detailed output
logging.basicConfig(level=logging.DEBUG)

# For less detailed output
logging.basicConfig(level=logging.INFO)
```

Logger hierarchy:
- `ConfigParser`: Messages related to configuration parsing
- `DataContext`: Messages related to data management and checkpointing
- `PipelineExecutor`: Messages related to pipeline execution
- `Module.<name>`: Messages from specific modules (e.g., `Module.data_loading`)

## Advanced Configuration

### Configuration File Structure

The configuration file controls the entire pipeline. Here's a more detailed explanation of its structure:

```yaml
pipeline:
  name: my_pipeline                # Name of the pipeline
  output_dir: ./output            # Directory for output files
  r_memory_limit_gb: 8            # Memory limit for R scripts
  figure_defaults:                # Default figure settings
    width: 8
    height: 6
  checkpointing:                  # Checkpointing configuration
    enabled: true
    modules_to_checkpoint: all    # or a list of module names
    max_checkpoints: 5
  
modules:                          # List of modules to execute
  - name: module1                 # Name of this module instance
    type: ModuleType              # Type of module (class name)
    params:                       # Module-specific parameters
      param1: value1
      param2: value2
      
  - name: module2
    type: AnotherModuleType
    params:
      param1: value1
```

### Dynamic Module Loading

The pipeline automatically tries to import module classes based on their type. For example, if the configuration specifies `type: DataLoading`, the executor will try to import:

```python
from sc_pipeline.modules.dataloading import DataLoading
```
Currently, the pipeline only supports loading modules that exist within the modules directory. The PipelineExecutor class supports a `register_module_type` method that can be passed a `module_type` parameter used as a key to return the `module_class` parameter. This would allow registering modules outside of the module folder, but it is as of yet unimplemented. For now, custom modules are best placed in the module directory.

### Report Generation

The pipeline includes a `ReportGenerator` module that can generate Markdown and HTML reports from the results of other modules. To add content to the report, modules can add figures using the `data_context.add_figure()` method:

```python
data_context.add_figure(
    module_name=self.name,
    title="Figure Title",
    description="Description of the figure",
    image_path=img_path,
    caption="Figure caption"
)
```
See the DataContext section above for details

## Best Practices

1. **Immutability**: Modules should avoid modifying input data directly, instead creating modified copies or using layers
2. **Error Handling**: Always catch exceptions and log appropriate messages
3. **Validation**: Validate inputs and parameters before processing
4. **Documentation**: Provide detailed docstrings and parameter descriptions
5. **Logging**: Use the logger to provide informative messages at appropriate levels

## Conclusion

The Anatalyst core provides a flexible and extensible framework for building analysis pipelines that leverage both R and Python libraries. By understanding these core components, you can effectively develop new modules and customize the pipeline for specific research needs.

For more information on specific modules, see the [Module Documentation](./modules.md).
