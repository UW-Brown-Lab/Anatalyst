# Creating Custom Modules for the Single-Cell Pipeline

## Overview

Anatalyst is designed to be extensible, allowing you to create custom modules that integrate seamlessly with the existing workflow or create an entirely new workflow. This guide will walk you through the process of creating a new analysis module.

## Module Structure

Each module is a Python class that inherits from `AnalysisModule`. The basic structure includes:

- Initialization method
- Parameter schema
- `run` method
- Optional helper methods

## Step-by-Step Guide

### 1. Create the Module File

Create a new file in `sc_pipeline/modules/` with a descriptive name, e.g., `myanalysis.py`.

### 2. Import Required Modules

```python
# At a minimum
import logging
from sc_pipeline.core.module import AnalysisModule
```

### 3. Define the Module Class

```python
class MyAnalysis(AnalysisModule):
    """
    Description of your custom module's purpose.
    """

    # Define the parameter schema
    PARAMETER_SCHEMA = {
        'param1': {
            'type': str,  # Parameter type
            'required': True,  # Whether the parameter is mandatory
            'description': 'Description of the parameter'
        },
        'optional_param': {
            'type': int,
            'default': 10,  # Default value if not provided
            'description': 'An optional parameter'
        }
    }

    def __init__(self, name, params):
        # Call the parent class constructor
        super().__init__(name, params)
        
        # Set up logging
        self.logger = logging.getLogger(f"Module.{name}")
        
        # Define required inputs and outputs
        self.required_inputs = ["data"]  # Inputs this module needs
        self.outputs = ["data"]  # Outputs this module will produce
```

### 4. Implement the `run` Method

```python
def run(self, data_context):
    """
    Main method to execute the module's analysis.
    
    Args:
        data_context: Shared data context containing pipeline data

    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Retrieve the AnnData object
        adata = data_context.get("data")
        
        # Access module parameters
        param1 = self.params.get('param1')
        optional_param = self.params.get('optional_param', 10)
        
        # Perform your analysis
        # Example: Do something with the data
        
        # Optional: Create visualization
        if self.params.get('create_plots', True):
            self._create_plots(adata, data_context)
        
        # Update the data context
        data_context.set("data", adata)
        
        return True
    
    except Exception as e:
        self.logger.error(f"Error in analysis: {e}", exc_info=True)
        return False
```

### 5. Add Optional Visualization Method

```python
def _create_plots(self, adata, data_context):
    """
    Create and save visualization figures.
    
    Args:
        adata: AnnData object
        data_context: Shared data context
    """
    try:
        # Create a matplotlib figure
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # Create your plot
        # sc.pl.something(adata, ax=ax)
        
        # Save the figure using the parent class method
        img_path = self.save_figure(data_context, self.name, fig)
        
        # Add figure to the report
        data_context.add_figure(
            module_name=self.name,
            title="My Analysis Plot",
            description="Description of the plot",
            image_path=img_path
        )
    
    except Exception as e:
        self.logger.warning(f"Error creating plots: {e}")
```

## Module Best Practices

- Use logging for tracking module progress and errors
- Handle exceptions gracefully
- Provide clear documentation and parameter descriptions
- Create visualizations when possible
- Minimize data transformation, preferring to add information to the AnnData object

## Using Your Custom Module

To use your new module in a pipeline configuration:

```yaml
modules:
  - name: my_custom_analysis
    type: MyAnalysis
    params:
      param1: "example_value"
      optional_param: 20
```


## Notes

- Modules should be idempotent (can be run multiple times without side effects)
- Prefer adding information to AnnData's `.obs`, `.var`, `.uns`, or as layers
- Keep modules focused on a single type of analysis
