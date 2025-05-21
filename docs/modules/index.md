## Module Dependencies

Some modules depend on the output of other modules. View the details of each module's indivual page to make sure they are implemented in an order that makes sense!


## Module Interface

All modules inherit from the `AnalysisModule` base class and implement the following key methods:

- `__init__(name, params)`: Initializes the module with a name and parameters
- `run(data_context)`: Executes the module's functionality on the data context
- `validate_inputs(data_context)`: Ensures all required inputs are available
- `validate_outputs(data_context)`: Ensures all expected outputs were created

## Creating Custom Modules

You can create custom modules by subclassing `AnalysisModule`. See the [Creating New Modules](../customization/new-modules.md) guide for detailed instructions.

## Configuring Modules

Modules are configured through the pipeline configuration file using YAML syntax. All modules have a name and type. The name is customizable and designed to make logs and checkpoints clear to the user. The type is case-insensitive name of the module file; for example, `DataLoading` maps to `modules/dataloading.py`. Additionally, Each module has its own set of parameters, detailed in the individual module documentation.

Example configuration for a module:

```yaml
- name: filtering # Name that appears in logs and checkpoints
  type: Filtering # maps to modules/filtering.py
  params: # Parameters accepted by filtering module
    filters:
      n_genes_by_counts:
        type: numeric
        min: 200
        max: 6000
      pct_counts_mt:
        type: numeric
        max: 20
    create_plots: true
```

Visit the individual module pages for detailed documentation on each module's parameters and functionality.