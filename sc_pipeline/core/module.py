# sc_pipeline/core/module.py
import os

class AnalysisModule:
    """Base class for all analysis modules in the pipeline"""

    # Define empty parameter schema as default
    PARAMETER_SCHEMA: dict[str, dict[str, any]] = {}

    def __init__(self, name, params=None):
        self.name = name
        self.params = self._validate_params(params or {})
        self.required_inputs = []
        self.outputs = []

    def validate_inputs(self, data_context):
        """Check if all required inputs are available in the data context."""
        for input_name in self.required_inputs:
            if input_name not in data_context:
                return False
        return True
        
    def validate_outputs(self, data_context):
        """Check if all expected outputs were created."""
        for output_name in self.outputs:
            if output_name not in data_context:
                return False
        return True
    
    def _validate_params(self, provided_params):
        """Validate parameters against schema and apply defaults."""
        # If no schema is defined, just return the provided parameters ¯\_(ツ)_/¯
        if not hasattr(self, 'PARAMETER_SCHEMA'):
            return provided_params
        
        validated_params = {}

        # Work through parameters, apply default values where applicable, and validate against defined types
        for param_name, param_spec in self.PARAMETER_SCHEMA.items():
            if param_name in provided_params:
                if param_name in provided_params:
                    param_value = provided_params[param_name]

                    if 'type' in param_spec and param_value is not None:
                        expected_type = param_spec['type']

                        # Handle List parameter type
                        if expected_type == list and 'element_type' in param_spec:
                            if not isinstance(param_value, list):
                                self.logger.warning(
                                    f"Parameter '{param_name}' should be a list, got {type(param_value).__name__}"
                                )
                            else:
                                element_type = param_spec['element_type']
                                for i, elem in enumerate(param_value):
                                    if not isinstance(elem, element_type):
                                        self.logger.warning(
                                        f"Element {i} of parameter '{param_name}' should be {element_type.__name__}, "
                                        f"got {type(elem).__name__}"
                                    )

                        # Handle everything else
                        elif not isinstance(param_value, expected_type):
                            self.logger.warning(
                                f"Parameter '{param_name}' should be {expected_type.__name__}, "
                                f"got {type(param_value).__name__}"
                            )
                    
                    validated_params[param_name] = param_value

                # Apply default if available and no param value was provided
                elif 'default' in param_spec:
                    validated_params[param_name] = param_spec['default']
                
                elif param_spec.get('required', False):
                    self.logger.error(f"Required parameter '{param_name}' not provided for module '{self.name}'")

            # Check for other unknown parameters being provided
            unknown_params = set(provided_params.keys()) - set(self.PARAMETER_SCHEMA.keys())
            if unknown_params:
                self.logger.warning(f"Unknown parameters provided to module '{self.name}': {unknown_params}")
            
                # Include unknown parameters in the final set (allows for extensibility or something I guess)
                for param in unknown_params:
                    validated_params[param] = provided_params[param]

            return validated_params

    def run(self, data_context):
        """
        Execute the module's analysis.

        Args:
            data_context: The shared data context containing all pipeline data.

        Returns:
            bool: True if successful, False otherwise
        """
        raise NotImplementedError("Each module must implement its own run method.")
    
    def get_metadata(self):
        """Return metadata about this module."""
        return {
            "name": self.name,
            "params": self.params,
            "required_inputs": self.required_inputs,
            "outputs": self.outputs
        }

    def save_figure(self, module_name, fig, figsize=None, name=None, output_dir=None, dpi=300):
        """ Save a matplotlib figure to the report images directory

        Args:
            module_name (str): Name of the module generating the figure.
            fig (Matplotlib figure object): The figure to be saved.
            figsize (tuple): Tuple of (width, height) in inches. If none, uses standardized size.
            name (str): Name for the figure file (without extension)
            output_dir (str): Base output directory for the pipeline. If none, uses pipeline output directory
            dpi (int, optional): DPI for saving the figure. Defaults to 300.

        Returns:
            str: Path to the saved image file
        """
        if output_dir is None:
            from sc_pipeline.core.executor import current_pipeline
            if hasattr(current_pipeline, 'data_context'):
                output_dir = current_pipeline.data_context.get('OUTPUT_DIR', os.getcwd())
            else:
                output_dir = os.getcwd()

        if figsize is None:
            from sc_pipeline.core.executor import current_pipeline
            if hasattr(current_pipeline, 'data_context'):
                fig_settings = current_pipeline.data_context.get('FIGURE_SETTINGS')
                figsize = fig_settings['width'], fig_settings['height']
            else:
                self.logger.warn(f"Could not retrieve pipeline data context. Proceeding with 8in x 6in...")
                figsize = (8,6)

        image_dir = os.path.join(output_dir, 'images')
        os.makedirs(image_dir, exist_ok=True)
        
        # Generate a sequential number for this module
        if not hasattr(self, '_figure_counter'):
            self._figure_counter = 0
        self._figure_counter += 1

        if name is None:
            name = self._figure_counter

        image_path = os.path.join(image_dir, f"{module_name}_{name}.png")

        fig.set_size_inches(figsize)
        fig.savefig(image_path, dpi=dpi, bbox_inches='tight')
        return image_path