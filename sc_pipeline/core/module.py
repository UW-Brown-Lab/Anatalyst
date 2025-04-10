# sc_pipeline/core/module.py
import os

class AnalysisModule:
    """Base class for all analysis modules in the pipeline"""

    def __init__(self, name, params=None):
        self.name = name
        self.params = params or {}
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

    def save_figure(self, module_name, fig, figsize, name, output_dir, dpi=300):
        """ Save a matplotlib figure to the report images directory

        Args:
            module_name (str): Name of the module generating the figure.
            fig (Matplotlib figure object): The figure to be saved.
            figsize (tuple): Tuple of (width, height) in inches. If none, uses standardized size.
            name (str): Name for the figure file (without extension)
            output_dir (str): Base output directory for the pipeline
            dpi (int, optional): DPI for saving the figure. Defaults to 300.

        Returns:
            str: Path to the saved image file
        """

        image_dir = os.path.join(output_dir, 'images')
        os.makedirs(image_dir, exist_ok=True)
        
        image_path = os.path.join(image_dir, f"{module_name}_{name}.png")

        fig.set_size_inches(figsize)
        fig.savefig(image_path, dpi=dpi, bbox_inches='tight')
        
        return image_path