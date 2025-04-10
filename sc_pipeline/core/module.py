# sc_pipeline/core/module.py

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