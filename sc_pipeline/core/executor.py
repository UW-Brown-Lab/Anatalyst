# sc_pipeline/core/executor.py

import os
import logging
import importlib

from .config import ConfigParser
from .data_context import DataContext

class PipelineExecutor:
    """Executes the pipeline based on the configuration."""

    def __init__(self, config_file):
        self.logger = logging.getLogger("PipelineExecutor")

        # Parse the configuration
        parser = ConfigParser()
        self.config = parser.parse(config_file)

        output_dir = self.config['pipeline']['output_dir']

        # Checkpoint settings
        checkpoint_config = self.config['pipeline'].get('checkpointing',{})
        checkpoint_dir = checkpoint_config.get('checkpoint_dir') or os.path.join(output_dir, 'checkpoints')
        self.modules_to_checkpoint = checkpoint_config.get('modules_to_checkpoint', 'all')
        max_checkpoints = checkpoint_config.get('max_checkpoints', 1)

        # Create the data context
        self.data_context = DataContext(checkpoint_dir, max_checkpoints)
        self._set_global_settings()

        # Module registry
        self.module_registry = {}

    def register_module_type(self, module_type, module_class):
        """Register a module type with its implementing class."""
        self.module_registry[module_type] = module_class

    def _get_module_class(self, module_type):
        """Get the class for a module type."""
        if module_type in self.module_registry:
            return self.module_registry[module_type]
        
        # Try to import the module dynamically
        try:
            module_path = f"sc_pipeline.modules.{module_type.lower()}"
            module = importlib.import_module(module_path)

            # Assume the class name is the same as the module type
            class_name = module_type
            if hasattr(module, class_name):
                return getattr(module, class_name)
            
            raise ValueError(f"Module class {class_name} not found in {module_path}")
        
        except ImportError as e:
            self.logger.error(f"Failed to import module {module_type}: {e}")

    def _set_global_settings(self):
        """Store global pipeline settings in the data context"""
        r_memory_limit_gb = self.config['pipeline'].get('r_memory_limit_gb', 8)
        self.data_context.set('R_MEMORY_LIMIT_GB', r_memory_limit_gb)
        self.logger.info(f"Setting global R memory limit to {r_memory_limit_gb} GB")
        self.data_context.set('OUTPUT_DIR', self.config['pipeline']['output_dir'])
        self.data_context.set('FIGURE_SETTINGS', self.config['pipeline']['figure_defaults'])
        self.data_context.set('CONFIG', self.config)

    def run(self, start_from=None):
        """
        Execute the pipeline from the configuration.

        Args:
            start_from: Optional checkpoint to start from

        Returns:
            bool: True if successful, False otherwise
        """
        
        self.logger.info(f"Starting pipeline: {self.config['pipeline']['name']}")

        # Create output directory if it doesn't exist
        output_dir = self.config['pipeline']['output_dir']
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Resume from checkpoint if specified
        if start_from:
            self.logger.info(f"Resuming from checkpoint: {start_from}")
            if not self.data_context.load_checkpoint(start_from):
                self.logger.error(f"Failed to load checkpoint: {start_from}")
                return False
            
        # Get the module configurations
        module_configs = self.config['modules']

        # If starting from a checkpoint, skip modules up to and including the checkpoint
        start_index = 0
        if start_from:
            for i, module_config in enumerate(module_configs):
                if module_config['name'] == start_from:
                    # Start with the NEXT module after the checkpoint
                    start_index = i + 1
                    self.logger.info(f"Resuming execution from module #{start_index+1} (after '{start_from}')")
                    break
            
            # Check if we found the module
            if start_index == 0 and start_from is not None:
                self.logger.error(f"Checkpoint module '{start_from}' not found in configuration")
                return False
                
            # Check if we're at the end of the pipeline
            if start_index >= len(module_configs):
                self.logger.info(f"All modules have been executed. Pipeline is complete.")
                return True
        
        # Execute each module in sequence
        for i in range(start_index, len(module_configs)):
            module_config = module_configs[i]
            module_name = module_config['name']
            module_type = module_config['type']
            module_params = module_config.get('params', {})

            self.logger.info(f"Executing module: {module_name} ({module_type})")

            try:
                # Get the module class and instantiate it
                module_class = self._get_module_class(module_type)
                module = module_class(module_name, module_params)

                # Validate inputs
                if not module.validate_inputs(self.data_context):
                    missing = [input_name for input_name in module.required_inputs if input_name not in self.data_context]
                    self.logger.error(f"Module {module_name} missing required inputs: {missing}")
                    return False
                
                # Run the module
                success = module.run(self.data_context)

                # Check for checkpointing
                should_checkpoint = (
                    self.modules_to_checkpoint == 'all' or
                    module_name in self.modules_to_checkpoint
                )


                if success:
                    self.logger.info(f"Creating checkpoint after module: {module_name}")
                    self.data_context.save_checkpoint(module_name)
                    if should_checkpoint:
                        self.logger.info(f"Creating checkpoint after module: {module_name}")
                        self.data_context.save_checkpoint(module_name)
                else:
                    self.logger.error(f"Module {module_name} failed.")
                    return False
                
                # Validate outputs
                if not module.validate_outputs(self.data_context):
                    missing = [output_name for output_name in module.outputs if output_name not in self.data_context]
                    self.logger.error(f"Module {module_name} failed to produce outputs: {missing}")
                    return False
            except Exception as e:
                self.logger.error(f"Error in module {module_name}: {e}", exc_info=True)
                return False
            
        self.logger.info("Pipeline completed successfully.")
        return True