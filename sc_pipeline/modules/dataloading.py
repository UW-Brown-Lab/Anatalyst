# sc_pipeline/modules/dataloading.py

import os
import logging
import scanpy as sc
from sc_pipeline.core.module import AnalysisModule

class DataLoading(AnalysisModule):
    """
    Module for loading aligned data.
    Currently supports only 10X sample_filtered_feature_bc_matrix.h5
    """

    PARAMETER_SCHEMA = {
        'file_path': {
            'type': str,
            'required': True,
            'description': 'Path to the 10X h5 file'
        }
    }

    def __init__(self, name, params):
        super().__init__(name,params)
        self.logger = logging.getLogger(f"Module.{name}")

        # Define required inputs and outputs
        self.required_inputs = []
        self.outputs = ["data"]
    
    def run(self, data_context):
        """
        Load aligned data outputted by cellranger. Currently only supports filtered_feature_bc_matrix.h5
        Args:
            data_context: The shared data context

        Returns:
            bool: True if successful, False otherwise
        """
        file_path = self.params.get('file_path')

        if not file_path:
            self.logger.error("No file path specified")
            return False
        
        if not os.path.exists(file_path):
            self.logger.error(f"File not found: {file_path}")
            return False
        
        try:
            self.logger.info(f"Loading data from: {file_path}")

            # Load the data
            adata = sc.read_10x_h5(file_path)

            # Make variable names unique
            adata.var_names_make_unique()

            # Store in data context
            data_context.set("data", adata)

            self.logger.info(f"Loaded data with shape: {adata.shape}")
            return True
        
        except Exception as e:
            self.logger.error(f"Error loading data: {e}", exc_info=True)
            return False