# sc_pipeline/modules/qc_metrics.py

import logging
import numpy as np
import pandas as pd
import scanpy as sc
from sc_pipeline.core.module import AnalysisModule

class QCMetrics(AnalysisModule):
    """Module for calculating quality control metrics."""
    
    def __init__(self, name, params):
        super().__init__(name, params)
        self.logger = logging.getLogger(f"Module.{name}")
        
        # Define required inputs and outputs
        self.required_inputs = ["data"]
        self.outputs = ["data"]  # We'll update the same AnnData object
        
    def run(self, data_context):
        """
        Calculate quality control metrics.
        
        Args:
            data_context: The shared data context
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            # Get the data
            adata = data_context.get("loaded_data")
            
            self.logger.info("Calculating QC metrics")
            
            # Calculate number of genes and UMIs per cell
            sc.pp.calculate_qc_metrics(
                adata, 
                inplace=True, 
                percent_top=None, 
                log1p=False
            )
            
            # Calculate mitochondrial gene percentage
            mito_pattern = self.params.get('mito_pattern', '^MT-')
            if mito_pattern:
                adata.var['mt'] = adata.var_names.str.match(mito_pattern)
                sc.pp.calculate_qc_metrics(
                    adata, 
                    qc_vars=['mt'], 
                    percent_top=None, 
                    log1p=False, 
                    inplace=True
                )
                
            # Calculate ribosomal gene percentage
            ribo_pattern = self.params.get('ribo_pattern', '^RP[SL]')
            if ribo_pattern:
                adata.var['ribo'] = adata.var_names.str.match(ribo_pattern)
                sc.pp.calculate_qc_metrics(
                    adata, 
                    qc_vars=['ribo'], 
                    percent_top=None, 
                    log1p=False, 
                    inplace=True
                )
                
            # Update the data context
            data_context.set("loaded_data", adata)
            
            self.logger.info("QC metrics calculated successfully")
            return True
            
        except Exception as e:
            self.logger.error(f"Error calculating QC metrics: {e}", exc_info=True)
            return False