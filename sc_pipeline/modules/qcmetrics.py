# sc_pipeline/modules/qcmetrics.py

import logging
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from sc_pipeline.core.data_context import DataContext
from sc_pipeline.core.module import AnalysisModule

class QCMetrics(AnalysisModule):
    """Module for calculating quality control metrics."""

    PARAMETER_SCHEMA = {
        'mito_pattern': {
            'type': str,
            'default': '^MT-',
            'description': 'Regex pattern to identify mitochondrial genes'
        },
        'ribo_pattern': {
            'type': str,
            'default': '^RP[SL]',
            'description': 'Regex pattern to identify ribosomal genes'
        },
        'create_plots': {
            'type': bool,
            'default': True,
            'description': 'Generate violin plots for key QC metrics'
        },
        'plot_title': {
            'type':str,
            'description': 'Figure title of created plots in exported report'
        }
    }
    
    def __init__(self, name, params):
        super().__init__(name, params)
        self.logger = logging.getLogger(f"Module.{name}")
        
        # Define required inputs and outputs
        self.required_inputs = ["data"]
        self.outputs = ["data"]  # We'll update the same AnnData object
        
    def run(self, data_context: DataContext):
        """
        Calculate quality control metrics.
        
        Args:
            data_context: The shared data context
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            # Get the data
            adata = data_context.get("data")
            
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
            data_context.set("data", adata)
            self.logger.info("QC metrics calculated successfully")

            if self.params.get("create_plots"):
                sc.pl.violin(
                    adata,
                    ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"],
                    jitter=0.4,
                    multi_panel=True,
                    show=False
                )
                vln_plot = plt.gcf()

                fig_title = self.params.get("plot_title", f"{self.name}: Violin Plot")
                desc = "Plots show the distribution of cells by their number of unique genes, number of UMIs, percentage of features that are mitochondrial, and percentage of features that are ribosomal."
                img_path = self.save_figure(data_context, self.name, vln_plot)
                data_context.add_figure(
                    module_name= self.name,
                    title= fig_title,
                    description= desc,
                    image_path= img_path
                )

            return True
            
        except Exception as e:
            self.logger.error(f"Error calculating QC metrics: {e}", exc_info=True)
            return False