# sc_pipeline/modules/pearsonnormalization.py

import logging
import scanpy as sc
import matplotlib.pyplot as plt
from sc_pipeline.core.module import AnalysisModule

class PearsonResidualsNormalization(AnalysisModule):
    """Module for normalizing data using Pearson residuals method."""

    PARAMETER_SCHEMA = {
        'n_top_genes': {
            'type': int,
            'default': 2000,
            'description': 'Number of highly variable genes to keep'
        },
        'theta': {
            'type': float,
            'default': 100,
            'description': 'Overdispersion parameter for the negative binomial model'
        },
        'chunk_size': {
            'type': int,
            'default': 1000,
            'description': 'Chunk size for processing large datasets'
        },
        'use_highly_variable': {
            'type': bool,
            'default': True,
            'description': 'Whether to subset to highly variable genes after normalization'
        },
        'layer_key': {
            'type': str,
            'default': None,
            'description': 'If provided, use this layer for normalization instead of .X'
        },
        'create_plots': {
            'type': bool,
            'default': True,
            'description': 'Generate plots showing the effects of normalization'
        }
    }
    
    def __init__(self, name, params):
        super().__init__(name, params)
        self.logger = logging.getLogger(f"Module.{name}")
        
        # Define required inputs and outputs
        self.required_inputs = ["data"]
        self.outputs = ["data"]
        
    def run(self, data_context):
        """
        Apply Pearson residuals normalization to the data.
        
        Args:
            data_context: The shared data context
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            # Get the data
            adata = data_context.get("data")
            
            # Get parameters
            n_top_genes = self.params.get('n_top_genes', 2000)
            theta = self.params.get('theta', 100)
            chunk_size = self.params.get('chunk_size', 1000)
            use_highly_variable = self.params.get('use_highly_variable', True)
            layer_key = self.params.get('layer_key', None)
            
            # Determine which matrix to use
            if layer_key and layer_key in adata.layers:
                self.logger.info(f"Using layer '{layer_key}' for normalization")
                X = adata.layers[layer_key]
            else:
                self.logger.info("Using main matrix for normalization")
                X = adata.X
            
            self.logger.info(f"Computing Pearson residuals with theta={theta}, n_top_genes={n_top_genes}")
            
            # Apply Pearson residuals normalization
            sc.experimental.pp.normalize_pearson_residuals(
                adata,
                theta=theta,
                chunk_size=chunk_size,
                n_top_genes=n_top_genes if use_highly_variable else None,
                layer=layer_key
            )
            
            # Log the results
            self.logger.info(f"Normalization completed, data shape: {adata.shape}")
            if use_highly_variable and 'highly_variable' in adata.var:
                n_hvg = adata.var['highly_variable'].sum()
                self.logger.info(f"Identified {n_hvg} highly variable genes")
            
            # Create visualization if requested
            if self.params.get('create_plots', True):
                self._create_plots(adata, data_context)
            
            # Update the data context
            data_context.set("data", adata)
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error in Pearson residuals normalization: {e}", exc_info=True)
            return False
            
    def _create_plots(self, adata, data_context):
        """Create plots showing the effects of normalization."""
        try:
            # Plot highly variable genes if they exist
            if 'highly_variable' in adata.var:
                fig = plt.figure(figsize=(8, 6))
                sc.pl.highly_variable_genes(adata, show=False)
                img_path = self.save_figure(data_context, self.name, fig, name="highly_variable_genes")
                
                data_context.add_figure(
                    module_name=self.name,
                    title="Highly Variable Genes",
                    description="Dispersion vs. mean expression highlighting highly variable genes selected for downstream analysis.",
                    image_path=img_path,
                    caption=f"Selected {adata.var['highly_variable'].sum()} highly variable genes for downstream analysis."
                )
            
            # Plot the distribution of normalized values
            fig, ax = plt.subplots(figsize=(8, 6))
            sc.pl.highest_expr_genes(adata, n_top=20, show=False, ax=ax)
            img_path = self.save_figure(data_context, self.name, fig, name="top_expressed_genes")
            
            data_context.add_figure(
                module_name=self.name,
                title="Top Expressed Genes After Normalization",
                description="Expression levels of the most highly expressed genes across cells after normalization.",
                image_path=img_path
            )
            
        except Exception as e:
            self.logger.warning(f"Error creating normalization plots: {e}")