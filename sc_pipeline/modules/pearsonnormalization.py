# sc_pipeline/modules/pearsonnormalization.py

import logging
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
from sc_pipeline.core.module import AnalysisModule
from sc_pipeline.utils.adata_utils import set_active_layer, save_layer

class PearsonNormalization(AnalysisModule):
    """
    Module for normalizing data using Pearson residuals method.
    
    This module uses Scanpy's recipe_pearson_residuals method which combines
    highly variable gene selection, Pearson residuals normalization, and PCA
    in a single step.
    """

    PARAMETER_SCHEMA = {
        'theta': {
            'type': float,
            'default': 100,
            'description': 'Overdispersion parameter for the negative binomial model'
        },
        'clip': {
            'type': float,
            'default': None,
            'description': 'Clipping threshold for residuals. If None, uses sqrt(n_cells)'
        },
        'n_top_genes': {
            'type': int,
            'default': 2000,
            'description': 'Number of highly variable genes to select'
        },
        'batch_key': {
            'type': str,
            'default': None,
            'description': 'If provided, performs batch-aware HVG selection using this key in adata.obs'
        },
        'n_comps': {
            'type': int,
            'default': 50,
            'description': 'Number of principal components to compute'
        },
        'random_state': {
            'type': int,
            'default': 0,
            'description': 'Random seed for PCA'
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
            theta = self.params.get('theta', 100)
            clip = self.params.get('clip', None)
            n_top_genes = self.params.get('n_top_genes', 2000)
            batch_key = self.params.get('batch_key', None)
            n_comps = self.params.get('n_comps', 50)
            random_state = self.params.get('random_state', 0)
            layer_key = self.params.get('layer_key', None)
            
            # Log the parameters
            self.logger.info(f"Computing Pearson residuals with theta={theta}, clip={clip} , n_top_genes={n_top_genes}, n_comps={n_comps}")

            if layer_key and layer_key in adata.layers:
                adata = set_active_layer(adata, layer_key) 
            
            # Apply Pearson residuals normalization using the recipe approach
            sc.experimental.pp.recipe_pearson_residuals(
                adata,
                theta=theta,
                clip=clip,
                n_top_genes=n_top_genes,
                batch_key=batch_key,
                n_comps=n_comps,
                random_state=random_state,
                check_values=True,
                inplace=True
            )
            
            # Add residuals as a layer
            pearson_data = adata.uns['pearson_residuals_normalization']['pearson_residuals_df']
            adata = save_layer(adata, name="pearson_residuals", data=pearson_data, make_active=True) 
            
            # Create visualization if requested
            if self.params.get('create_plots', True):
                self._create_plots(adata, data_context)
            
            # Update the data context
            data_context.set("data", adata)
            
            self.logger.info("Pearson residuals normalization completed successfully")
            return True
            
        except Exception as e:
            self.logger.error(f"Error in Pearson residuals normalization: {e}", exc_info=True)
            return False
            
    def _create_plots(self, adata, data_context):
        """Create plots showing the effects of normalization."""
        try:
            # Plot highly variable genes
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
            
            # Plot PCA variance explained
            if 'pca' in adata.uns and 'variance_ratio' in adata.uns['pca']:
                fig, ax = plt.subplots(1, 2, figsize=(12, 5))
                
                # Variance explained by each PC
                variance_ratio = adata.uns['pca']['variance_ratio']
                n_components = min(20, len(variance_ratio))
                ax[0].bar(range(1, n_components + 1), variance_ratio[:n_components] * 100)
                ax[0].set_xlabel('Principal Component')
                ax[0].set_ylabel('Variance Explained (%)')
                ax[0].set_title('Variance Explained by Each PC')
                
                # Cumulative variance explained
                cumulative_variance = np.cumsum(variance_ratio)
                ax[1].plot(range(1, len(cumulative_variance) + 1), cumulative_variance * 100)
                ax[1].set_xlabel('Number of Principal Components')
                ax[1].set_ylabel('Cumulative Variance Explained (%)')
                ax[1].set_title('Cumulative Variance Explained')
                ax[1].axhline(y=90, color='r', linestyle='--')
                ax[1].text(len(cumulative_variance)/2, 91, '90% Variance', color='r')
                
                plt.tight_layout()
                img_path = self.save_figure(data_context, self.name, fig, name="pca_variance")
                
                # Calculate number of PCs needed for 90% variance
                n_pcs_90 = np.argmax(cumulative_variance >= 0.9) + 1 if any(cumulative_variance >= 0.9) else len(cumulative_variance)
                
                data_context.add_figure(
                    module_name=self.name,
                    title="PCA Variance Explained",
                    description="Left: Variance explained by each principal component. Right: Cumulative variance explained.",
                    image_path=img_path,
                    caption=f"{n_pcs_90} PCs required to explain 90% of variance."
                )
            
            # Plot first two PCA components
            if 'X_pca' in adata.obsm:
                fig = plt.figure(figsize=(8, 6))
                sc.pl.pca(adata, color=None, show=False)
                img_path = self.save_figure(data_context, self.name, fig, name="pca_scatter")
                
                data_context.add_figure(
                    module_name=self.name,
                    title="PCA Projection",
                    description="Projection of cells onto the first two principal components after Pearson residuals normalization.",
                    image_path=img_path
                )
            
        except Exception as e:
            self.logger.warning(f"Error creating normalization plots: {e}")