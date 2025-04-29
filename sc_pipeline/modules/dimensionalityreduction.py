# sc_pipeline/modules/dimensionalityreduction.py

import logging
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
from sc_pipeline.core.module import AnalysisModule

class DimensionalityReduction(AnalysisModule):
    """Module for performing dimensionality reduction (PCA, UMAP, t-SNE) on normalized data."""

    PARAMETER_SCHEMA = {
        'n_pcs': {
            'type': int,
            'default': 50,
            'description': 'Number of principal components to compute'
        },
        'use_highly_variable': {
            'type': bool,
            'default': True,
            'description': 'Whether to use only highly variable genes for dimensionality reduction'
        },
        'pca_only': {
            'type': bool,
            'default': False,
            'description': 'If True, only perform PCA without UMAP or t-SNE'
        },
        'n_neighbors': {
            'type': int,
            'default': 15,
            'description': 'Number of neighbors for UMAP and t-SNE'
        },
        'random_state': {
            'type': int,
            'default': 42,
            'description': 'Random seed for reproducibility'
        },
        'run_tsne': {
            'type': bool,
            'default': False,
            'description': 'Whether to run t-SNE (computationally expensive)'
        },
        'run_umap': {
            'type': bool,
            'default': True,
            'description': 'Whether to run UMAP'
        },
        'create_plots': {
            'type': bool,
            'default': True,
            'description': 'Generate plots showing the dimensionality reduction results'
        },
        'pca_key': {
            'type': str,
            'default': 'X_pca',
            'description': 'Key to store PCA results in adata.obsm'
        },
        'umap_key': {
            'type': str,
            'default': 'X_umap',
            'description': 'Key to store UMAP results in adata.obsm'
        },
        'tsne_key': {
            'type': str,
            'default': 'X_tsne',
            'description': 'Key to store t-SNE results in adata.obsm'
        },
        'neighbors_key': {
            'type': str,
            'default': 'neighbors',
            'description': 'Key to store neighborhood graph in adata.uns'
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
        Perform dimensionality reduction on the data.
        
        Args:
            data_context: The shared data context
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            # Get the data
            adata = data_context.get("data")
            
            # Get parameters
            n_pcs = self.params.get('n_pcs', 50)
            use_highly_variable = self.params.get('use_highly_variable', True)
            n_neighbors = self.params.get('n_neighbors', 15)
            random_state = self.params.get('random_state', 42)
            pca_only = self.params.get('pca_only', False)
            run_tsne = self.params.get('run_tsne', False)
            run_umap = self.params.get('run_umap', True)
            
            # Get custom keys for storing results
            pca_key = self.params.get('pca_key', 'X_pca')
            umap_key = self.params.get('umap_key', 'X_umap')
            tsne_key = self.params.get('tsne_key', 'X_tsne')
            neighbors_key = self.params.get('neighbors_key', 'neighbors')
            
            # Store key names for reference
            adata.uns['dim_reduction_keys'] = {
                'pca': pca_key,
                'umap': umap_key,
                'tsne': tsne_key,
                'neighbors': neighbors_key
            }
            
            # Extract PCA key without 'X_' prefix for scanpy compatibility
            pca_key_short = pca_key.replace('X_', '')
            
            # Check if we have highly variable genes if requested
            if use_highly_variable and 'highly_variable' not in adata.var:
                self.logger.warning("Highly variable genes not found. Computing highly variable genes.")
                
                # Check for and handle infinity or NaN values
                has_problematic_values = False
                if sp.issparse(adata.X):
                    # For sparse matrix
                    if np.any(~np.isfinite(adata.X.data)):
                        has_problematic_values = True
                else:
                    # For dense array
                    if np.any(~np.isfinite(adata.X)):
                        has_problematic_values = True
                
                if has_problematic_values:
                    self.logger.warning("Data contains infinity or NaN values. Using 'seurat' flavor for highly variable genes.")
                    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
                else:
                    sc.pp.highly_variable_genes(adata, n_top_genes=2000)

            # Run PCA
            self.logger.info(f"Running PCA with n_comps={n_pcs}, storing as '{pca_key}'")
            sc.pp.pca(
                adata, 
                n_comps=n_pcs, 
                use_highly_variable=use_highly_variable, 
                svd_solver='arpack',
                random_state=random_state,
                copy=False,
                key=pca_key_short  # PCA key doesn't include 'X_' prefix
            )
            
            # Compute PCA variance ratio
            self._compute_variance_explained(adata, pca_key_short)
            
            # Compute neighborhood graph if needed for UMAP or t-SNE
            if not pca_only and (run_umap or run_tsne):
                self.logger.info(f"Computing neighborhood graph with n_neighbors={n_neighbors}, storing as '{neighbors_key}'")
                sc.pp.neighbors(
                    adata, 
                    n_neighbors=n_neighbors, 
                    n_pcs=n_pcs, 
                    random_state=random_state,
                    use_rep=pca_key,
                    key=neighbors_key
                )
            
            # Run UMAP if requested
            if not pca_only and run_umap:
                self.logger.info(f"Computing UMAP embedding, storing as '{umap_key}'")
                # Extract key without 'X_' prefix for compatibility
                umap_key_short = umap_key.replace('X_', '')
                sc.tl.umap(
                    adata, 
                    random_state=random_state, 
                    neighbors_key=neighbors_key,
                    copy=False,
                    key=umap_key_short
                )
            
            # Run t-SNE if requested
            if not pca_only and run_tsne:
                self.logger.info(f"Computing t-SNE embedding, storing as '{tsne_key}'")
                # Extract key without 'X_' prefix for compatibility
                tsne_key_short = tsne_key.replace('X_', '')
                sc.tl.tsne(
                    adata, 
                    n_pcs=n_pcs, 
                    random_state=random_state,
                    neighbors_key=neighbors_key,
                    copy=False,
                    key=tsne_key_short
                )
            
            # Create plots if requested
            if self.params.get('create_plots', True):
                self._create_plots(adata, data_context)
            
            # Update the data context
            data_context.set("data", adata)
            
            self.logger.info("Dimensionality reduction completed successfully")
            return True
            
        except Exception as e:
            self.logger.error(f"Error in dimensionality reduction: {e}", exc_info=True)
            return False
            
    def _compute_variance_explained(self, adata, pca_key='pca'):
        """
        Compute and store the variance explained by PCA components.
        
        Args:
            adata: AnnData object
            pca_key: Key for PCA in adata.uns
        """
        if pca_key in adata.uns and 'variance_ratio' in adata.uns[pca_key]:
            # Calculate cumulative variance
            variance_ratio = adata.uns[pca_key]['variance_ratio']
            cumulative_variance = np.cumsum(variance_ratio)
            
            # Store in uns
            adata.uns[pca_key]['cumulative_variance'] = cumulative_variance
            
            # Find number of PCs for variance thresholds
            for threshold in [0.5, 0.75, 0.9, 0.95]:
                n_pcs = np.argmax(cumulative_variance >= threshold) + 1
                self.logger.info(f"{n_pcs} PCs explain {threshold:.0%} of variance")
    
    def _create_plots(self, adata, data_context):
        """
        Create plots showing the dimensionality reduction results.
        
        Args:
            adata: AnnData object
            data_context: The shared data context
        """
        try:
            # Get the keys for the dimensionality reductions
            dim_keys = adata.uns.get('dim_reduction_keys', {})
            pca_key = dim_keys.get('pca', 'X_pca')
            umap_key = dim_keys.get('umap', 'X_umap')
            tsne_key = dim_keys.get('tsne', 'X_tsne')
            
            # Get the PCA key without 'X_' prefix for uns lookup
            pca_key_short = pca_key.replace('X_', '')
            
            # Plot PCA variance explained
            if pca_key_short in adata.uns and 'variance_ratio' in adata.uns[pca_key_short]:
                fig, ax = plt.subplots(1, 2, figsize=(12, 5))
                
                # Variance explained by each PC
                variance_ratio = adata.uns[pca_key_short]['variance_ratio']
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
                
                data_context.add_figure(
                    module_name=self.name,
                    title="PCA Variance Explained",
                    description="Left: Variance explained by each principal component. Right: Cumulative variance explained.",
                    image_path=img_path,
                    caption=f"{np.argmax(cumulative_variance >= 0.9) + 1} PCs required to explain 90% of variance."
                )
            
            # Plot PCA projection (first two components)
            if pca_key in adata.obsm:
                fig = plt.figure(figsize=(8, 6))
                # Use the key directly for plotting
                sc.pl.embedding(
                    adata, 
                    basis=pca_key, 
                    color=None, 
                    show=False, 
                    title=f"PCA ({pca_key})"
                )
                img_path = self.save_figure(data_context, self.name, fig, name=f"{pca_key.lower()}_projection")
                
                data_context.add_figure(
                    module_name=self.name,
                    title=f"PCA Projection ({pca_key})",
                    description="Projection of cells onto the first two principal components.",
                    image_path=img_path
                )
            
            # Plot UMAP if computed
            if umap_key in adata.obsm:
                fig = plt.figure(figsize=(8, 6))
                # Use the key directly for plotting
                sc.pl.embedding(
                    adata, 
                    basis=umap_key, 
                    color=None, 
                    show=False, 
                    title=f"UMAP ({umap_key})"
                )
                img_path = self.save_figure(data_context, self.name, fig, name=f"{umap_key.lower()}_projection")
                
                data_context.add_figure(
                    module_name=self.name,
                    title=f"UMAP Projection ({umap_key})",
                    description="UMAP embedding of cells showing global structure of the dataset.",
                    image_path=img_path
                )
            
            # Plot t-SNE if computed
            if tsne_key in adata.obsm:
                fig = plt.figure(figsize=(8, 6))
                # Use the key directly for plotting
                sc.pl.embedding(
                    adata, 
                    basis=tsne_key, 
                    color=None, 
                    show=False, 
                    title=f"t-SNE ({tsne_key})"
                )
                img_path = self.save_figure(data_context, self.name, fig, name=f"{tsne_key.lower()}_projection")
                
                data_context.add_figure(
                    module_name=self.name,
                    title=f"t-SNE Projection ({tsne_key})",
                    description="t-SNE embedding of cells showing global structure of the dataset.",
                    image_path=img_path
                )
                
        except Exception as e:
            self.logger.warning(f"Error creating dimensionality reduction plots: {e}")