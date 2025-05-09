# sc_pipeline/modules/dimensionalityreduction.py

import logging
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
from sc_pipeline.core.module import AnalysisModule
from sc_pipeline.utils.adata_utils import set_active_layer

class DimensionalityReduction(AnalysisModule):
    """
    Module for performing dimensionality reduction (PCA, UMAP, t-SNE) on normalized data.
    
    This module is designed to work with both standard scanpy pipelines and with the
    output from PearsonResidualsNormalization (which performs PCA as part of its process).
    """

    PARAMETER_SCHEMA = {
        # PCA parameters
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
        
        # Neighborhood parameters
        'n_neighbors': {
            'type': int,
            'default': 15,
            'description': 'Number of neighbors for neighborhood graph'
        },
        
        # Random seed
        'random_state': {
            'type': int,
            'default': 0,
            'description': 'Random seed for reproducibility'
        },
        
        # Which reductions to run
        'compute_pca': {
            'type': bool,
            'default': True,
            'description': 'Whether to compute/recompute PCA'
        },
        'compute_umap': {
            'type': bool,
            'default': True,
            'description': 'Whether to compute UMAP'
        },
        'compute_tsne': {
            'type': bool,
            'default': True,
            'description': 'Whether to compute t-SNE (computationally expensive)'
        },
        
        # Keys for storing results
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
        },
        
        # Recomputation control
        'force_recompute': {
            'type': bool,
            'default': False,
            'description': 'Force recomputation of existing dimensionality reductions'
        },
        
        # Layer selection
        'layer_key': {
            'type': str,
            'default': None,
            'description': 'If provided, use this layer for computations instead of .X'
        },
        
        # Visualization
        'create_plots': {
            'type': bool,
            'default': True,
            'description': 'Generate plots showing the dimensionality reduction results'
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
            random_state = self.params.get('random_state', 0)
            force_recompute = self.params.get('force_recompute', False)
            
            # Which reductions to compute
            compute_pca = self.params.get('compute_pca', True)
            compute_umap = self.params.get('compute_umap', True)
            compute_tsne = self.params.get('compute_tsne', True)
            
            # Get layer key
            layer_key = self.params.get('layer_key', None)
            
            # Get keys for storing results
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
            
            # Set layer if specified
            if layer_key is not None:
                if layer_key in adata.layers:
                    self.logger.info(f"Using layer '{layer_key}' for dimensionality reduction")
                    adata = set_active_layer(adata, layer_key)
                else:
                    self.logger.warning(f"Layer '{layer_key}' not found, using current .X matrix")
            
            # --- STEP 1: PCA COMPUTATION ---
            
            # Get the PCA key without 'X_' prefix for scanpy compatibility
            pca_key_short = pca_key.replace('X_', '')
            
            # Check if we need to compute PCA
            pca_exists = pca_key in adata.obsm
            pearson_pca_exists = 'X_pca' in adata.obsm and 'pearson_residuals_normalization' in adata.uns
            
            if compute_pca and (force_recompute or not (pca_exists or pearson_pca_exists)):
                # We need to compute PCA
                
                # Check if we have highly variable genes
                if use_highly_variable and 'highly_variable' not in adata.var:
                    self.logger.info("Computing highly variable genes for PCA")
                    self._compute_highly_variable_genes(adata)
                
                # Run the PCA
                self.logger.info(f"Computing PCA with n_comps={n_pcs}, storing as '{pca_key}'")
                sc.pp.pca(
                    adata, 
                    n_comps=n_pcs, 
                    use_highly_variable=use_highly_variable,
                    random_state=random_state,
                    copy=False,
                    key=pca_key_short
                )
                
                # Compute variance explained
                self._compute_variance_explained(adata, pca_key_short)
                
            elif pearson_pca_exists and pca_key != 'X_pca':
                # We can reuse Pearson residuals PCA
                self.logger.info(f"Using existing Pearson residuals PCA, copying to '{pca_key}'")
                
                # Copy the PCA results
                adata.obsm[pca_key] = adata.obsm['X_pca'].copy()
                
                # Copy the PCA information if possible
                if 'pca' in adata.uns:
                    adata.uns[pca_key_short] = adata.uns['pca'].copy()
                    
                    # Ensure variance explained is computed
                    if 'cumulative_variance' not in adata.uns[pca_key_short]:
                        self._compute_variance_explained(adata, pca_key_short)
                else:
                    self.logger.warning("Could not copy PCA variance from Pearson residuals")
                    
            elif pca_exists:
                # PCA already exists under our specified key
                self.logger.info(f"Using existing PCA from '{pca_key}'")
                
                # Ensure variance explained is computed if not already
                if pca_key_short in adata.uns and 'variance_ratio' in adata.uns[pca_key_short]:
                    if 'cumulative_variance' not in adata.uns[pca_key_short]:
                        self._compute_variance_explained(adata, pca_key_short)
                else:
                    self.logger.warning(f"PCA exists in obsm but no variance information found in adata.uns['{pca_key_short}']")
            
            # --- STEP 2: NEIGHBORHOOD COMPUTATION ---
            
            # Check if we need to compute neighborhoods for UMAP/tSNE
            neighbors_needed = compute_umap or compute_tsne
            neighbors_exist = neighbors_key in adata.uns
            
            if neighbors_needed and (force_recompute or not neighbors_exist):
                self.logger.info(f"Computing neighborhood graph with n_neighbors={n_neighbors}, storing as '{neighbors_key}'")
                sc.pp.neighbors(
                    adata, 
                    n_neighbors=n_neighbors, 
                    n_pcs=n_pcs, 
                    random_state=random_state,
                    use_rep=pca_key,
                    key=neighbors_key
                )
            elif neighbors_needed and neighbors_exist:
                self.logger.info(f"Using existing neighborhood graph from '{neighbors_key}'")
            
            # --- STEP 3: UMAP COMPUTATION ---
            
            if compute_umap:
                umap_exists = umap_key in adata.obsm
                
                if force_recompute or not umap_exists:
                    umap_key_short = umap_key.replace('X_', '')
                    self.logger.info(f"Computing UMAP embedding, storing as '{umap_key}'")
                    sc.tl.umap(
                        adata, 
                        random_state=random_state, 
                        neighbors_key=neighbors_key,
                        copy=False,
                        key=umap_key_short
                    )
                else:
                    self.logger.info(f"Using existing UMAP embedding from '{umap_key}'")
            
            # --- STEP 4: t-SNE COMPUTATION ---
            
            if compute_tsne:
                tsne_exists = tsne_key in adata.obsm
                
                if force_recompute or not tsne_exists:
                    tsne_key_short = tsne_key.replace('X_', '')
                    self.logger.info(f"Computing t-SNE embedding, storing as '{tsne_key}'")
                    sc.tl.tsne(
                        adata, 
                        n_pcs=n_pcs, 
                        random_state=random_state,
                        neighbors_key=neighbors_key,
                        copy=False,
                        key=tsne_key_short
                    )
                else:
                    self.logger.info(f"Using existing t-SNE embedding from '{tsne_key}'")
            
            # --- STEP 5: CREATE VISUALIZATIONS ---
            
            if self.params.get('create_plots', True):
                self._create_plots(adata, data_context)
            
            # Update the data context
            data_context.set("data", adata)
            self.logger.info("Dimensionality reduction completed successfully")
            return True
            
        except Exception as e:
            self.logger.error(f"Error in dimensionality reduction: {e}", exc_info=True)
            return False
    
    def _compute_highly_variable_genes(self, adata):
        """Safely compute highly variable genes, handling potential data issues."""
        try:
            # Check for problematic values in the data
            has_problematic_values = False
            
            if sp.issparse(adata.X):
                # For sparse matrix
                if np.any(~np.isfinite(adata.X.data)):
                    has_problematic_values = True
            else:
                # For dense array
                if np.any(~np.isfinite(adata.X)):
                    has_problematic_values = True
            
            # Use appropriate flavor based on data quality
            if has_problematic_values:
                self.logger.warning("Data contains infinity or NaN values. Using 'seurat' flavor for highly variable genes.")
                sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
            else:
                sc.pp.highly_variable_genes(adata, n_top_genes=2000)
                
        except Exception as e:
            self.logger.error(f"Error computing highly variable genes: {e}")
            # Fall back to a safer method
            self.logger.info("Falling back to 'seurat' flavor for highly variable genes")
            sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
            
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
                if any(cumulative_variance >= threshold):
                    n_pcs = np.argmax(cumulative_variance >= threshold) + 1
                    self.logger.info(f"{n_pcs} PCs explain {threshold:.0%} of variance")
                else:
                    self.logger.info(f"All {len(cumulative_variance)} PCs explain less than {threshold:.0%} of variance")
    
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
                
                # Calculate number of PCs needed for 90% variance
                if any(cumulative_variance >= 0.9):
                    n_pcs_90 = np.argmax(cumulative_variance >= 0.9) + 1
                    caption = f"{n_pcs_90} PCs required to explain 90% of variance."
                else:
                    caption = f"All {len(cumulative_variance)} PCs explain less than 90% of variance."
                
                data_context.add_figure(
                    module_name=self.name,
                    title="PCA Variance Explained",
                    description="Left: Variance explained by each principal component. Right: Cumulative variance explained.",
                    image_path=img_path,
                    caption=caption
                )
            
            # Plot PCA projection (first two components)
            if pca_key in adata.obsm:
                fig = plt.figure(figsize=(8, 6))
                sc.pl.embedding(
                    adata, 
                    basis=pca_key, 
                    color=None, 
                    show=False, 
                    title=f"PCA Projection"
                )
                img_path = self.save_figure(data_context, self.name, fig, name="pca_projection")
                
                data_context.add_figure(
                    module_name=self.name,
                    title="PCA Projection",
                    description="Projection of cells onto the first two principal components.",
                    image_path=img_path
                )
            
            # Plot UMAP if computed
            if umap_key in adata.obsm:
                fig = plt.figure(figsize=(8, 6))
                sc.pl.embedding(
                    adata, 
                    basis=umap_key, 
                    color=None, 
                    show=False, 
                    title="UMAP Projection"
                )
                img_path = self.save_figure(data_context, self.name, fig, name="umap_projection")
                
                data_context.add_figure(
                    module_name=self.name,
                    title="UMAP Projection",
                    description="UMAP embedding of cells showing global structure of the dataset.",
                    image_path=img_path
                )
            
            # Plot t-SNE if computed
            if tsne_key in adata.obsm:
                fig = plt.figure(figsize=(8, 6))
                sc.pl.embedding(
                    adata, 
                    basis=tsne_key, 
                    color=None, 
                    show=False, 
                    title="t-SNE Projection"
                )
                img_path = self.save_figure(data_context, self.name, fig, name="tsne_projection")
                
                data_context.add_figure(
                    module_name=self.name,
                    title="t-SNE Projection",
                    description="t-SNE embedding of cells showing global structure of the dataset.",
                    image_path=img_path
                )
                
        except Exception as e:
            self.logger.warning(f"Error creating dimensionality reduction plots: {e}")