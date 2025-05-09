# sc_pipeline/modules/doubletdetection.py

import logging
import scanpy as sc
from sc_pipeline.core.module import AnalysisModule
from sc_pipeline.utils.adata_utils import set_active_layer

class DoubletDetection(AnalysisModule):
    """Module for detecting doublets using scanpy's scrublet wrapper."""

    PARAMETER_SCHEMA = {
        'expected_doublet_rate': {
            'type': float,
            'default': 0.05,
            'description': 'The expected doublet rate'
        },
        'min_counts': {
            'type': int,
            'default': 3,
            'description': 'Minimum number of counts for a gene to be used'
        },
        'min_cells': {
            'type': int,
            'default': 3,
            'description': 'Minimum number of cells for a gene to be used'
        },
        'n_neighbors': {
            'type': int,
            'default': 30,
            'description': 'Number of nearest neighbors for predicting doublets'
        },
        'sim_doublet_ratio': {
            'type': float,
            'default': 2.0,
            'description': 'Ratio of simulated doublets to original cells'
        },
        'create_plots': {
            'type': bool,
            'default': True,
            'description': 'Whether to create histogram of doublet scores'
        },
        'layer_key':{
            'type':str,
            'description':'If provided, use this layer instead of .X',
            'default':None
        }

    }
    
    def __init__(self, name, params):
        super().__init__(name, params)
        self.logger = logging.getLogger(f"Module.{name}")
        
        # Define required inputs and outputs
        self.required_inputs = ["data"]
        self.outputs = ["data"]  # We'll update the same AnnData object
        
    def run(self, data_context):
        """
        Detect potential doublets using scanpy's scrublet wrapper.
        
        Args:
            data_context: The shared data context
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            # Get the data
            adata = data_context.get("data")
            
            self.logger.info("Running doublet detection")
            
            # Set active layer if specified:
            layer_key = self.params.get('layer_key',None)
            if layer_key:
                adata = set_active_layer(adata, layer_key)

            # Run scrublet via scanpy
            sc.pp.scrublet(
                adata,
                expected_doublet_rate=self.params.get('expected_doublet_rate', 0.05),
                min_counts=self.params.get('min_counts', 3),
                min_cells=self.params.get('min_cells', 3),
                n_neighbors=self.params.get('n_neighbors', 30),
                sim_doublet_ratio=self.params.get('sim_doublet_ratio', 2.0)
            )
            
            # Log results
            n_doublets = adata.obs['predicted_doublet'].sum()
            doublet_rate = n_doublets / len(adata.obs)
            self.logger.info(f"Detected {n_doublets} doublets ({doublet_rate:.2%} of cells)")
            
            # Create plot if requested
            if self.params.get('create_plots', True):
                fig = sc.pl.scrublet_score_distribution(adata, return_fig=True)
                img_path = self.save_figure(self.name, fig)
                
                data_context.add_figure(
                    module_name=self.name,
                    title="Doublet Score Distribution",
                    description="Distribution of doublet scores with threshold shown as a vertical line.",
                    image_path=img_path,
                    caption=f"Detected {n_doublets} doublets ({doublet_rate:.2%} of cells)"
                )
            
            # Update the data context
            data_context.set("data", adata)
            
            self.logger.info("Doublet detection completed successfully")
            return True
            
        except Exception as e:
            self.logger.error(f"Error in doublet detection: {e}", exc_info=True)
            return False