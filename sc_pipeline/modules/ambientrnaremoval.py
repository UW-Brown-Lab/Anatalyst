# sc_pipeline/modules/ambient_removal.py

import os
import logging
import tempfile
import scanpy as sc
import scipy.sparse as sp
from sc_pipeline.core.module import AnalysisModule
from sc_pipeline.utils.r_bridge import RBridge

class AmbientRNARemoval(AnalysisModule):
    """Module for removing ambient RNA contamination using SoupX."""
    
    def __init__(self, name, params):
        super().__init__(name, params)
        self.logger = logging.getLogger(f"Module.{name}")
        
        # Define required inputs and outputs
        self.required_inputs = ["loaded_data"]
        self.outputs = ["loaded_data"]  # We now modify the same object instead of creating a new one
        
    def run(self, data_context):
        """
        Remove ambient RNA contamination using SoupX.
        
        Args:
            data_context: The shared data context
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            # Get the data
            adata = data_context.get("loaded_data")
            
            # Get the raw and filtered counts paths
            raw_counts_path = self.params.get('raw_counts_path')
            filtered_counts_path = self.params.get('filtered_counts_path')
            
            if not raw_counts_path:
                self.logger.error("No raw counts path specified")
                return False
                
            if not filtered_counts_path:
                self.logger.error("No filtered counts path specified")
                return False
                
            if not os.path.exists(raw_counts_path):
                self.logger.error(f"Raw counts file not found: {raw_counts_path}")
                return False
                
            if not os.path.exists(filtered_counts_path):
                self.logger.error(f"Filtered counts file not found: {filtered_counts_path}")
                return False
                
            # Create R bridge
            r_bridge = RBridge()
            
            # Prepare arguments for the R script
            args = {
                'raw_counts_path': raw_counts_path,
                'filtered_counts_path': filtered_counts_path,
                'modality': self.params.get('modality', 'Gene Expression'),
                'ndims': self.params.get('ndims', 30),
                'resolution': self.params.get('resolution', 0.8)
            }
            
            # Prepare any data to pass to R
            data = {}
            
            # Add optional cluster information if available
            if 'clusters' in adata.obs:
                self.logger.info("Using existing cluster assignments")
                data['clusters'] = adata.obs['clusters'].to_dict()
            
            # Run the R script
            self.logger.info("Running SoupX with Seurat clustering...")
            results = r_bridge.run_r_script('run_soupx.R', args, data)
            
            # Check if execution was successful
            if not results.get('success', False):
                self.logger.error(f"SoupX failed: {results.get('error', 'Unknown error')}")
                return False
            
            # Log messages from R
            for msg in results.get('messages', []):
                self.logger.info(f"SoupX: {msg}")
            
            # Get the corrected counts
            corrected_counts = results.get('corrected_counts')
            if corrected_counts is None:
                self.logger.error("No corrected counts returned from SoupX")
                return False
            
            # Ensure the corrected counts are sparse
            if not sp.issparse(corrected_counts):
                corrected_counts = sp.csr_matrix(corrected_counts)
            
            # Save the original counts as a layer if not already done
            if 'original' not in adata.layers:
                self.logger.info("Storing original counts in 'original' layer")
                adata.layers['original'] = adata.X.copy()
            
            # Add the SoupX-corrected counts as a new layer
            self.logger.info("Adding SoupX-corrected counts as 'soupx_corrected' layer")
            adata.layers['soupx_corrected'] = corrected_counts
            
            # Replace the main matrix with the corrected counts if specified
            if self.params.get('replace_main_matrix', True):
                self.logger.info("Replacing main matrix with SoupX-corrected counts")
                adata.X = adata.layers['soupx_corrected']
            
            # Add SoupX as a processing step in the unstructured data
            adata.uns['ambient_rna_removed'] = True
            
            # Add some metadata about the process
            contamination = results.get('contamination_fraction', 0)
            self.logger.info(f"Estimated contamination fraction: {contamination}")
            
            adata.uns['ambient_removal'] = {
                'method': 'SoupX',
                'modality': self.params.get('modality', 'Gene Expression'),
                'estimated_contamination': contamination,
                'clustering': 'seurat' if 'clusters' not in data else 'provided'
            }
            
            # No need to update data_context since we modified the existing object
            
            self.logger.info("Ambient RNA removal completed successfully")
            return True
                
        except Exception as e:
            self.logger.error(f"Error in ambient RNA removal: {e}", exc_info=True)
            return False