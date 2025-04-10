# sc_pipeline/modules/ambient_removal.py

import os
import logging
import tempfile
from bleach import clean
import scanpy as sc
import scipy.io
import scipy.sparse as sp
import json
from sc_pipeline.core.module import AnalysisModule
from sc_pipeline.utils.r_bridge import RBridge

class AmbientRNARemoval(AnalysisModule):
    """Module for removing ambient RNA contamination using SoupX."""
    
    def __init__(self, name, params):
        super().__init__(name, params)
        self.logger = logging.getLogger(f"Module.{name}")
        
        # Define required inputs and outputs
        self.required_inputs = ["data"]
        self.outputs = ["data"]  # We now modify the same object instead of creating a new one
        
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
            adata = data_context.get("data")
            
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
            r_memory_limit_gb = data_context.get("R_MEMORY_LIMIT_GB", 8)
            self.logger.info(f"Using R memory limit of {r_memory_limit_gb} GB")
            r_bridge = RBridge(cleanup=True, memory_limit_gb=r_memory_limit_gb)
            
            # Prepare arguments for the R script
            args = {
                'raw_counts_path': raw_counts_path,
                'filtered_counts_path': filtered_counts_path,
                'output_matrix':'corrected_counts.mtx',
                'output_metadata': 'soupx_metadata.json',
                'modality': self.params.get('modality', 'Gene Expression'),
                'ndims': self.params.get('ndims', 30),
                'resolution': self.params.get('resolution', 0.8)
            }
            
            # Run the R script
            self.logger.info("Running SoupX with Seurat clustering...")
            success, stdout, stderr = r_bridge.run_r_script('run_soupx.R', args)
            
            # Check if execution was successful
            if not success:
                self.logger.error(f"SoupX failed: {stderr}")
                return False
            
            # Log output from R
            for line in stdout.splitlines():
                self.logger.info(f"SoupX: {line}")
            
            # Get the output filepaths
            corrected_mtx_path = r_bridge.get_workspace_path('corrected_counts.mtx')
            metadata_path = r_bridge.get_workspace_path('soupx_metadata.json')
            features_path = r_bridge.get_workspace_path('features.txt')
            barcodes_path = r_bridge.get_workspace_path('barcodes.txt')

            if not os.path.exists(corrected_mtx_path):
                self.logger.error(f"Can not located SoupX-corrected matrix at path: {corrected_mtx_path}")
                return False
        
            if not os.path.exists(metadata_path):
                self.logger.error(f"Can not located SoupX metadata at path: {metadata_path}")
                return False
            
            if not os.path.exists(features_path):
                self.logger.error(f"Can not located SoupX features at path: {features_path}")
                return False
            
            if not os.path.exists(barcodes_path):
                self.logger.error(f"Can not located SoupX barcodes at path: {barcodes_path}")
                return False
            
            # Read the metadata
            with open(metadata_path, 'r') as f:
                soupx_metadata = json.load(f)
            
            # Read the corrected matrix
            corrected_counts = scipy.io.mmread(corrected_mtx_path)
            corrected_counts = sp.csr_matrix(corrected_counts)

            # Log the current shape of the matrix
            self.logger.info(f"Original matrix shape from SoupX: {corrected_counts.shape}")
            self.logger.info(f"Expected shape for AnnData: {adata.shape}")
            
            # Check if we need to transpose the matrix
            if corrected_counts.shape[0] != adata.shape[0]:
                self.logger.info("Matrix dimensions don't match AnnData - transposing")
                corrected_counts = corrected_counts.transpose()
                self.logger.info(f"New shape after transposition: {corrected_counts.shape}")

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
                adata.X = adata.layers['soupx_corrected'].copy()
            
            # Add SoupX as a processing step in the unstructured data
            adata.uns['ambient_rna_removed'] = True
            
            # Add some metadata about the process
            contamination = soupx_metadata['contamination_fraction']
            self.logger.info(f"Estimated contamination fraction: {contamination}")
            
            adata.uns['ambient_removal'] = {
                'method': 'SoupX',
                'modality': self.params.get('modality', 'Gene Expression'),
                'estimated_contamination': contamination,
                'clustering': f"Seurat at Res: {self.params.get('resolution', 0.8)} and ndims: {self.params.get('ndims', 30)}"
            }
            
            # No need to update data_context since we modified the existing object
            
            self.logger.info("Ambient RNA removal completed successfully")
            return True
                
        except Exception as e:
            self.logger.error(f"Error in ambient RNA removal: {e}", exc_info=True)
            return False