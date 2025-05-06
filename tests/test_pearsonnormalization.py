# tests/test_pearson_normalization.py

import os
import sys
import unittest
import scanpy as sc
import numpy as np
import scipy.sparse
import tempfile
import shutil
import logging
import random

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sc_pipeline.core.data_context import DataContext
from sc_pipeline.modules.pearsonnormalization import PearsonResidualsNormalization
from sc_pipeline.modules.dataloading import DataLoading

class PearsonResidualsNormalizationTest(unittest.TestCase):
    """Tests for the PearsonResidualsNormalization module."""
    
    def setUp(self):
        # Set random seed for reproducibility
        random.seed(42)
        np.random.seed(42)
        
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        
        # Create a test data context
        self.data_context = DataContext(self.test_dir)
        
        # Set up logging
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        
        # Define path to test data
        self.filtered_h5_path = os.path.join(
            os.path.dirname(__file__), 
            'test_data', 
            '10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.h5'
        )
        
        # Skip if test file doesn't exist
        if not os.path.exists(self.filtered_h5_path):
            self.skipTest("Test data file not found")
        
        # First load the data using the DataLoading module
        data_loader = DataLoading('data_loading', {'file_path': self.filtered_h5_path})
        success = data_loader.run(self.data_context)
        if not success:
            self.fail("Failed to load test data")
        
    def tearDown(self):
        # Clean up the temporary directory
        shutil.rmtree(self.test_dir)
        
    def test_pearson_normalization(self):
        """Test the PearsonResidualsNormalization module with default parameters."""
        
        # Create the module with default parameters
        module = PearsonResidualsNormalization('pearson_normalization', {
            'theta': 100.0,
            'use_highly_variable': True,
            'create_plots': False  # Disable plotting for testing
        })
        
        # Run the module
        success = module.run(self.data_context)
        
        # Check that it ran successfully
        self.assertTrue(success, "Pearson normalization module failed to run")
        
        # Get the processed data
        adata = self.data_context.get('data')
        
        # Check that pearson residual layer was created
        self.assertIn('pearson_residuals', adata.layers, "Pearson residuals layer not created")
        
        # Check that the pearson residuals have expected properties
        pearson_residuals = adata.layers['pearson_residuals']
        
        # Verify that pearson residuals contain finite values
        if scipy.sparse.issparse(pearson_residuals):
            self.assertFalse(np.any(~np.isfinite(pearson_residuals.data)), "Pearson residuals contain non-finite values")
        else:
            self.assertFalse(np.any(~np.isfinite(pearson_residuals)), "Pearson residuals contain non-finite values")
        
        # Check if the highly variable genes were computed
        if adata.var.get('highly_variable') is not None:
            self.assertTrue(adata.var['highly_variable'].any(), "No highly variable genes detected")
            
        # Verify that the normalization changes the data
        if scipy.sparse.issparse(adata.X) and scipy.sparse.issparse(pearson_residuals):
            self.assertFalse(np.allclose(adata.X.data, pearson_residuals.data), 
                            "Pearson residuals are identical to original data")
        
    def test_pearson_normalization_with_layer(self):
        """Test the PearsonResidualsNormalization module with a specific layer."""
        
        # Get the initial data
        adata = self.data_context.get('data')
        
        # Create a layer for testing
        adata.layers['test_layer'] = adata.X.copy()
        
        # Create the module with layer parameter
        module = PearsonResidualsNormalization('pearson_normalization_layer', {
            'theta': 50,  # Different theta for testing
            'use_highly_variable': True,
            'layer_key': 'test_layer',
            'create_plots': False
        })
        
        # Run the module
        success = module.run(self.data_context)
        
        # Check that it ran successfully
        self.assertTrue(success, "Pearson normalization module with layer failed to run")
        
        # Get the processed data
        adata = self.data_context.get('data')
        
        # Check that pearson residual layer was created
        self.assertIn('pearson_residuals', adata.layers, "Pearson residuals layer not created")
        
        # Check that the normalization was applied to the correct layer
        if 'test_layer' in adata.layers and 'pearson_residuals' in adata.layers:
            # Verify that the normalization was actually applied
            # (Exact parameters depend on how Scanpy's normalize_pearson_residuals works,
            # so we just check that the values are different)
            # Compare a sample of the data rather than converting whole matrices
            # This is more memory-efficient for large datasets
            if scipy.sparse.issparse(adata.layers['test_layer']) and scipy.sparse.issparse(adata.layers['pearson_residuals']):
                # Sample some non-zero elements
                test_indices = adata.layers['test_layer'].nonzero()
                sample_size = min(100, len(test_indices[0]))
                if sample_size > 0:
                    sample_indices = (test_indices[0][:sample_size], test_indices[1][:sample_size])
                    test_values = adata.layers['test_layer'][sample_indices]
                    pearson_values = adata.layers['pearson_residuals'][sample_indices]
                    self.assertFalse(np.allclose(test_values, pearson_values), 
                                    "Sampled Pearson residuals are identical to original data")
            else:
                # For dense matrices, we can compare a random subset
                n_cells, n_genes = adata.shape
                if n_cells > 0 and n_genes > 0:
                    sample_rows = np.random.choice(n_cells, min(10, n_cells), replace=False)
                    sample_cols = np.random.choice(n_genes, min(10, n_genes), replace=False)
                    test_sample = adata.layers['test_layer'][sample_rows][:, sample_cols]
                    pearson_sample = adata.layers['pearson_residuals'][sample_rows][:, sample_cols]
                    self.assertFalse(np.array_equal(test_sample, pearson_sample),
                                    "Sampled Pearson residuals are identical to original data")
        else:
            self.fail("Required layers not found in the AnnData object")

if __name__ == '__main__':
    unittest.main()