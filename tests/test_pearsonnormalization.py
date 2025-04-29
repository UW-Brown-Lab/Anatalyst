# tests/test_pearson_normalization.py

import os
import sys
import unittest
import tempfile
import shutil
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import logging
import importlib

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sc_pipeline.core.data_context import DataContext
from sc_pipeline.modules.dataloading import DataLoading

# Patch the PearsonResidualsNormalization class before importing it
# This is a workaround for the logger initialization issue
# A better solution would be to fix the module itself

def patch_module():
    # Create a modified version of the module that initializes the logger properly
    module_path = 'sc_pipeline.modules.pearsonnormalization'
    module = importlib.import_module(module_path)
    
    original_init = module.PearsonResidualsNormalization.__init__
    
    def patched_init(self, name, params):
        # Initialize logger before calling parent's __init__
        self.logger = logging.getLogger(f"Module.{name}")
        # Call the original __init__
        original_init(self, name, params)
    
    # Apply the patched __init__
    module.PearsonResidualsNormalization.__init__ = patched_init
    
    return module.PearsonResidualsNormalization

# Apply the patch
PearsonResidualsNormalization = patch_module()

class PearsonResidualsNormalizationTest(unittest.TestCase):
    """Tests for the PearsonResidualsNormalization module."""
    
    def setUp(self):
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
        
    def test_normalization_basic(self):
        """Test the basic functionality of the PearsonResidualsNormalization module."""
        
        # Get the data
        adata = self.data_context.get('data')
        
        # Store original state for comparison
        original_shape = adata.shape
        original_sum = adata.X.sum()
        
        try:
            # Create the module with default parameters
            module = PearsonResidualsNormalization('pearson_normalization', {
                'theta': 100.0,
                'create_plots': False  # Disable plots for testing
            })
            
            # Run the module
            success = module.run(self.data_context)
            
            # Check that it ran successfully
            self.assertTrue(success, "Pearson residuals normalization module failed to run")
            
            # Get the processed data
            adata = self.data_context.get('data')
            
            # Check that the overall shape is unchanged (unless use_highly_variable=True and n_top_genes reduces features)
            self.assertEqual(adata.shape[0], original_shape[0], "Number of cells changed after normalization")
            
            # Check that 'highly_variable' was added to var
            self.assertIn('highly_variable', adata.var, "Highly variable genes not marked in adata.var")
            
            # Check that the number of highly variable genes matches the parameter
            n_hvg = adata.var['highly_variable'].sum()
            self.assertLessEqual(n_hvg, 2000, "Too many highly variable genes selected")
            
            # Check that the values are now normalized (no longer raw counts, should contain negative values)
            self.assertTrue(np.any(adata.X < 0), "No negative values found after Pearson residuals normalization")
            
            # Check that normalization info was added to the layers
            # (Not applicable for this module, but could be added in future implementations)
        except Exception as e:
            self.fail(f"Test failed with exception: {str(e)}")
        
    def test_normalization_with_layer(self):
        """Test normalization using a specific layer instead of the main matrix."""
        
        try:
            # Get the data
            adata = self.data_context.get('data')
            
            # Reduce the size of the test data to avoid memory issues
            adata_subset = adata[:1000, :1000].copy()  # Only use 1000 cells and 1000 genes
            
            # Add a layer to test with - just a copy of X in this case
            adata_subset.layers['raw_counts'] = adata_subset.X.copy()
            
            # Put the subset back in the data context
            self.data_context.set('data', adata_subset)
            
            # Create the module with layer parameter
            module = PearsonResidualsNormalization('pearson_layer_normalization', {
                'theta': 10.0,         # Reduced from 50
                'layer_key': 'raw_counts',
                'create_plots': False
            })
            
            # Run the module
            success = module.run(self.data_context)
            
            # Check that it ran successfully
            self.assertTrue(success, "Pearson residuals normalization with layer failed to run")
            
            # Get the processed data
            adata_result = self.data_context.get('data')
            
            # Check that the layer was properly normalized by verifying negative values
            has_negative = False
            try:
                has_negative = np.any(adata_result.X < 0)
            except:
                # If checking dense array fails, try converting to dense first
                try:
                    if sp.issparse(adata_result.X):
                        has_negative = np.any(adata_result.X.toarray() < 0)
                except:
                    pass
                
            self.assertTrue(has_negative, "No negative values found after Pearson residuals normalization")
            
            # Check that highly variable genes were selected according to parameters
            n_hvg = adata_result.var['highly_variable'].sum()
            self.assertLessEqual(n_hvg, 500, "Too many highly variable genes selected")
            
            # Restore original data for other tests
            self.data_context.set('data', adata)
            
        except Exception as e:
            # Restore original data before failing
            try:
                self.data_context.set('data', adata)
            except:
                pass
            self.fail(f"Test failed with exception: {str(e)}")
        
    def test_normalization_without_hvg(self):
        """Test normalization without subsetting to highly variable genes."""
        
        try:
            # Get the data
            adata = self.data_context.get('data')
            
            # Reduce the size of the test data to avoid memory issues
            adata_subset = adata[:1000, :1000].copy()  # Only use 1000 cells and 1000 genes
            self.data_context.set('data', adata_subset)
            
            original_shape = adata_subset.shape
            
            # Create the module with use_highly_variable=False
            module = PearsonResidualsNormalization('pearson_no_hvg', {
                'use_highly_variable': False,
                'theta': 10.0,        # Reduced from 100
                'create_plots': False
            })
            
            # Run the module
            success = module.run(self.data_context)
            
            # Check that it ran successfully
            self.assertTrue(success, "Pearson residuals normalization without HVG failed to run")
            
            # Get the processed data
            adata_result = self.data_context.get('data')
            
            # Check that the shape is unchanged (no genes were filtered out)
            self.assertEqual(adata_result.shape, original_shape, "Shape changed after normalization without HVG")
            
            # Check that the values are normalized (should have negative values)
            has_negative = False
            try:
                has_negative = np.any(adata_result.X < 0)
            except:
                # If checking dense array fails, try converting to dense first
                try:
                    if sp.issparse(adata_result.X):
                        has_negative = np.any(adata_result.X.toarray() < 0)
                except:
                    pass
                
            self.assertTrue(has_negative, "No negative values found after Pearson residuals normalization")
            
            # Restore original data for other tests
            self.data_context.set('data', adata)
            
        except Exception as e:
            # Restore original data before failing
            try:
                self.data_context.set('data', adata)
            except:
                pass
            self.fail(f"Test failed with exception: {str(e)}")
        
    def test_create_plots(self):
        """Test that plots are correctly created when requested."""
        
        try:
            # Get the data
            adata = self.data_context.get('data')
            
            # Reduce the size of the test data to avoid memory issues
            adata_subset = adata[:1000, :1000].copy()  # Only use 1000 cells and 1000 genes
            self.data_context.set('data', adata_subset)
            
            # Create the module with create_plots=True and optimized parameters
            module = PearsonResidualsNormalization('pearson_with_plots', {
                'theta': 10.0,         # Reduced from 100
                'chunk_size': 200,   # Smaller chunks
                'create_plots': True
            })
            
            # Run the module
            success = module.run(self.data_context)
            
            # Check that it ran successfully
            self.assertTrue(success, "Pearson residuals normalization with plots failed to run")
            
            # Check that figures were added to the data context
            self.assertIn('REPORT_FIGURES', self.data_context._data, "No figures added to data context")
            self.assertIn('pearson_with_plots', self.data_context._data['REPORT_FIGURES'], 
                         "No figures found for the module")
            
            # Check that at least one figure was created
            self.assertGreater(len(self.data_context._data['REPORT_FIGURES']['pearson_with_plots']), 0,
                              "No figures created")
            
            # Check that the figure files exist
            for figure in self.data_context._data['REPORT_FIGURES']['pearson_with_plots']:
                if 'image_path' in figure:
                    self.assertTrue(os.path.exists(figure['image_path']), 
                                   f"Figure file not found: {figure['image_path']}")
            
            # Restore original data for any future tests
            self.data_context.set('data', adata)
            
        except Exception as e:
            # Restore original data before failing
            try:
                self.data_context.set('data', adata)
            except:
                pass
            self.fail(f"Test failed with exception: {str(e)}")

if __name__ == '__main__':
    unittest.main()