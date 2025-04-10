# tests/test_ambient_removal.py

import os
import sys
import unittest
import scanpy as sc
import numpy as np
import tempfile
import shutil
import logging

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sc_pipeline.core.data_context import DataContext
from sc_pipeline.modules.ambientrnaremoval import AmbientRNARemoval
from sc_pipeline.modules.dataloading import DataLoading

class AmbientRNARemovalTest(unittest.TestCase):
    """Tests for the AmbientRNARemoval module."""
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        
        # Create a test data context
        self.data_context = DataContext(self.test_dir)

        # Test logging
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        
        # Define paths to test data
        self.filtered_h5_path = os.path.join(
            os.path.dirname(__file__), 
            'test_data', 
            '10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.h5'
        )
        
        self.raw_h5_path = os.path.join(
            os.path.dirname(__file__), 
            'test_data', 
            '10k_PBMC_3p_nextgem_Chromium_Controller_raw_feature_bc_matrix.h5'
        )
        
        # Skip if test files don't exist
        if not os.path.exists(self.filtered_h5_path) or not os.path.exists(self.raw_h5_path):
            self.skipTest("Test data files not found")
        
        # First load the data using the DataLoading module
        data_loader = DataLoading('data_loading', {'file_path': self.filtered_h5_path})
        success = data_loader.run(self.data_context)
        if not success:
            self.fail("Failed to load test data")
        
    def tearDown(self):
        # Clean up the temporary directory
        shutil.rmtree(self.test_dir)
        
    def test_ambient_removal(self):
        """Test the AmbientRNARemoval module with real data."""
        
        # Create the module with our test files
        module = AmbientRNARemoval('ambient_removal', {
            'raw_counts_path': self.raw_h5_path,
            'filtered_counts_path': self.filtered_h5_path,
            'ndims': 20,  # Use fewer dimensions for speed
            'resolution': 0.5,
            'replace_main_matrix': True
        })
        
        # Run the module
        success = module.run(self.data_context)
        
        # Check that it ran successfully
        self.assertTrue(success, "Ambient RNA removal module failed to run")
        
        # Get the processed data
        adata = self.data_context.get('data')
        
        # Check that layers were created
        self.assertIn('original', adata.layers, "Original layer not created")
        self.assertIn('soupx_corrected', adata.layers, "SoupX corrected layer not created")
        
        # Check that the main matrix was replaced
        self.assertTrue(np.array_equal(adata.X.toarray(), adata.layers['soupx_corrected'].toarray()), 
                        "Main matrix was not replaced with corrected counts")
        
        # Check that the contamination estimate was added to uns
        self.assertIn('ambient_removal', adata.uns, "Ambient removal metadata not added")
        self.assertIn('estimated_contamination', adata.uns['ambient_removal'], 
                     "Contamination estimate not found")
        
        # Check that the contamination value is reasonable (typically between 0.01 and 0.3)
        contamination = adata.uns['ambient_removal']['estimated_contamination']
        self.assertGreaterEqual(contamination, 0.01, "Contamination estimate abnormally low")
        self.assertLessEqual(contamination, 0.4, "Contamination estimate abnormally high")
        
        # Check that the counts were actually modified (i.e., original and corrected are different)
        # This just checks that at least some values are different
        diff = np.sum(np.abs(adata.layers['original'].toarray() - adata.layers['soupx_corrected'].toarray()))
        self.assertGreater(diff, 0, "No difference between original and corrected counts")
        
        # Check that the corrected matrix has fewer counts overall (ambient RNA removal should reduce counts)
        original_sum = adata.layers['original'].sum()
        corrected_sum = adata.layers['soupx_corrected'].sum()
        self.assertLess(corrected_sum, original_sum, 
                       "Corrected matrix doesn't have fewer counts than original")

if __name__ == '__main__':
    unittest.main()