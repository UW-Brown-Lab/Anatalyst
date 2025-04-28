# tests/test_filtering.py

import os
import sys
import unittest
import tempfile
import shutil
import logging
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sc_pipeline.core.data_context import DataContext
from sc_pipeline.modules.filtering import Filtering, FilterManager
from sc_pipeline.modules.dataloading import DataLoading
from sc_pipeline.modules.qcmetrics import QCMetrics

class FilteringTest(unittest.TestCase):
    """Tests for the Filtering module."""
    
    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        
        # Create a test data context
        self.data_context = DataContext(self.test_dir)

        # Test logging
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        
        # Define path to test data
        self.h5_path = os.path.join(
            os.path.dirname(__file__), 
            'test_data', 
            '10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.h5'
        )
        
        # Skip if test file doesn't exist
        if not os.path.exists(self.h5_path):
            self.skipTest("Test data file not found")
        
        # First load the data using the DataLoading module
        data_loader = DataLoading('data_loading', {'file_path': self.h5_path})
        success = data_loader.run(self.data_context)
        if not success:
            self.fail("Failed to load test data")
            
        # Then calculate QC metrics
        qc_module = QCMetrics('qc_metrics', {
            'mito_pattern': '^MT-',
            'ribo_pattern': '^RP[SL]',
            'create_plots': False
        })
        success = qc_module.run(self.data_context)
        if not success:
            self.fail("Failed to calculate QC metrics")
        
    def tearDown(self):
        # Clean up the temporary directory
        shutil.rmtree(self.test_dir)
        
    def test_filter_manager(self):
        """Test the FilterManager class."""
        # Get the data
        adata = self.data_context.get('data')
        
        # Save original cell count
        original_cells = adata.n_obs
        
        # Create a filter manager
        filter_manager = FilterManager()
        
        # Define test filters
        test_filters = {
            'n_genes_by_counts': {
                'type': 'numeric',
                'min': 200,  # Set a minimum threshold that should filter some cells
                'max': None,
                'outlier_detection': False
            },
            'pct_counts_mt': {
                'type': 'numeric',
                'min': None,
                'max': 10,  # Set a maximum threshold that should filter some cells
                'outlier_detection': False
            }
        }
        
        # Apply filters
        filtered_adata = filter_manager.apply_filters(adata, test_filters, inplace=False)
        
        # Verify cells were filtered
        self.assertLess(filtered_adata.n_obs, original_cells, 
                       "No cells were filtered when some should have been")
        
        # Verify the right cells were filtered based on thresholds
        self.assertTrue(all(filtered_adata.obs['n_genes_by_counts'] >= 200),
                       "Some cells with n_genes_by_counts < 200 were not filtered")
        
        self.assertTrue(all(filtered_adata.obs['pct_counts_mt'] <= 10),
                       "Some cells with pct_counts_mt > 10 were not filtered")
        
    def test_outlier_detection(self):
        """Test the outlier detection functionality."""
        # Get the data
        adata = self.data_context.get('data')
        
        # Create a filter manager
        filter_manager = FilterManager()
        
        # Define test filters with outlier detection
        test_filters = {
            'n_genes_by_counts': {
                'type': 'numeric',
                'outlier_detection': True,
                'remove_outliers': True,
                'n_mads': 3.0
            }
        }
        
        # Apply filters
        filtered_adata = filter_manager.apply_filters(adata, test_filters, inplace=False)
        
        # Calculate expected outliers manually
        values = adata.obs['n_genes_by_counts']
        median = np.median(values)
        mad = median_abs_deviation(values, nan_policy='omit')
        expected_outliers = (values < median - 3.0 * mad) | (values > median + 3.0 * mad)
        expected_outlier_count = np.sum(expected_outliers)
        
        # Calculate actual filtered cell count
        actual_filtered = adata.n_obs - filtered_adata.n_obs
        
        # Verify the outlier detection worked correctly
        self.assertEqual(actual_filtered, expected_outlier_count,
                        f"Expected {expected_outlier_count} outliers, but got {actual_filtered}")
        
    def test_boolean_filter(self):
        """Test the boolean filter functionality."""
        # Get the data
        adata = self.data_context.get('data')
        
        # Create a boolean column for testing
        adata.obs['test_bool'] = adata.obs['pct_counts_mt'] > 5
        
        # Create a filter manager
        filter_manager = FilterManager()
        
        # Define test filters
        test_filters = {
            'test_bool': {
                'type': 'boolean',
                'keep': False  # Keep cells where test_bool is False
            }
        }
        
        # Apply filters
        filtered_adata = filter_manager.apply_filters(adata, test_filters, inplace=False)
        
        # Count cells that should be kept
        expected_kept = np.sum(~adata.obs['test_bool'])
        
        # Verify the filter worked correctly
        self.assertEqual(filtered_adata.n_obs, expected_kept,
                        f"Expected {expected_kept} cells, but got {filtered_adata.n_obs}")
        
    def test_categorical_filter(self):
        """Test the categorical filter functionality."""
        # Get the data
        adata = self.data_context.get('data')
        
        # Create a categorical column for testing
        adata.obs['test_cat'] = pd.Categorical(
            np.where(adata.obs['pct_counts_mt'] > 7.5, 'high_mt', 
                    np.where(adata.obs['pct_counts_mt'] > 2.5, 'medium_mt', 'low_mt'))
        )
        
        # Create a filter manager
        filter_manager = FilterManager()
        
        # Define test filters
        test_filters = {
            'test_cat': {
                'type': 'categorical',
                'keep': ['low_mt', 'medium_mt']  # Filter out 'high_mt' cells
            }
        }
        
        # Apply filters
        filtered_adata = filter_manager.apply_filters(adata, test_filters, inplace=False)
        
        # Count cells that should be kept
        expected_kept = np.sum(adata.obs['test_cat'].isin(['low_mt', 'medium_mt']))
        
        # Verify the filter worked correctly
        self.assertEqual(filtered_adata.n_obs, expected_kept,
                        f"Expected {expected_kept} cells, but got {filtered_adata.n_obs}")
        
        # Verify no high_mt cells remain
        self.assertEqual(np.sum(filtered_adata.obs['test_cat'] == 'high_mt'), 0,
                        "Some 'high_mt' cells were not filtered out")
        
    def test_filtering_module(self):
        """Test the complete Filtering module with default parameters."""
        # Create the module
        module = Filtering('filtering', {
            'filters': {
                'n_genes_by_counts': {
                    'type': 'numeric',
                    'min': 200,
                    'max': None,
                    'outlier_detection': True,
                    'remove_outliers': True
                },
                'pct_counts_mt': {
                    'type': 'numeric',
                    'min': None,
                    'max': 10,
                    'outlier_detection': True,
                    'remove_outliers': True,
                    'two_sided': False  # Only detect high outliers for MT%
                }
            },
            'create_plots': False
        })
        
        # Get original count
        original_cells = self.data_context.get('data').n_obs
        
        # Run the module
        success = module.run(self.data_context)
        
        # Check that it ran successfully
        self.assertTrue(success, "Filtering module failed to run")
        
        # Get the filtered data
        adata = self.data_context.get('data')
        
        # Verify cells were filtered
        self.assertLess(adata.n_obs, original_cells, 
                       "No cells were filtered when some should have been")
        
        # Verify stats were added to the AnnData object
        self.assertIn('filtering', adata.uns, 
                     "Filtering stats not added to AnnData.uns")
        
        self.assertIn('original_cells', adata.uns['filtering'], 
                     "original_cells not in filtering stats")
        
        self.assertEqual(adata.uns['filtering']['original_cells'], original_cells,
                        "original_cells count is incorrect")
        
        # Verify the right cells were filtered based on thresholds
        self.assertTrue(all(adata.obs['n_genes_by_counts'] >= 200),
                       "Some cells with n_genes_by_counts < 200 were not filtered")
        
        self.assertTrue(all(adata.obs['pct_counts_mt'] <= 10),
                       "Some cells with pct_counts_mt > 10 were not filtered")
        
if __name__ == '__main__':
    unittest.main()