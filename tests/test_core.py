# tests/test_core.py

import os
import sys
import unittest
import tempfile
import shutil
import yaml
import scanpy as sc
import numpy as np
import pandas as pd

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sc_pipeline.core.module import AnalysisModule
from sc_pipeline.core.data_context import DataContext
from sc_pipeline.core.config import ConfigParser
from sc_pipeline.core.executor import PipelineExecutor

class TestModule(AnalysisModule):
    """A simple test module for testing the core functionality."""
    
    def __init__(self, name, params):
        super().__init__(name, params)
        self.required_inputs = params.get('required_inputs', [])
        self.outputs = params.get('outputs', [])
        
    def run(self, data_context):
        # Simply copy an input to an output or create a new value
        for output in self.outputs:
            if output in self.required_inputs:
                # Copy the input to the output
                data_context.set(output, data_context.get(output))
            else:
                # Create a new value
                data_context.set(output, self.params.get('value', 'test_value'))
        return True

class CoreTest(unittest.TestCase):
    """Tests for the core pipeline components."""
    
    def setUp(self):
        # Create a temporary directory for test outputs
        self.test_dir = tempfile.mkdtemp()
        
    def tearDown(self):
        # Clean up the temporary directory
        shutil.rmtree(self.test_dir)
        
    def test_data_context(self):
        """Test the DataContext functionality."""
        # Create a data context
        context = DataContext(self.test_dir)
        
        # Test setting and getting values
        context.set('test_key', 'test_value')
        self.assertEqual(context.get('test_key'), 'test_value')
        
        # Test checkpoint functionality
        context.save_checkpoint('test_checkpoint')
        
        # Modify the data
        context.set('test_key', 'modified_value')
        self.assertEqual(context.get('test_key'), 'modified_value')
        
        # Load the checkpoint
        context.load_checkpoint('test_checkpoint')
        self.assertEqual(context.get('test_key'), 'test_value')
        
    def test_config_parser(self):
        """Test the ConfigParser functionality."""
        # Create a test config
        config = {
            'pipeline': {
                'name': 'test_pipeline',
                'output_dir': self.test_dir
            },
            'modules': [
                {
                    'name': 'test_module',
                    'type': 'TestModule',
                    'params': {
                        'value': 'test_value'
                    }
                }
            ]
        }
        
        # Write the config to a file
        config_path = os.path.join(self.test_dir, 'test_config.yaml')
        with open(config_path, 'w') as f:
            yaml.dump(config, f)
            
        # Parse the config
        parser = ConfigParser()
        parsed_config = parser.parse(config_path)
        
        # Verify the parsed config
        self.assertEqual(parsed_config['pipeline']['name'], 'test_pipeline')
        self.assertEqual(len(parsed_config['modules']), 1)
        self.assertEqual(parsed_config['modules'][0]['name'], 'test_module')
        
    def test_pipeline_executor(self):
        """Test the PipelineExecutor functionality."""
        # Create a test config
        config = {
            'pipeline': {
                'name': 'test_pipeline',
                'output_dir': self.test_dir
            },
            'modules': [
                {
                    'name': 'module1',
                    'type': 'TestModule',
                    'params': {
                        'outputs': ['value1'],
                        'value': 'test_value1'
                    }
                },
                {
                    'name': 'module2',
                    'type': 'TestModule',
                    'params': {
                        'required_inputs': ['value1'],
                        'outputs': ['value2'],
                        'value': 'test_value2'
                    }
                }
            ]
        }
        
        # Write the config to a file
        config_path = os.path.join(self.test_dir, 'test_config.yaml')
        with open(config_path, 'w') as f:
            yaml.dump(config, f)
            
        # Create a pipeline executor
        executor = PipelineExecutor(config_path)
        
        # Register the test module
        executor.register_module_type('TestModule', TestModule)
        
        # Run the pipeline
        success = executor.run()
        
        # Verify the pipeline ran successfully
        self.assertTrue(success)
        
        # Check that the outputs were created
        self.assertEqual(executor.data_context.get('value1'), 'test_value1')
        self.assertEqual(executor.data_context.get('value2'), 'test_value2')
        
    def test_minimal_pipeline(self):
        """Test a minimal pipeline with real data."""
        # Skip this test if test data is not available
        test_data_path = os.path.join(os.path.dirname(__file__), 'test_data', '10k_PBMC_3p_nextgem_Chromium_Controller_filtered_feature_bc_matrix.h5')
        if not os.path.exists(test_data_path):
            self.skipTest("Test data not available")
        
        # Create a test config that uses our DataLoading and QCMetrics modules
        config = {
            'pipeline': {
                'name': 'minimal_test_pipeline',
                'output_dir': self.test_dir
            },
            'modules': [
                {
                    'name': 'data_loading',
                    'type': 'DataLoading',
                    'params': {
                        'file_path': test_data_path
                    }
                },
                {
                    'name': 'qc_metrics',
                    'type': 'QCMetrics',
                    'params': {
                        'mito_pattern': '^MT-'
                    }
                }
            ]
        }
        
        # Write the config to a file
        config_path = os.path.join(self.test_dir, 'minimal_test_config.yaml')
        with open(config_path, 'w') as f:
            yaml.dump(config, f)
            
        # Create a pipeline executor
        executor = PipelineExecutor(config_path)
        
        # Run the pipeline
        success = executor.run()
        
        # Verify the pipeline ran successfully
        self.assertTrue(success)
        
        # Check that the outputs were created
        adata = executor.data_context.get('loaded_data')
        self.assertIsNotNone(adata)
        self.assertIn('n_genes_by_counts', adata.obs.columns)
        self.assertIn('pct_counts_mt', adata.obs.columns)

if __name__ == '__main__':
    unittest.main()