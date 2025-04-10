# sc_pipeline/core/config.py

import yaml
import logging
from pathlib import Path

class ConfigParser:
    """Parses and validates pipeline configuration files."""

    def __init__(self):
        self.logger = logging.getLogger("ConfigParser")

    def parse(self, config_file):
        """
        Parse a YAML configuration file.

        Args:
            config_file: Path to the YAML configuration file

        Returns:
            dict: Parsed Configuration
        """
        config_path = Path(config_file)

        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_file}")
        
        try:
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)

                # Validate the configuration
                self._validate_config(config)
                return config
            
        except yaml.YAMLError as e:
            self.logger.error(f"Error parsing YAML file: {e}")
            raise

    def _validate_config(self, config):
        """Validate the configuration structure."""
        # Check that the required top-level sections exist
        required_sections = ['pipeline','modules']
        for section in required_sections:
            if section not in config:
                raise ValueError(f"Missing required section in config: {section}")
            
        # Check pipeline section:
        if 'name' not in config['pipeline']:
            raise ValueError("Pipeline must have a name.")
        
        if 'output_dir' not in config['pipeline']:
            raise ValueError("Pipeline must specify an output directory")
        

        # Check modules section:
        if not isinstance(config['modules'], list):
            raise ValueError("Modules must be a list")
        
        for idx, module in enumerate(config['modules']):
            if 'name' not in module:
                raise ValueError(f"Module at index {idx} missing name")
            if 'type' not in module:
                raise ValueError(f"Module {module.get('name', f'at index {idx}')} missing type")