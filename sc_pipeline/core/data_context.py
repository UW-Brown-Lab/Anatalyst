# sc_pipeline/core/data_context.py

import os
import pickle
import logging

class DataContext:
    """
    Stores and manages data shared between pipeline modules.
    Provides checkpoint capabilities
    """

    def __init__(self, checkpoint_dir=None):
        self._data = {}
        self.checkpoint_dir = checkpoint_dir

        if checkpoint_dir and not os.path.exists(checkpoint_dir):
            os.makedirs(checkpoint_dir)

        self.logger = logging.getLogger("DataContext")

    def __contains__(self, key):
        return key in self._data
    
    def get(self, key, default=None):
        """Get data by key."""
        return self._data.get(key, default)
    
    def set(self, key, value):
        """Store data by key."""
        self._data[key] = value

    def keys(self):
        """Return all available data keys."""
        return self._data.keys()
    
    def save_checkpoint(self, checkpoint_name):
        """Save the current state to a checkpoint file."""
        if not self.checkpoint_dir:
            self.logger.warning("No checkpoint directory specified, skipping checkpoint.")
            return False
        
        checkpoint_path = os.path.join(self.checkpoint_dir, f"{checkpoint_name}.pkl")

        try:
            with open(checkpoint_path, 'wb') as f:
                pickle.dump(self._data, f)
            self.logger.info(f"Saved checkpoint: {checkpoint_path}")
            return True
        except Exception as e:
            self.logger.error(f"Failed to save checkpoint: {e}")
            return False
        
    def load_checkpoint(self, checkpoint_name):
        """Load data from a checkpoint file."""
        if not self.checkpoint_dir:
            self.logger.warning("No checkpoint directory specified.")
            return False
        
        checkpoint_path = os.path.join(self.checkpoint_dir, f"{checkpoint_name}.pkl")

        if not os.path.exists(checkpoint_path):
            self.logger.error(f"Checkpoint file not found: {checkpoint_path}")
            return False

        try:
            with open(checkpoint_path, 'rb') as f:
                self._data = pickle.load(f)
            self.logger.info(f"Loaded checkpoint: {checkpoint_path}")
            return True
        except Exception as e:
            self.logger.error(f"Failed to load checkpoint: {e}")
            return False