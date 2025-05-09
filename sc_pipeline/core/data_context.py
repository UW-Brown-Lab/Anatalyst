# sc_pipeline/core/data_context.py

import os
import pickle
import logging

class DataContext:
    """
    Stores and manages data shared between pipeline modules.
    Provides checkpoint capabilities
    """

    def __init__(self, checkpoint_dir=None, max_checkpoints=1):
        self._data = {}
        self.checkpoint_dir = checkpoint_dir
        self.max_checkpoints = max_checkpoints

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
            self.cleanup_old_checkpoints()
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
        
    def add_figure(self, module_name, title=None, description=None, image_path=None, caption=None):
        """
        Adds figure to data context to ultimately be compiled into a report by ReportGenerator Module at end of pipeline run.

        Args:
            module_name (str): Module name that generated the figure.
            title (str, optional): Title of the figure. Defaults to None.
            description (str, optional): Description of the figure to be displayed below the title. Defaults to None.
            image_path (str, optional): Path of image to be displayed in the figure. Defaults to None.
            caption (str, optional): Caption to be displayed under the image. Defaults to None.
        """
        
        if 'REPORT_FIGURES' not in self._data:
            self._data['REPORT_FIGURES'] = {}

        if module_name not in self._data['REPORT_FIGURES']:
            self._data['REPORT_FIGURES'][module_name] = []

        self._data['REPORT_FIGURES'][module_name].append({
            'title': title,
            'description': description,
            'image_path': image_path,
            'caption': caption
        })

    def cleanup_old_checkpoints(self):
        """Remove old checkpoints if we have more than max_checkpoints"""
        if not self.checkpoint_dir:
            return
            
        checkpoints = []
        for filename in os.listdir(self.checkpoint_dir):
            if filename.endswith('.pkl'):
                checkpoint_path = os.path.join(self.checkpoint_dir, filename)
                created_time = os.path.getctime(checkpoint_path)
                checkpoints.append((checkpoint_path, created_time))
        
        # Sort by creation time (oldest first)
        checkpoints.sort(key=lambda x: x[1])
        
        # Remove oldest checkpoints if we have too many
        while len(checkpoints) > self.max_checkpoints:
            oldest_path, _ = checkpoints.pop(0)
            try:
                os.remove(oldest_path)
                self.logger.info(f"Removed old checkpoint: {oldest_path}")
            except Exception as e:
                self.logger.warning(f"Failed to remove old checkpoint {oldest_path}: {e}")