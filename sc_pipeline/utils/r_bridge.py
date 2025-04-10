# sc_pipeline/utils/r_bridge.py

import os
import tempfile
import subprocess
import logging
import pickle

class RBridge:
    """Bridge for calling R functions from Python"""

    def __init__(self, r_script_dir=None):
        self.logger = logging.getLogger("RBridge")
        self.r_script_dir = r_script_dir

    def run_r_script(self, script_name, args=None, data=None):
        """
        Run an R script with optional arguments and data.
        
        Args:
            script_name: Name of the R script to run
            args: Dictionary of arguments to pass to the R script
            data: Dictionary of data to pass to the R script
            
        Returns:
            dict: Results from the R script
        """
        if self.r_script_dir:
            script_path = os.path.join(self.r_script_dir, script_name)
        else:
            # Try to find the script in the package
            import sc_pipeline
            package_dir = os.path.dirname(sc_pipeline.__file__)
            script_path = os.path.join(package_dir, 'r_scripts', script_name)

        if not os.path.exists(script_path):
            raise FileNotFoundError(f"R script not found: {script_path}")
        
        # Create temporary files for input and output
        with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as input_file, tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as output_file:
            input_path = input_file.name
            output_path = output_file.name

            # Prepare input data
            input_data = {
                'args': args or {},
                'data': data or {}
            }

            with open(input_path, "wb") as f:
                pickle.dump(input_data, f)
            
            # Build the R command
            r_cmd = [
                'Rscript',
                script_path,
                '--input', input_path,
                '--output', output_path
            ]

            # Run the R script
            try:
                self.logger.info(f"Running R script: {script_name}")
                result = subprocess.run(r_cmd, check=True, capture_output=True, text=True)

                # Check for errors
                if result.returncode != 0:
                    self.logger.error(f"R script failed: {result.stderr}")
                    raise RuntimeError(f"R script failed: {result.stderr}")
                
                # Load the results
                with open(output_path, 'rb') as f:
                    results = pickle.load(f)

                return results
            
            except subprocess.CalledProcessError as e:
                self.logger.error(f"Error running R script: {e.stderr}")
                raise

            finally:
                # Clean up temporary files
                for path in [input_path, output_path]:
                    if os.path.exists(path):
                        os.unlink(path)