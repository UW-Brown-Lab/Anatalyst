# sc_pipeline/utils/r_bridge.py

import os
import tempfile
import subprocess
import logging
import shutil

class RBridge:
    """Bridge for calling R functions from Python with a workspace directory for file exchange"""

    def __init__(self, r_script_dir=None, workspace_dir=None, cleanup=True):
        self.logger = logging.getLogger("RBridge")
        self.r_script_dir = r_script_dir
        self.cleanup = cleanup
        
        # Create workspace directory if not provided
        if workspace_dir:
            self.workspace_dir = workspace_dir
            os.makedirs(workspace_dir, exist_ok=True)
        else:
            self.workspace_dir = tempfile.mkdtemp(prefix="r_bridge_workspace_")
            self.logger.info(f"Created temporary workspace: {self.workspace_dir}")

    def run_r_script(self, script_name, args=None):
        """
        Run an R script with the specified arguments.
        
        Args:
            script_name: Name of the R script to run
            args: Dictionary of arguments to pass to the R script
            
        Returns:
            tuple: (success, stdout, stderr)
        """
        try:
            # Locate the R script
            if self.r_script_dir:
                script_path = os.path.join(self.r_script_dir, script_name)
            else:
                # Try to find the script in the package
                import sc_pipeline
                package_dir = os.path.dirname(sc_pipeline.__file__)
                script_path = os.path.join(package_dir, 'r_scripts', script_name)

            if not os.path.exists(script_path):
                raise FileNotFoundError(f"R script not found: {script_path}")
            
            # Build the R command
            r_cmd = ['Rscript', script_path]
            
            # Add workspace directory as first argument
            r_cmd.append(self.workspace_dir)
            
            # Add any additional arguments
            if args:
                for key, value in args.items():
                    r_cmd.append(f"--{key}")
                    r_cmd.append(str(value))
            
            # Run the R script
            self.logger.info(f"Running R script: {script_name}")
            self.logger.debug(f"Command: {' '.join(r_cmd)}")
            
            result = subprocess.run(r_cmd, check=True, capture_output=True, text=True)
            
            self.logger.info(f"R script {script_name} completed successfully")
            return (True, result.stdout, result.stderr)
        
        except FileNotFoundError as e:
            self.logger.error(f"Error locating R script: {e}")
            return (False, "", str(e))
        
        except subprocess.CalledProcessError as e:
            self.logger.error(f"R script execution failed: {e.stderr}")
            return (False, e.stdout, e.stderr)
            
        except Exception as e:
            self.logger.error(f"Unexpected error: {e}", exc_info=True)
            return (False, "", str(e))
    
    def get_workspace_path(self, filename):
        """Get the full path to a file in the workspace directory"""
        return os.path.join(self.workspace_dir, filename)
    
    def cleanup_workspace(self):
        """Clean up the workspace directory if temporary"""
        if self.cleanup and self.workspace_dir and os.path.exists(self.workspace_dir):
            self.logger.info(f"Cleaning up workspace: {self.workspace_dir}")
            shutil.rmtree(self.workspace_dir)
    
    def __del__(self):
        """Destructor to clean up the workspace when the bridge is deleted"""
        if hasattr(self, 'cleanup') and self.cleanup:
            self.cleanup_workspace()