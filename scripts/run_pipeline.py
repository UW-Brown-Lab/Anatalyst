# sc_pipeline/scripts/run_pipeline.py
import argparse
from ast import parse
import logging
import os
import sys

# TEMPORARY TO ALLOW RUNNING WHOLE PIPELINE
#   "nothing more permanent than a temporary solution"
# Add the project root directory to the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
#################################################


from sc_pipeline.core.executor import PipelineExecutor

def setup_logging(log_file=None, log_level=logging.INFO):
    """Set up logging configuration"""
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    # Configure root logger
    logging.basicConfig(
        level=log_level,
        format=log_format,
        handlers=[
            logging.StreamHandler(), # Log to console
        ]
    )

    # Add file handler if log file is specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter(log_format))
        logging.getLogger().addHandler(file_handler)


def main():
    """Main entry point for the pipeline"""
    parser = argparse.ArgumentParser(description="Wesley's Single-cell RNA-seq Analysis pipeline")
    parser.add_argument("--config", required=True, help="Path to the configuration file")
    parser.add_argument("--checkpoint", help="Start from a specified checkpoint")
    parser.add_argument("--log-file", help="Path to the log file")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG","INFO","WARNING","ERROR"], help="Logging level to display")

    args = parser.parse_args()

    # Set up logging
    log_level = getattr(logging, args.log_level)
    setup_logging(args.log_file, log_level)
    logger = logging.getLogger("main")
    logger.info("Starting pipeline execution")

    try:
        # Run the pipeline
        executor = PipelineExecutor(args.config)
        success = executor.run(start_from=args.checkpoint)

        if success:
            logger.info("Pipeline completed successfully")
            return 0
        else:
            logger.error("Pipeline failed")
            return 1
    
    except Exception as e:
        logger.error(f"Pipeline execution failed: {e}", exc_info=True)
        return 1
    
if __name__ == "__main__":
    sys.exit(main())