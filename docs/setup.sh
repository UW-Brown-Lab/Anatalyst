#!/bin/bash
# Master setup script for Read the Docs build

# Make the script exit on any error
set -e

echo "Setting up environment for documentation build..."

# Make R dependencies setup script executable
chmod +x setup_r_deps.sh

# Run the R dependencies setup if R is available
if command -v R &> /dev/null; then
    echo "Installing R dependencies..."
    # Create .Renviron file to set R library path
    echo "R_LIBS_USER='$HOME/R/library'" > $HOME/.Renviron
    ./setup_r_deps.sh
else
    echo "R is not installed, skipping R dependencies."
    echo "Documentation will use mock modules for R-dependent components."
fi

echo "Setup completed successfully!"