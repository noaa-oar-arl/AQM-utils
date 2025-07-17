#!/bin/bash
# setup_environment.sh
# Script to set up Python virtual environment for GRIB2 processing

set -e  # Exit on any error

VENV_NAME="grib2_env"
PYTHON_VERSION="python3"

echo "Setting up Python virtual environment for GRIB2 processing..."

# Check if python3 is available
if ! command -v $PYTHON_VERSION &> /dev/null; then
    echo "Error: $PYTHON_VERSION is not available. Please install Python 3.7 or later."
    exit 1
fi

# Create virtual environment
echo "Creating virtual environment: $VENV_NAME"
$PYTHON_VERSION -m venv $VENV_NAME

# Activate virtual environment
echo "Activating virtual environment..."
source $VENV_NAME/bin/activate

# Upgrade pip
echo "Upgrading pip..."
pip install --upgrade pip

# Install requirements
echo "Installing requirements from requirements.txt..."
pip install -r requirements.txt

echo ""
echo "Virtual environment setup complete!"
echo ""
echo "To activate the environment in the future, run:"
echo "  source $VENV_NAME/bin/activate"
echo ""
echo "To test the GRIB2 reader, run:"
echo "  python grib2_to_cube.py"
echo ""
echo "To deactivate the environment when done:"
echo "  deactivate"
