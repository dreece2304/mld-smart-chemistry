#!/bin/bash
# Automated MLD Environment Setup Script
# Fixes common environment issues and sets up permanent configuration

set -e  # Exit on any error

echo "🔧 MLD ENVIRONMENT SETUP"
echo "========================"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_status() {
    case $2 in
        "pass") echo -e "${GREEN}✅ $1${NC}" ;;
        "fail") echo -e "${RED}❌ $1${NC}" ;;
        "warn") echo -e "${YELLOW}⚠️  $1${NC}" ;;
        "info") echo -e "${BLUE}ℹ️  $1${NC}" ;;
        *) echo "$1" ;;
    esac
}

# Check if conda is available
if ! command -v conda &> /dev/null; then
    print_status "Conda not found! Please install Miniconda or Anaconda first." "fail"
    exit 1
fi

print_status "Conda found" "pass"

# Check if mld_modeling environment exists
if conda env list | grep -q "mld_modeling"; then
    print_status "mld_modeling environment exists" "pass"
else
    print_status "Creating mld_modeling environment..." "info"
    conda create -n mld_modeling python=3.10 -y
    print_status "Created mld_modeling environment" "pass"
fi

# Check current environment
if [[ "$CONDA_DEFAULT_ENV" == "mld_modeling" ]]; then
    print_status "Already in mld_modeling environment" "pass"
else
    print_status "Please run: conda activate mld_modeling" "warn"
    print_status "Then re-run this script" "info"
    exit 1
fi

# Set up GPAW_SETUP_PATH permanently
GPAW_PATH="$HOME/.local/share/gpaw/gpaw-setups-24.11.0"
echo "export GPAW_SETUP_PATH=\"$GPAW_PATH\"" >> ~/.bashrc

# Also set for current session
export GPAW_SETUP_PATH="$GPAW_PATH"

print_status "Set GPAW_SETUP_PATH permanently" "pass"

# Check if GPAW datasets exist
if [ -d "$GPAW_PATH" ]; then
    print_status "GPAW datasets found" "pass"
else
    print_status "Installing GPAW datasets..." "info"
    gpaw install-data ~/.local/share/gpaw
    if [ -d "$GPAW_PATH" ]; then
        print_status "GPAW datasets installed successfully" "pass"
    else
        print_status "Failed to install GPAW datasets" "fail"
        exit 1
    fi
fi

# Install/upgrade core packages
print_status "Installing core packages..." "info"

pip install --upgrade pip

# Core scientific packages
pip install "ase>=3.22.0"
pip install "gpaw>=22.8.0" 
pip install "numpy>=1.21.0"
pip install "scipy>=1.7.0"
pip install "matplotlib>=3.5.0"
pip install "tqdm>=4.0.0"

# Essential for MLD
pip install "pymatgen>=2023.1.0"
pip install "spglib>=2.0.0"

print_status "Core packages installed" "pass"

# Create directories if they don't exist
mkdir -p structures
mkdir -p results
mkdir -p logs

print_status "Created necessary directories" "pass"

# Test basic functionality
print_status "Testing basic functionality..." "info"

python -c "
import ase, gpaw, numpy, scipy, matplotlib, tqdm
from ase import Atoms
from gpaw import GPAW, PW
print('✅ All core packages import successfully')

# Test GPAW setup
import os
if os.path.exists(os.environ.get('GPAW_SETUP_PATH', '')):
    print('✅ GPAW datasets accessible')
else:
    print('❌ GPAW datasets not accessible')
    exit(1)

# Test basic ASE functionality
h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
h2.center(vacuum=5.0)
print('✅ ASE Atoms creation works')

# Test GPAW calculator creation
calc = GPAW(mode=PW(200), xc='PBE', kpts=(1,1,1), txt=None)
print('✅ GPAW calculator creation works')
"

if [ $? -eq 0 ]; then
    print_status "Basic functionality test passed" "pass"
else
    print_status "Basic functionality test failed" "fail"
    exit 1
fi

# Final summary
echo ""
echo "🎉 ENVIRONMENT SETUP COMPLETE!"
echo "=============================="
echo ""
print_status "Environment: mld_modeling" "pass"
print_status "GPAW datasets: installed" "pass"
print_status "Core packages: installed" "pass"
print_status "Basic tests: passed" "pass"
echo ""
echo "Next steps:"
echo "1. Run: python validate_setup.py --full"
echo "2. Run: python quick_run.py --test"
echo "3. Start MLD simulations!"
echo ""
print_status "Setup completed successfully!" "pass"