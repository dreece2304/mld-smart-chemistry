#!/bin/bash

# Complete installation script for MLD modeling with radiation effects
echo "=== MLD Modeling Environment Setup ==="
echo "Installing ASE + GPAW with bulk structure and radiation modeling capabilities"

# Check Python version
python_version=$(python3 --version 2>&1 | grep -oP '\d+\.\d+' | head -1)
echo "Python version: $python_version"

# Create conda environment
echo "Creating conda environment..."
conda create -n mld_modeling python=3.10 -y
conda activate mld_modeling

# Install core scientific packages
echo "Installing core packages..."
conda install -c conda-forge \
    numpy scipy matplotlib pandas \
    ase gpaw \
    pymatgen spglib phonopy \
    scikit-learn jupyter ipython \
    mpi4py openmpi \
    -y

# Install GPAW setups
echo "Installing GPAW pseudopotentials..."
gpaw install-data --register

# Install additional packages via pip
echo "Installing specialized packages..."
pip install -r requirements.txt

# Install radiation damage tools
echo "Installing radiation damage simulation tools..."

# LAMMPS for MD simulations
echo "Installing LAMMPS..."
conda install -c conda-forge lammps -y

# Install Monte Carlo packages for defect kinetics
pip install kmc-python pykmc

# Install optical property calculation tools
echo "Installing optical property tools..."
pip install exciting-python

# Install machine learning tools for structure prediction
echo "Installing ML tools..."
pip install tensorflow torch scikit-learn matminer

# Set up GPAW for parallel calculations
echo "Configuring GPAW for parallel execution..."
cat > ~/.gpaw/config.py << EOF
scalapack = True
fftw = True
libraries = ['blas', 'lapack', 'fftw3', 'scalapack-openmpi']
mpicompiler = 'mpicc'
mpilinker = 'mpicc'
EOF

# Verify installation
echo "Verifying installation..."
python3 -c "
import ase
import gpaw
import pymatgen
import numpy as np
print('✓ ASE version:', ase.__version__)
print('✓ GPAW version:', gpaw.__version__)
print('✓ PyMatGen version:', pymatgen.__version__)
print('✓ NumPy version:', np.__version__)
print('Installation successful!')
"

echo ""
echo "=== Next Steps ==="
echo "1. Activate environment: conda activate mld_modeling"
echo "2. Test GPAW: python3 -c 'from gpaw import GPAW; print(\"GPAW ready\")'"
echo "3. Run initial calculations: cd calculations && python3 optimize_structures.py"
echo ""
echo "For radiation damage modeling, additional setup may be needed:"
echo "- Configure cluster for large-scale MD simulations"
echo "- Set up Monte Carlo defect kinetics calculations"
echo "- Install specialized radiation transport codes if needed"