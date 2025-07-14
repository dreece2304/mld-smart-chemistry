#!/bin/bash
# Quick installation script for MLD simulation environment
# Run this after completing the system dependencies installation

set -e  # Exit on any error

echo "=== MLD Simulation Environment Quick Install ==="
echo "This script installs Python packages only."
echo "Make sure you've installed system dependencies first!"
echo ""

# Check if conda environment exists
if ! conda info --envs | grep -q "mld_modeling"; then
    echo "Creating conda environment..."
    conda create -n mld_modeling python=3.10 -y
fi

# Activate environment
echo "Activating conda environment..."
eval "$(conda shell.bash hook)"
conda activate mld_modeling

# Verify environment
echo "Python version: $(python --version)"
echo "Environment: $(which python)"

# Install packages in order
echo ""
echo "Installing core scientific packages..."
pip install --upgrade pip
pip install numpy scipy matplotlib pandas

echo ""
echo "Installing ASE..."
pip install ase

echo ""
echo "Installing GPAW (this may take a while)..."
pip install gpaw

echo ""
echo "Installing GPAW datasets..."
gpaw install-data

echo ""
echo "Installing structure analysis packages..."
pip install pymatgen spglib seekpath phonopy

echo ""
echo "Installing additional analysis tools..."
pip install scikit-learn seaborn plotly jupyter

echo ""
echo "Installing optional packages (some may fail, that's OK)..."
pip install pyscf || echo "PySCF installation failed (optional)"
pip install ovito || echo "Ovito installation failed (optional)"
pip install nglview || echo "NGLView installation failed (optional)"
pip install py3dmol || echo "Py3DMol installation failed (optional)"
pip install mdanalysis || echo "MDAnalysis installation failed (optional)"

echo ""
echo "Configuring GPAW..."
mkdir -p ~/.gpaw
cat > ~/.gpaw/config.py << 'EOF'
# Basic GPAW configuration
libraries = ['blas', 'lapack']
define_macros = []
EOF

echo ""
echo "Testing installation..."
python -c "import numpy, scipy, matplotlib, pandas; print('✓ Basic packages OK')"
python -c "import ase; print('✓ ASE version:', ase.__version__)"

if python -c "import gpaw; print('✓ GPAW version:', gpaw.__version__)" 2>/dev/null; then
    echo "✓ GPAW installed successfully"
else
    echo "⚠ GPAW installation may have issues"
fi

if python -c "import pymatgen; print('✓ PyMatGen version:', pymatgen.__version__)" 2>/dev/null; then
    echo "✓ PyMatGen installed successfully"
else
    echo "⚠ PyMatGen installation may have issues"
fi

echo ""
echo "Testing project workflow..."
cd "$(dirname "$0")"
if python test_workflow.py > /dev/null 2>&1; then
    echo "✓ Project workflow test passed"
else
    echo "⚠ Project workflow test failed (may need system dependencies)"
fi

echo ""
echo "=== Installation Summary ==="
echo "Conda environment: mld_modeling"
echo "To activate: conda activate mld_modeling"
echo "To test: python test_workflow.py"
echo ""
echo "If you see any warnings above, consult INSTALLATION_GUIDE.md"
echo "for detailed troubleshooting steps."