#!/bin/bash
# OVITO Installation Script for MLD Visualizations

echo "🎨 Installing OVITO for MLD Visualizations"
echo "=========================================="

# Method 1: Try conda-forge first (most reliable)
echo "🔧 Attempting conda-forge installation..."
conda remove shiboken6 pyside6 -y 2>/dev/null
conda install -c conda-forge ovito -y

# Test installation
if python -c "import ovito; print(f'✅ OVITO {ovito.__version__} installed successfully')" 2>/dev/null; then
    echo "🎉 OVITO installation successful!"
    exit 0
fi

echo "⚠️  Conda installation failed, trying pip method..."

# Method 2: Pip with force reinstall
pip uninstall ovito pyside6 shiboken6 -y 2>/dev/null
pip install --force-reinstall --no-deps ovito

# Test again
if python -c "import ovito; print(f'✅ OVITO {ovito.__version__} installed successfully')" 2>/dev/null; then
    echo "🎉 OVITO installation successful!"
    exit 0
fi

echo "❌ OVITO installation failed"
echo "💡 Alternative visualization tools available:"
echo "   - ASE GUI: pip install ase[gui]"
echo "   - NGLView: pip install nglview"  
echo "   - Py3DMol: pip install py3dmol"
echo "   - Plotly: pip install plotly"

echo ""
echo "🧪 Testing current visualization capabilities..."
python visualization_tools.py --check