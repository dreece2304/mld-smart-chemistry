# Installation & Setup Guide

## Quick Start

### Prerequisites Check
```bash
# Verify conda environment
conda activate mld_modeling

# Test core functionality
python tests/test_suite.py --level 1
```

### Immediate Usage
```bash
# Quick functionality test
python examples/quick_run.py --test

# Basic MLD simulation
python src/simulation/dft_fast_mld.py
```

---

## Complete Installation

### 1. System Dependencies (Ubuntu/Debian)

```bash
# Update system
sudo apt update && sudo apt upgrade -y

# Essential build tools
sudo apt install -y build-essential cmake git wget curl
sudo apt install -y gcc g++ gfortran
sudo apt install -y python3-dev python3-pip

# Scientific computing libraries
sudo apt install -y libopenmpi-dev openmpi-bin
sudo apt install -y libblas-dev liblapack-dev libscalapack-mpi-dev
sudo apt install -y libfftw3-dev libfftw3-mpi-dev
sudo apt install -y libhdf5-dev libhdf5-mpi-dev
sudo apt install -y pkg-config

# GPAW dependencies
sudo apt install -y libxc-dev
sudo apt install -y libelpa-dev
sudo apt install -y libscalapack-openmpi-dev
sudo apt install -y libblacs-openmpi-dev

# Graphics libraries (for visualization)
sudo apt install -y libgl1-mesa-glx libglu1-mesa
```

### 2. Conda Environment Setup

```bash
# Create environment
conda create -n mld_modeling python=3.9

# Activate environment
conda activate mld_modeling

# Install core packages
conda install -y numpy scipy matplotlib
conda install -y -c conda-forge ase
conda install -y -c conda-forge gpaw
conda install -y -c conda-forge mpi4py

# Install additional packages
pip install tqdm rich plotly py3dmol
```

### 3. GPAW Configuration

```bash
# Download GPAW datasets
python -m gpaw install-data

# Set environment variable
export GPAW_SETUP_PATH=~/.local/share/gpaw/gpaw-setups-24.11.0

# Add to your shell profile
echo 'export GPAW_SETUP_PATH=~/.local/share/gpaw/gpaw-setups-24.11.0' >> ~/.bashrc
```

### 4. Validation

```bash
# Run comprehensive test suite
python tests/test_suite.py

# Test specific levels
python tests/test_suite.py --level 1-3

# Validate installation
python scripts/setup/validate_installation.py
```

---

## Optional Visualization Tools

### For Enhanced 3D Visualization

```bash
# In mld_modeling environment
conda activate mld_modeling

# Essential visualization tools
pip install py3dmol plotly nglview

# For Jupyter notebooks
pip install jupyter
jupyter nbextension enable --py nglview
```

### Professional OVITO (Advanced)

**Option 1 - Desktop Application (Recommended)**:
1. Download from: https://www.ovito.org/download/
2. Install OVITO Basic (free version)
3. Use for high-quality structure visualization

**Option 2 - Python Integration**:
```bash
# Try official OVITO conda channel
conda install -c https://conda.ovito.org ovito -y
```

---

## SSH Server Setup

### Connection Setup
```bash
# Basic connection
ssh -p 2222 dreece23@172.25.145.180

# With SSH config (~/.ssh/config)
Host lab-desktop
    HostName 172.25.145.180
    Port 2222
    User dreece23
    ForwardX11 yes
```

### Remote Execution
```bash
# Background execution
nohup python src/simulation/dft_small_mld.py > output.log 2>&1 &

# Monitor progress
python scripts/monitoring/monitor_dft_progress.py output.log

# Check status
tail -f output.log
```

### Connection Troubleshooting
```bash
# Test connection
ssh -v -p 2222 dreece23@172.25.145.180

# Check VPN connection (required)
ping 172.25.145.180

# Port forwarding check (lab desktop PowerShell as Admin)
netsh interface portproxy show all
```

---

## Environment Variables

### Required Settings
```bash
# GPAW setup path
export GPAW_SETUP_PATH=~/.local/share/gpaw/gpaw-setups-24.11.0

# MPI settings (if using parallel)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

### Optional Performance Settings
```bash
# For large calculations
export GPAW_PYTHON_EGG_CACHE=/tmp
export PYTHONPATH="/path/to/mld-smart-chemistry:$PYTHONPATH"
```

---

## Common Issues & Solutions

### GPAW Installation Problems
```bash
# Missing libxc
sudo apt install -y libxc-dev

# Rebuild GPAW
pip uninstall gpaw
pip install gpaw --no-cache-dir
```

### MPI Issues
```bash
# Test MPI
mpirun -np 2 python -c "from mpi4py import MPI; print(MPI.COMM_WORLD.Get_rank())"

# If MPI fails
sudo apt install --reinstall libopenmpi-dev openmpi-bin
```

### Memory Issues
```bash
# Monitor memory usage
python scripts/monitoring/monitor_dft_progress.py --memory

# Reduce memory usage
export GPAW_SETUP_PATH_PRIORITY=1
```

### Graphics Issues (WSL2)
```bash
# Test headless mode
python -c "import matplotlib; matplotlib.use('Agg'); print('Headless OK')"

# Install WSL2 graphics
sudo apt install -y mesa-utils
```

---

## Performance Optimization

### For Large Calculations
```bash
# Use more cores
export OMP_NUM_THREADS=4

# Parallel execution
mpirun -np 4 python src/simulation/dft_small_mld.py
```

### Memory Management
```bash
# Reduce memory usage
python src/simulation/dft_fast_mld.py  # Use fast version

# Monitor resources
htop  # Install with: sudo apt install htop
```

---

## Verification Commands

### Quick Tests
```bash
# Basic imports
python -c "import ase, gpaw; print('✅ Core packages OK')"

# Test calculation
python examples/quick_run.py --test

# Full validation
python tests/test_suite.py --level 1
```

### System Status
```bash
# Check conda environment
conda info --envs

# Check Python packages
pip list | grep -E "(ase|gpaw|numpy|scipy)"

# Check system resources
free -h
df -h
```

---

## Support

### Diagnostic Commands
```bash
# Generate system report
python scripts/setup/validate_installation.py --report

# Check installation logs
tail -f ~/.conda/envs/mld_modeling/conda-meta/history

# Test specific functionality
python tests/test_suite.py --level 1 --verbose
```

### Getting Help
1. Check [troubleshooting section](#common-issues--solutions)
2. Run diagnostic validation: `python scripts/setup/validate_installation.py`
3. Check system resources: `htop` and `free -h`
4. Review installation logs
5. Test with minimal example: `python examples/quick_run.py --test`

---

**Last Updated**: July 2025