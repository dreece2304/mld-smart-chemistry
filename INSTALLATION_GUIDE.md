# Complete Installation Guide for MLD Simulation Environment

## Overview
This guide will help you install all necessary components for the MLD simulation project, including system dependencies, DFT codes, and Python packages.

## Step 1: System Dependencies (Ubuntu/Debian)

First, update your package manager and install basic development tools:

```bash
sudo apt update && sudo apt upgrade -y

# Essential build tools
sudo apt install -y build-essential cmake git wget curl
sudo apt install -y gcc g++ gfortran
sudo apt install -y libopenmpi-dev openmpi-bin
sudo apt install -y libblas-dev liblapack-dev libscalapack-mpi-dev
sudo apt install -y libfftw3-dev libfftw3-mpi-dev
sudo apt install -y python3-dev python3-pip
sudo apt install -y pkg-config

# Additional libraries for scientific computing
sudo apt install -y libhdf5-dev libhdf5-mpi-dev
sudo apt install -y libboost-all-dev
sudo apt install -y libeigen3-dev
sudo apt install -y libxc-dev  # Exchange-correlation functionals
```

## Step 2: Install GPAW Dependencies

GPAW requires specific libraries that need to be installed system-wide:

```bash
# Install GPAW system dependencies
sudo apt install -y libxc-dev
sudo apt install -y libelpa-dev  # For faster linear algebra
sudo apt install -y libxcfun-dev

# For parallel calculations
sudo apt install -y libscalapack-openmpi-dev
sudo apt install -y libblacs-openmpi-dev

# Python development headers
sudo apt install -y python3-numpy-dev
sudo apt install -y python3-scipy-dev
```

## Step 3: Install Optional but Recommended Tools

```bash
# For visualization and analysis
sudo apt install -y povray  # For ASE visualization
sudo apt install -y imagemagick
sudo apt install -y gnuplot

# For molecular dynamics (if using LAMMPS)
sudo apt install -y libeigen3-dev libvoro++-dev

# For database and workflow management
sudo apt install -y mongodb  # For atomate2/fireworks
sudo apt install -y redis-server  # For job queuing
```

## Step 4: Verify Conda Environment

Check that your conda environment is properly set up:

```bash
conda activate mld_modeling
python --version  # Should show Python 3.10.x
which python     # Should point to conda environment
```

## Step 5: Install Python Packages (Core)

Install in this specific order to avoid conflicts:

```bash
# First, basic scientific packages
pip install numpy scipy matplotlib pandas

# ASE (Atomic Simulation Environment)
pip install ase

# Install GPAW with all dependencies
pip install gpaw
gpaw install-data  # Install pseudopotentials and basis sets
```

## Step 6: Install Structure Analysis Packages

```bash
# Crystal structure analysis
pip install pymatgen
pip install spglib
pip install seekpath

# Phonon calculations
pip install phonopy
pip install phono3py
```

## Step 7: Install Radiation Modeling Packages

```bash
# Electronic structure for radiation effects
pip install pyscf

# Monte Carlo for defect kinetics
pip install kmc-python

# Molecular dynamics
pip install mdanalysis
```

## Step 8: Install Visualization and Analysis Tools

```bash
# 3D visualization
pip install ovito-pro  # or ovito for free version
pip install nglview
pip install py3dmol

# Machine learning
pip install scikit-learn
pip install tensorflow
pip install matminer

# Plotting and analysis
pip install seaborn
pip install plotly
pip install jupyter
```

## Step 9: Install Optional Advanced Packages

```bash
# Workflow management
pip install fireworks
pip install atomate2

# Database tools
pip install pymongo
pip install gridfs

# Parallel computing
pip install dask
pip install joblib

# Additional DFT interfaces
pip install pymatgen-analysis-defects
pip install sumo  # Electronic structure plotting
```

## Step 10: Configure GPAW

Create GPAW configuration file:

```bash
mkdir -p ~/.gpaw
cat > ~/.gpaw/config.py << 'EOF'
# GPAW configuration
compiler = 'gcc'
mpicompiler = 'mpicc'
mpilinker = 'mpicc'
libraries = ['blas', 'lapack', 'fftw3']
library_dirs = ['/usr/lib/x86_64-linux-gnu']
include_dirs = ['/usr/include']

# Enable ScaLAPACK for parallel calculations
scalapack = True
define_macros = [('GPAW_NO_UNDERSCORE_CBLACS', '1'),
                 ('GPAW_NO_UNDERSCORE_CSCALAPACK', '1')]

# FFTW settings
fftw = True

# XC library
libxc = True
EOF
```

## Step 11: Test Installation

Run these tests to verify everything works:

```bash
# Test basic packages
python -c "import numpy, scipy, matplotlib, pandas; print('Basic packages: OK')"

# Test ASE
python -c "import ase; print('ASE version:', ase.__version__)"

# Test GPAW
python -c "import gpaw; print('GPAW version:', gpaw.__version__)"

# Test structure analysis
python -c "import pymatgen, spglib; print('Structure analysis: OK')"

# Test our project
cd /home/dreece23/MLD_ASE_GPAW
python test_workflow.py
```

## Step 12: Install LAMMPS (Optional for MD)

If you want molecular dynamics capabilities:

```bash
# Download and compile LAMMPS
cd /tmp
git clone -b stable https://github.com/lammps/lammps.git lammps
cd lammps
mkdir build && cd build

cmake ../cmake -DBUILD_SHARED_LIBS=yes -DPKG_PYTHON=yes
make -j4
sudo make install

# Install Python interface
cd ../python
python setup.py install
```

## Step 13: Final Configuration

Set up environment variables in your `.bashrc`:

```bash
echo 'export GPAW_SETUP_PATH=$HOME/.local/share/gpaw' >> ~/.bashrc
echo 'export OMP_NUM_THREADS=1' >> ~/.bashrc  # Important for MPI
echo 'export OMPI_MCA_btl_vader_single_copy_mechanism=none' >> ~/.bashrc
source ~/.bashrc
```

## Troubleshooting Common Issues

### If GPAW installation fails:
```bash
# Try installing with conda instead
conda activate mld_modeling
conda install -c conda-forge gpaw
```

### If you get MPI errors:
```bash
# Reinstall OpenMPI
sudo apt remove --purge openmpi-bin openmpi-common libopenmpi-dev
sudo apt install openmpi-bin openmpi-common libopenmpi-dev
```

### If FFTW is not found:
```bash
sudo apt install libfftw3-dev libfftw3-mpi-dev
# Then reinstall GPAW
pip uninstall gpaw
pip install gpaw --no-cache-dir
```

## Verification Checklist

After installation, verify these work:

- [ ] `python -c "import ase, gpaw, pymatgen"`
- [ ] `gpaw test` (runs GPAW test suite)
- [ ] `python test_workflow.py` (our project test)
- [ ] MPI test: `mpirun -np 2 python -c "import gpaw; print('MPI OK')"`

## Estimated Installation Time

- System packages: 15-30 minutes
- Python packages: 30-60 minutes  
- GPAW compilation: 15-30 minutes
- Total: 1-2 hours depending on internet speed

## Storage Requirements

- Base installation: ~2 GB
- GPAW datasets: ~1 GB
- Example calculations: ~500 MB
- Total: ~4 GB

## Next Steps

Once everything is installed:

1. Run `python test_workflow.py` to verify the setup
2. Try the basic MLD simulation: `python calculations/mld_deposition_workflow.py`
3. Explore the analysis tools: `python analysis/bulk_structure_analysis.py`

## Support

If you encounter issues:
1. Check the error messages carefully
2. Verify all system dependencies are installed
3. Try installing problematic packages individually
4. Check GPAW and ASE documentation for specific errors

## Alternative Installation Methods

If the above doesn't work, try these alternatives:

### Using Docker:
```bash
# Pull pre-built container with everything installed
docker pull materialsproject/pymatgen
```

### Using Singularity:
```bash
# For HPC environments
singularity pull docker://materialsproject/pymatgen
```

### Using pre-compiled binaries:
```bash
# Download pre-compiled GPAW from official site
wget https://wiki.fysik.dtu.dk/gpaw/install.html
```