#!/bin/bash
# System dependencies installation script
# Run this with sudo privileges

echo "=== Installing System Dependencies for MLD Simulation ==="
echo "This script requires sudo privileges."
echo ""

# Check if running as root or with sudo
if [[ $EUID -ne 0 ]]; then
   echo "This script must be run with sudo privileges"
   echo "Usage: sudo bash install_system_deps.sh"
   exit 1
fi

# Update package manager
echo "Updating package manager..."
apt update

# Essential development tools
echo "Installing development tools..."
apt install -y build-essential cmake git wget curl
apt install -y gcc g++ gfortran
apt install -y pkg-config

# MPI and parallel computing
echo "Installing MPI and parallel libraries..."
apt install -y libopenmpi-dev openmpi-bin
apt install -y libblas-dev liblapack-dev
apt install -y libfftw3-dev libfftw3-mpi-dev

# Python development
echo "Installing Python development tools..."
apt install -y python3-dev python3-pip
apt install -y python3-numpy-dev python3-scipy-dev

# Scientific libraries
echo "Installing scientific computing libraries..."
apt install -y libhdf5-dev libhdf5-mpi-dev
apt install -y libboost-all-dev
apt install -y libeigen3-dev

# GPAW specific dependencies
echo "Installing GPAW dependencies..."
apt install -y libxc-dev || echo "libxc-dev not available, will install via pip"
apt install -y libscalapack-openmpi-dev || echo "ScaLAPACK not available"

# Optional visualization tools
echo "Installing visualization tools..."
apt install -y povray imagemagick gnuplot || echo "Some visualization tools not available"

# Check installation
echo ""
echo "=== Installation Summary ==="
echo "Checking key components..."

if command -v gcc &> /dev/null; then
    echo "✓ GCC: $(gcc --version | head -n1)"
else
    echo "✗ GCC not found"
fi

if command -v mpirun &> /dev/null; then
    echo "✓ MPI: $(mpirun --version | head -n1)"
else
    echo "✗ MPI not found"
fi

if pkg-config --exists blas; then
    echo "✓ BLAS library found"
else
    echo "✗ BLAS library not found"
fi

if pkg-config --exists lapack; then
    echo "✓ LAPACK library found"
else
    echo "✗ LAPACK library not found"
fi

if pkg-config --exists fftw3; then
    echo "✓ FFTW3 library found"
else
    echo "✗ FFTW3 library not found"
fi

echo ""
echo "System dependencies installation complete!"
echo ""
echo "Next steps:"
echo "1. Exit sudo mode"
echo "2. Run: bash quick_install.sh"
echo "3. Or follow the detailed INSTALLATION_GUIDE.md"