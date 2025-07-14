#!/bin/bash
# Complete Clean MLD Environment Installation
# Optimized for WSL2 Ubuntu with proper dependency management

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

print_status() {
    case $2 in
        "pass") echo -e "${GREEN}✅ $1${NC}" ;;
        "fail") echo -e "${RED}❌ $1${NC}" ;;
        "warn") echo -e "${YELLOW}⚠️  $1${NC}" ;;
        "info") echo -e "${BLUE}ℹ️  $1${NC}" ;;
        "step") echo -e "${CYAN}🔧 $1${NC}" ;;
        *) echo "$1" ;;
    esac
}

print_header() {
    echo -e "${MAGENTA}"
    echo "==========================================="
    echo "$1"
    echo "==========================================="
    echo -e "${NC}"
}

# Check if we're in WSL
check_wsl() {
    if grep -q Microsoft /proc/version; then
        print_status "Running in WSL2 - optimizing for headless operation" "info"
        export DISPLAY=""
        export MPLBACKEND=Agg
        export QT_QPA_PLATFORM=offscreen
    else
        print_status "Running on native Linux" "info"
    fi
}

# Verify conda is available
check_conda() {
    if ! command -v conda &> /dev/null; then
        print_status "Conda not found! Please install Miniconda or Anaconda first." "fail"
        exit 1
    fi
    
    print_status "Conda found: $(conda --version)" "pass"
}

# Check system dependencies
check_system_deps() {
    print_header "CHECKING SYSTEM DEPENDENCIES"
    
    # Essential build tools
    missing_packages=()
    
    if ! dpkg -l | grep -q build-essential; then
        missing_packages+=("build-essential")
    fi
    
    if ! dpkg -l | grep -q gfortran; then
        missing_packages+=("gfortran")
    fi
    
    if ! dpkg -l | grep -q libopenmpi-dev; then
        missing_packages+=("libopenmpi-dev")
    fi
    
    if [ ${#missing_packages[@]} -gt 0 ]; then
        print_status "Missing system packages: ${missing_packages[*]}" "warn"
        print_status "Please install with: sudo apt install ${missing_packages[*]}" "info"
        echo "Continue anyway? (y/N)"
        read -r response
        if [[ ! "$response" =~ ^[Yy]$ ]]; then
            exit 1
        fi
    else
        print_status "System dependencies satisfied" "pass"
    fi
}

# Remove old environment
remove_old_env() {
    print_header "REMOVING OLD ENVIRONMENT"
    
    if conda env list | grep -q "mld_modeling"; then
        print_status "Removing existing mld_modeling environment..." "step"
        conda env remove -n mld_modeling --yes || true
        print_status "Old environment removed" "pass"
    else
        print_status "No existing mld_modeling environment found" "info"
    fi
}

# Create fresh environment
create_environment() {
    print_header "CREATING FRESH ENVIRONMENT"
    
    print_status "Creating mld_modeling environment with Python 3.10..." "step"
    conda create -n mld_modeling python=3.10 -y
    
    print_status "Environment created successfully" "pass"
}

# Activate environment function
activate_env() {
    print_status "Activating mld_modeling environment..." "step"
    eval "$(conda shell.bash hook)"
    conda activate mld_modeling
    
    # Verify activation
    if [[ "$CONDA_DEFAULT_ENV" == "mld_modeling" ]]; then
        print_status "Environment activated: $CONDA_DEFAULT_ENV" "pass"
    else
        print_status "Failed to activate environment" "fail"
        exit 1
    fi
}

# Install core scientific packages
install_core_packages() {
    print_header "INSTALLING CORE SCIENTIFIC PACKAGES"
    
    print_status "Installing MPI and parallel computing..." "step"
    conda install -c conda-forge openmpi mpi4py -y
    
    print_status "Installing NumPy, SciPy, Matplotlib..." "step"
    conda install -c conda-forge numpy scipy matplotlib -y
    
    print_status "Installing ASE and GPAW..." "step"
    conda install -c conda-forge ase gpaw -y
    
    print_status "Installing additional utilities..." "step"
    pip install tqdm psutil
    
    print_status "Core packages installed successfully" "pass"
}

# Install optional packages
install_optional_packages() {
    print_header "INSTALLING OPTIONAL PACKAGES"
    
    print_status "Installing structure analysis packages..." "step"
    # Use conda for these to avoid conflicts
    conda install -c conda-forge pymatgen spglib -y || {
        print_status "Some optional packages failed - continuing" "warn"
    }
    
    print_status "Installing visualization packages..." "step"
    pip install plotly py3dmol || {
        print_status "Some visualization packages failed - continuing" "warn"
    }
    
    # NGLView for Jupyter (optional)
    pip install nglview || {
        print_status "NGLView installation failed - continuing" "warn"
    }
    
    print_status "Optional packages installation complete" "pass"
}

# Setup GPAW datasets
setup_gpaw() {
    print_header "SETTING UP GPAW DATASETS"
    
    print_status "Installing GPAW pseudopotential datasets..." "step"
    gpaw install-data ~/.local/share/gpaw
    
    # Set environment variable permanently
    GPAW_PATH="$HOME/.local/share/gpaw/gpaw-setups-24.11.0"
    
    if [ -d "$GPAW_PATH" ]; then
        print_status "GPAW datasets installed successfully" "pass"
        
        # Add to bashrc if not already there
        if ! grep -q "GPAW_SETUP_PATH" ~/.bashrc; then
            echo "export GPAW_SETUP_PATH=\"$GPAW_PATH\"" >> ~/.bashrc
            print_status "Added GPAW_SETUP_PATH to ~/.bashrc" "pass"
        fi
        
        # Set for current session
        export GPAW_SETUP_PATH="$GPAW_PATH"
        print_status "GPAW_SETUP_PATH set: $GPAW_SETUP_PATH" "pass"
    else
        print_status "GPAW datasets installation failed" "fail"
        exit 1
    fi
}

# Configure environment for WSL2
configure_wsl() {
    print_header "CONFIGURING WSL2 ENVIRONMENT"
    
    # Set headless backends
    if grep -q Microsoft /proc/version; then
        # Add WSL-specific settings to bashrc
        if ! grep -q "WSL MLD Settings" ~/.bashrc; then
            cat >> ~/.bashrc << 'EOF'

# WSL MLD Settings
export MPLBACKEND=Agg
export QT_QPA_PLATFORM=offscreen
export DISPLAY=""
EOF
            print_status "Added WSL-specific environment variables" "pass"
        fi
        
        # Set for current session
        export MPLBACKEND=Agg
        export QT_QPA_PLATFORM=offscreen
        export DISPLAY=""
        
        print_status "WSL2 configuration complete" "pass"
    else
        print_status "Not running in WSL - skipping WSL configuration" "info"
    fi
}

# Test installation
test_installation() {
    print_header "TESTING INSTALLATION"
    
    print_status "Testing core package imports..." "step"
    
    # Test core packages
    python -c "
import sys
print(f'Python: {sys.version.split()[0]}')

packages = ['numpy', 'scipy', 'matplotlib', 'ase', 'gpaw', 'mpi4py', 'tqdm']
for pkg in packages:
    try:
        module = __import__(pkg)
        version = getattr(module, '__version__', 'unknown')
        print(f'✅ {pkg}: {version}')
    except ImportError as e:
        print(f'❌ {pkg}: {e}')
        exit(1)

print('🎉 All core packages imported successfully!')
" || {
        print_status "Package import test failed" "fail"
        exit 1
    }
    
    print_status "Testing GPAW setup..." "step"
    python -c "
import os
import gpaw
gpaw_path = os.environ.get('GPAW_SETUP_PATH', '')
if not gpaw_path:
    print('❌ GPAW_SETUP_PATH not set')
    exit(1)
if not os.path.exists(gpaw_path):
    print(f'❌ GPAW path does not exist: {gpaw_path}')
    exit(1)
print(f'✅ GPAW setup path verified: {gpaw_path}')
" || {
        print_status "GPAW setup test failed" "fail"
        exit 1
    }
    
    print_status "Testing MPI functionality..." "step"
    python -c "
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
print(f'✅ MPI working: size={size}, rank={rank}')
" || {
        print_status "MPI test failed" "fail"
        exit 1
    }
    
    print_status "All tests passed successfully!" "pass"
}

# Main installation function
main() {
    print_header "MLD CLEAN INSTALLATION"
    echo "This script will:"
    echo "1. Remove existing mld_modeling environment"
    echo "2. Create fresh environment with Python 3.10"
    echo "3. Install all required packages"
    echo "4. Set up GPAW datasets"
    echo "5. Configure for WSL2 (if applicable)"
    echo "6. Test installation"
    echo ""
    echo "Continue? (Y/n)"
    read -r response
    if [[ "$response" =~ ^[Nn]$ ]]; then
        echo "Installation cancelled"
        exit 0
    fi
    
    # Run installation steps
    check_wsl
    check_conda
    check_system_deps
    remove_old_env
    create_environment
    activate_env
    install_core_packages
    install_optional_packages
    setup_gpaw
    configure_wsl
    test_installation
    
    # Final summary
    print_header "INSTALLATION COMPLETE"
    print_status "MLD environment successfully installed!" "pass"
    print_status "Environment: mld_modeling" "pass"
    print_status "Python: $(python --version 2>&1)" "pass"
    print_status "Conda environment path: $CONDA_PREFIX" "pass"
    
    echo ""
    echo -e "${GREEN}🎉 SUCCESS! Your MLD environment is ready.${NC}"
    echo ""
    echo "Next steps:"
    echo "1. conda activate mld_modeling"
    echo "2. python validate_clean_install.py"
    echo "3. python quick_run.py --test"
    echo ""
    echo "To use in new terminal sessions:"
    echo "  conda activate mld_modeling"
}

# Run main function
main "$@"