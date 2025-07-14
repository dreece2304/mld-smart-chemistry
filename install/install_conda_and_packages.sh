#!/bin/bash
# MLD Conda and Package Installer (No sudo required)
# Run AFTER completing manual system dependency installation

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

# Configuration
CONDA_ENV_NAME="mld_modeling"
PYTHON_VERSION="3.10"
INSTALL_LOG="mld_conda_install_$(date +%Y%m%d_%H%M%S).log"

print_status() {
    case $2 in
        "pass") echo -e "${GREEN}✅ $1${NC}" | tee -a "$INSTALL_LOG" ;;
        "fail") echo -e "${RED}❌ $1${NC}" | tee -a "$INSTALL_LOG" ;;
        "warn") echo -e "${YELLOW}⚠️  $1${NC}" | tee -a "$INSTALL_LOG" ;;
        "info") echo -e "${BLUE}ℹ️  $1${NC}" | tee -a "$INSTALL_LOG" ;;
        "step") echo -e "${CYAN}🔧 $1${NC}" | tee -a "$INSTALL_LOG" ;;
        *) echo "$1" | tee -a "$INSTALL_LOG" ;;
    esac
}

print_header() {
    echo -e "${CYAN}" | tee -a "$INSTALL_LOG"
    echo "================================================================" | tee -a "$INSTALL_LOG"
    echo "$1" | tee -a "$INSTALL_LOG"
    echo "================================================================" | tee -a "$INSTALL_LOG"
    echo -e "${NC}" | tee -a "$INSTALL_LOG"
}

check_prerequisites() {
    print_header "CHECKING PREREQUISITES"
    
    # Check if system dependencies are installed
    missing_deps=()
    
    for cmd in gcc gfortran cmake mpirun; do
        if ! command -v "$cmd" &> /dev/null; then
            missing_deps+=("$cmd")
        fi
    done
    
    if [ ${#missing_deps[@]} -ne 0 ]; then
        print_status "Missing system dependencies: ${missing_deps[*]}" "fail"
        print_status "Please run manual installation steps first:" "info"
        print_status "See install/INSTALLATION_STEPS.md for details" "info"
        exit 1
    fi
    
    print_status "System dependencies: All found" "pass"
    
    # Check available space
    available_space=$(df . | awk 'NR==2 {print $4}')
    available_gb=$((available_space / 1024 / 1024))
    
    if [ "$available_gb" -lt 5 ]; then
        print_status "Insufficient disk space: ${available_gb}GB available, 5GB+ required" "fail"
        exit 1
    else
        print_status "Disk space: ${available_gb}GB available" "pass"
    fi
}

install_miniconda() {
    print_header "INSTALLING MINICONDA"
    
    # Check if conda already exists
    if command -v conda &> /dev/null; then
        print_status "Conda already installed: $(conda --version)" "info"
        
        # Source conda for this script
        CONDA_BASE=$(conda info --base)
        source "$CONDA_BASE/etc/profile.d/conda.sh"
        return 0
    fi
    
    print_status "Downloading Miniconda installer..." "step"
    cd /tmp
    wget -q --show-progress https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    
    print_status "Installing Miniconda..." "step"
    bash miniconda.sh -b -p "$HOME/miniconda3"
    
    # Initialize conda
    print_status "Initializing conda..." "step"
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    
    # Configure conda
    conda config --set auto_activate_base false
    conda config --set solver libmamba  # Faster solver
    
    # Add to bashrc if not already there
    if ! grep -q "miniconda3/etc/profile.d/conda.sh" ~/.bashrc; then
        echo '' >> ~/.bashrc
        echo '# Initialize conda' >> ~/.bashrc
        echo 'source "$HOME/miniconda3/etc/profile.d/conda.sh"' >> ~/.bashrc
    fi
    
    print_status "Miniconda installation completed" "pass"
    
    # Cleanup
    rm -f /tmp/miniconda.sh
}

create_mld_environment() {
    print_header "CREATING MLD CONDA ENVIRONMENT"
    
    # Source conda
    if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    else
        print_status "Conda not found. Please run miniconda installation first." "fail"
        exit 1
    fi
    
    # Remove existing environment if it exists
    if conda env list | grep -q "$CONDA_ENV_NAME"; then
        print_status "Removing existing environment: $CONDA_ENV_NAME" "step"
        conda env remove -n "$CONDA_ENV_NAME" -y
    fi
    
    # Create new environment
    print_status "Creating conda environment: $CONDA_ENV_NAME" "step"
    conda create -n "$CONDA_ENV_NAME" python="$PYTHON_VERSION" -y | tee -a "$INSTALL_LOG"
    
    # Activate environment
    conda activate "$CONDA_ENV_NAME"
    
    print_status "Environment created and activated" "pass"
}

install_scientific_packages() {
    print_header "INSTALLING SCIENTIFIC PACKAGES"
    
    # Ensure environment is activated
    conda activate "$CONDA_ENV_NAME"
    
    # Install core scientific packages via conda-forge (more reliable)
    print_status "Installing core scientific packages..." "step"
    conda install -c conda-forge -y \
        numpy \
        scipy \
        matplotlib \
        tqdm \
        psutil \
        requests | tee -a "$INSTALL_LOG"
    
    # Install MPI packages
    print_status "Installing MPI packages..." "step"
    conda install -c conda-forge -y \
        openmpi \
        mpi4py | tee -a "$INSTALL_LOG"
    
    print_status "Core packages installed" "pass"
}

install_ase_gpaw() {
    print_header "INSTALLING ASE AND GPAW"
    
    conda activate "$CONDA_ENV_NAME"
    
    # Install ASE first
    print_status "Installing ASE (Atomic Simulation Environment)..." "step"
    conda install -c conda-forge ase -y | tee -a "$INSTALL_LOG"
    
    # Install GPAW
    print_status "Installing GPAW (DFT calculator)..." "step"
    conda install -c conda-forge gpaw -y | tee -a "$INSTALL_LOG"
    
    print_status "ASE and GPAW installed successfully" "pass"
}

install_visualization_packages() {
    print_header "INSTALLING VISUALIZATION PACKAGES"
    
    conda activate "$CONDA_ENV_NAME"
    
    # Install via pip for better compatibility
    print_status "Installing visualization packages..." "step"
    pip install \
        plotly \
        py3dmol \
        nglview | tee -a "$INSTALL_LOG"
    
    # Note: Skipping OVITO due to Qt conflicts in headless environments
    print_status "Visualization packages installed (headless-compatible)" "pass"
}

download_gpaw_datasets() {
    print_header "DOWNLOADING GPAW PAW DATASETS"
    
    conda activate "$CONDA_ENV_NAME"
    
    # Determine GPAW setup directory
    GPAW_SETUP_DIR="$HOME/miniconda3/envs/$CONDA_ENV_NAME/share/gpaw"
    
    print_status "Downloading GPAW PAW datasets..." "step"
    print_status "This may take 5-15 minutes depending on connection..." "info"
    
    # Try automatic download first
    if gpaw install-data "$GPAW_SETUP_DIR" 2>&1 | tee -a "$INSTALL_LOG"; then
        print_status "GPAW datasets downloaded successfully" "pass"
    else
        print_status "Automatic download failed, trying manual method..." "warn"
        
        # Manual download as fallback
        mkdir -p "$GPAW_SETUP_DIR"
        cd "$GPAW_SETUP_DIR"
        
        # Try latest version
        if wget -q --show-progress https://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-24.1.0.tar.gz; then
            tar -xzf gpaw-setups-24.1.0.tar.gz --strip-components=1
            rm -f gpaw-setups-24.1.0.tar.gz
            print_status "GPAW datasets installed manually" "pass"
        else
            print_status "Manual download also failed - datasets may need manual installation" "fail"
            print_status "Try: gpaw install-data $GPAW_SETUP_DIR" "info"
        fi
    fi
    
    # Set environment variable permanently
    if ! grep -q "GPAW_SETUP_PATH" ~/.bashrc; then
        echo '' >> ~/.bashrc
        echo '# GPAW setup path' >> ~/.bashrc
        echo "export GPAW_SETUP_PATH=\"$GPAW_SETUP_DIR\"" >> ~/.bashrc
    fi
    
    export GPAW_SETUP_PATH="$GPAW_SETUP_DIR"
    
    # Count datasets
    if [ -d "$GPAW_SETUP_DIR" ]; then
        dataset_count=$(ls "$GPAW_SETUP_DIR"/*.gz 2>/dev/null | wc -l)
        print_status "GPAW datasets: $dataset_count files available" "info"
    fi
}

configure_mld_environment() {
    print_header "CONFIGURING MLD ENVIRONMENT"
    
    # WSL-specific optimizations
    if grep -q Microsoft /proc/version; then
        print_status "Applying WSL2 optimizations..." "step"
        
        # Add WSL configuration to bashrc
        if ! grep -q "# MLD WSL Configuration" ~/.bashrc; then
            cat >> ~/.bashrc << 'EOF'

# MLD WSL Configuration
if grep -q Microsoft /proc/version; then
    # Force headless backends for compatibility
    export MPLBACKEND=Agg
    export QT_QPA_PLATFORM=offscreen
    export DISPLAY=""
    
    # Graphics optimizations
    export LIBGL_ALWAYS_INDIRECT=1
    export MESA_GL_VERSION_OVERRIDE=3.3
    
    # MPI optimizations for WSL2
    export OMPI_MCA_btl_vader_single_copy_mechanism=none
fi
EOF
        fi
        print_status "WSL2 optimizations added" "pass"
    fi
    
    # Set optimal thread count
    cpu_cores=$(nproc)
    optimal_threads=$((cpu_cores > 8 ? 8 : cpu_cores))
    
    if ! grep -q "OMP_NUM_THREADS" ~/.bashrc; then
        echo "export OMP_NUM_THREADS=$optimal_threads" >> ~/.bashrc
    fi
    
    print_status "Environment configuration completed" "pass"
}

install_mld_package() {
    print_header "INSTALLING MLD FRAMEWORK"
    
    conda activate "$CONDA_ENV_NAME"
    
    # Install MLD package in development mode
    if [ -f "setup.py" ]; then
        print_status "Installing MLD framework in development mode..." "step"
        pip install -e . | tee -a "$INSTALL_LOG"
        print_status "MLD framework installed" "pass"
    else
        print_status "No setup.py found - MLD modules will be imported directly" "info"
        
        # Add src to Python path in bashrc
        project_dir=$(pwd)
        if ! grep -q "PYTHONPATH.*mld-smart-chemistry" ~/.bashrc; then
            echo '' >> ~/.bashrc
            echo '# MLD framework Python path' >> ~/.bashrc
            echo "export PYTHONPATH=\"$project_dir/src:\$PYTHONPATH\"" >> ~/.bashrc
        fi
        print_status "MLD framework path configured" "pass"
    fi
}

run_basic_validation() {
    print_header "RUNNING BASIC VALIDATION"
    
    conda activate "$CONDA_ENV_NAME"
    
    print_status "Testing basic imports..." "step"
    
    # Test core packages
    python -c "
import numpy, scipy, matplotlib, ase, gpaw, tqdm
print('✅ Core packages: OK')
" | tee -a "$INSTALL_LOG"
    
    # Test MPI
    python -c "
from mpi4py import MPI
print(f'✅ MPI: {MPI.COMM_WORLD.Get_size()} process(es)')
" | tee -a "$INSTALL_LOG"
    
    # Test GPAW datasets
    python -c "
import os
from pathlib import Path
setup_path = os.environ.get('GPAW_SETUP_PATH', '')
if setup_path and Path(setup_path).exists():
    datasets = list(Path(setup_path).glob('*.gz'))
    print(f'✅ GPAW datasets: {len(datasets)} files')
else:
    print('⚠️  GPAW datasets: Path not set or empty')
" | tee -a "$INSTALL_LOG"
    
    print_status "Basic validation completed" "pass"
}

print_completion_summary() {
    print_header "INSTALLATION COMPLETED"
    
    print_status "MLD framework conda installation successful!" "pass"
    print_status "Installation log: $INSTALL_LOG" "info"
    
    echo "" | tee -a "$INSTALL_LOG"
    echo "🎉 Next steps:" | tee -a "$INSTALL_LOG"
    echo "1. Restart terminal or run: source ~/.bashrc" | tee -a "$INSTALL_LOG"
    echo "2. Activate environment: conda activate $CONDA_ENV_NAME" | tee -a "$INSTALL_LOG"
    echo "3. Full validation: python install/validate_installation.py" | tee -a "$INSTALL_LOG"
    echo "4. Test functionality: python scripts/step1_create_h2o.py" | tee -a "$INSTALL_LOG"
    echo "" | tee -a "$INSTALL_LOG"
    
    echo "📁 Environment details:" | tee -a "$INSTALL_LOG"
    echo "   Conda environment: $CONDA_ENV_NAME" | tee -a "$INSTALL_LOG"
    echo "   Python version: $PYTHON_VERSION" | tee -a "$INSTALL_LOG"
    echo "   GPAW setup path: $GPAW_SETUP_PATH" | tee -a "$INSTALL_LOG"
    echo "" | tee -a "$INSTALL_LOG"
    
    echo "⚠️  If validation fails:" | tee -a "$INSTALL_LOG"
    echo "   - Check that system dependencies were installed (see INSTALLATION_STEPS.md)" | tee -a "$INSTALL_LOG"
    echo "   - Try: gpaw install-data \$GPAW_SETUP_PATH" | tee -a "$INSTALL_LOG"
    echo "   - Run: python install/validate_installation.py --verbose" | tee -a "$INSTALL_LOG"
}

main() {
    print_header "MLD CONDA & PACKAGES INSTALLER"
    
    echo "This script installs conda, creates the MLD environment, and installs all packages."
    echo "This script does NOT require sudo privileges."
    echo ""
    echo "Prerequisites (must be done manually first):"
    echo "  - System dependencies installed (build-essential, MPI, etc.)"
    echo "  - See install/INSTALLATION_STEPS.md for manual steps"
    echo ""
    echo "Continue? (Y/n)"
    read -r response
    if [[ "$response" =~ ^[Nn]$ ]]; then
        echo "Installation cancelled"
        exit 0
    fi
    
    # Run installation steps
    check_prerequisites
    install_miniconda
    create_mld_environment
    install_scientific_packages
    install_ase_gpaw
    install_visualization_packages
    download_gpaw_datasets
    configure_mld_environment
    install_mld_package
    run_basic_validation
    print_completion_summary
    
    print_status "Total installation time: $SECONDS seconds" "info"
}

# Trap errors
trap 'echo "❌ Installation failed at line $LINENO. Check $INSTALL_LOG for details."' ERR

# Run main installation
main "$@"