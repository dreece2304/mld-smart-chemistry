#!/bin/bash
# Complete MLD Environment Setup - Master Installer
# Handles WSL2/Ubuntu installation from scratch with comprehensive validation

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
INSTALL_LOG="mld_install_$(date +%Y%m%d_%H%M%S).log"

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
    
    # Check OS
    if [[ "$OSTYPE" != "linux-gnu"* ]]; then
        print_status "This installer requires Linux (WSL2 or native Ubuntu)" "fail"
        exit 1
    fi
    
    # Check Ubuntu version
    if [ -f /etc/os-release ]; then
        . /etc/os-release
        print_status "Detected: $NAME $VERSION_ID" "info"
        
        # Check if version is supported
        if [[ "$VERSION_ID" < "20.04" ]]; then
            print_status "Ubuntu 20.04+ required. Current: $VERSION_ID" "fail"
            exit 1
        fi
    else
        print_status "Could not detect Ubuntu version" "warn"
    fi
    
    # Check if WSL
    if grep -q Microsoft /proc/version; then
        print_status "Running in WSL2 - will apply WSL optimizations" "info"
        export IS_WSL=true
    else
        print_status "Running on native Linux" "info"
        export IS_WSL=false
    fi
    
    # Check available space
    available_space=$(df . | awk 'NR==2 {print $4}')
    available_gb=$((available_space / 1024 / 1024))
    
    if [ "$available_gb" -lt 10 ]; then
        print_status "Insufficient disk space: ${available_gb}GB available, 10GB+ required" "fail"
        exit 1
    else
        print_status "Disk space: ${available_gb}GB available" "pass"
    fi
    
    # Check memory
    memory_gb=$(free -g | awk 'NR==2{print $2}')
    if [ "$memory_gb" -lt 4 ]; then
        print_status "Low memory: ${memory_gb}GB detected, 8GB+ recommended" "warn"
    else
        print_status "Memory: ${memory_gb}GB available" "pass"
    fi
    
    print_status "Prerequisites check completed" "pass"
}

install_system_dependencies() {
    print_header "INSTALLING SYSTEM DEPENDENCIES"
    
    # Update package list
    print_status "Updating package repositories..." "step"
    sudo apt update | tee -a "$INSTALL_LOG"
    
    # Essential build tools
    print_status "Installing build essentials..." "step"
    sudo apt install -y \
        build-essential \
        gfortran \
        cmake \
        pkg-config \
        curl \
        wget \
        git \
        python3-dev \
        python3-pip | tee -a "$INSTALL_LOG"
    
    # MPI and parallel computing
    print_status "Installing MPI and parallel libraries..." "step"
    sudo apt install -y \
        libopenmpi-dev \
        openmpi-bin \
        libscalapack-openmpi-dev \
        mpi-default-dev | tee -a "$INSTALL_LOG"
    
    # Scientific computing libraries
    print_status "Installing scientific libraries..." "step"
    sudo apt install -y \
        libfftw3-dev \
        liblapack-dev \
        libblas-dev \
        libxc-dev \
        libhdf5-serial-dev | tee -a "$INSTALL_LOG"
    
    # Graphics libraries (version-aware)
    print_status "Installing graphics libraries..." "step"
    
    ubuntu_version=$(lsb_release -rs)
    if [[ $(echo "$ubuntu_version >= 22.04" | bc -l) -eq 1 ]]; then
        # Ubuntu 22.04+
        sudo apt install -y \
            libgl1-mesa-dev \
            libglu1-mesa-dev \
            libglx-mesa0 || print_status "Some graphics packages failed" "warn"
    else
        # Ubuntu 20.04
        sudo apt install -y \
            libgl1-mesa-glx \
            libglu1-mesa-dev || print_status "Some graphics packages failed" "warn"
    fi
    
    print_status "System dependencies installed" "pass"
}

install_conda() {
    print_header "INSTALLING MINICONDA"
    
    # Check if conda already exists
    if command -v conda &> /dev/null; then
        print_status "Conda already installed: $(conda --version)" "info"
        return 0
    fi
    
    print_status "Downloading Miniconda..." "step"
    cd /tmp
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    
    print_status "Installing Miniconda..." "step"
    bash miniconda.sh -b -p "$HOME/miniconda3"
    
    # Initialize conda
    print_status "Initializing conda..." "step"
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda config --set auto_activate_base false
    
    # Add to bashrc if not already there
    if ! grep -q "miniconda3/etc/profile.d/conda.sh" ~/.bashrc; then
        echo '# Initialize conda' >> ~/.bashrc
        echo 'source "$HOME/miniconda3/etc/profile.d/conda.sh"' >> ~/.bashrc
    fi
    
    print_status "Conda installation completed" "pass"
    
    # Cleanup
    rm -f /tmp/miniconda.sh
}

create_mld_environment() {
    print_header "CREATING MLD CONDA ENVIRONMENT"
    
    # Source conda
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    
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
    
    # Install conda packages first (more reliable)
    print_status "Installing conda packages..." "step"
    conda install -c conda-forge -y \
        numpy \
        scipy \
        matplotlib \
        openmpi \
        mpi4py \
        tqdm \
        psutil | tee -a "$INSTALL_LOG"
    
    # Install ASE
    print_status "Installing ASE..." "step"
    conda install -c conda-forge ase -y | tee -a "$INSTALL_LOG"
    
    # Install GPAW
    print_status "Installing GPAW..." "step"
    conda install -c conda-forge gpaw -y | tee -a "$INSTALL_LOG"
    
    # Install additional packages via pip
    print_status "Installing additional packages..." "step"
    pip install \
        plotly \
        py3dmol \
        nglview \
        requests | tee -a "$INSTALL_LOG"
    
    print_status "MLD environment created successfully" "pass"
}

download_gpaw_datasets() {
    print_header "DOWNLOADING GPAW DATASETS"
    
    # Activate environment
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV_NAME"
    
    # Set GPAW setup path
    GPAW_SETUP_DIR="$HOME/miniconda3/envs/$CONDA_ENV_NAME/share/gpaw"
    
    print_status "Downloading GPAW PAW datasets..." "step"
    print_status "This may take 5-10 minutes..." "info"
    
    # Try to install datasets
    if gpaw install-data "$GPAW_SETUP_DIR" 2>&1 | tee -a "$INSTALL_LOG"; then
        print_status "GPAW datasets downloaded successfully" "pass"
    else
        print_status "GPAW dataset download failed - trying alternative method" "warn"
        
        # Alternative download method
        mkdir -p "$GPAW_SETUP_DIR"
        cd "$GPAW_SETUP_DIR"
        
        if wget -q https://wiki.fysik.dtu.dk/gpaw-files/gpaw-setups-24.1.0.tar.gz; then
            tar -xzf gpaw-setups-24.1.0.tar.gz
            mv gpaw-setups-24.1.0/* .
            rm -rf gpaw-setups-24.1.0 gpaw-setups-24.1.0.tar.gz
            print_status "GPAW datasets installed via alternative method" "pass"
        else
            print_status "GPAW dataset download failed" "fail"
            return 1
        fi
    fi
    
    # Set environment variable
    echo "export GPAW_SETUP_PATH=\"$GPAW_SETUP_DIR\"" >> ~/.bashrc
    export GPAW_SETUP_PATH="$GPAW_SETUP_DIR"
    
    # Count datasets
    dataset_count=$(ls "$GPAW_SETUP_DIR"/*.gz 2>/dev/null | wc -l)
    print_status "GPAW datasets: $dataset_count files installed" "info"
}

configure_environment() {
    print_header "CONFIGURING ENVIRONMENT"
    
    # WSL-specific configuration
    if [ "$IS_WSL" = true ]; then
        print_status "Applying WSL2 optimizations..." "step"
        
        # Add WSL configuration to bashrc
        if ! grep -q "# MLD WSL Configuration" ~/.bashrc; then
            cat >> ~/.bashrc << 'EOF'

# MLD WSL Configuration
if grep -q Microsoft /proc/version; then
    # Force headless backends
    export MPLBACKEND=Agg
    export QT_QPA_PLATFORM=offscreen
    export DISPLAY=""
    
    # Graphics optimizations
    export LIBGL_ALWAYS_INDIRECT=1
    export MESA_GL_VERSION_OVERRIDE=3.3
    
    # MPI optimizations
    export OMPI_MCA_btl_vader_single_copy_mechanism=none
fi
EOF
        fi
    fi
    
    # Set optimal thread count
    cpu_cores=$(nproc)
    optimal_threads=$((cpu_cores > 8 ? 8 : cpu_cores))
    echo "export OMP_NUM_THREADS=$optimal_threads" >> ~/.bashrc
    
    print_status "Environment configuration completed" "pass"
}

run_validation() {
    print_header "VALIDATING INSTALLATION"
    
    # Activate environment
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV_NAME"
    
    # Run validation script
    if [ -f "install/validate_installation.py" ]; then
        print_status "Running comprehensive validation..." "step"
        python install/validate_installation.py | tee -a "$INSTALL_LOG"
    else
        print_status "Validation script not found - running basic tests" "warn"
        
        # Basic import tests
        print_status "Testing basic imports..." "step"
        python -c "import numpy, scipy, matplotlib, ase, gpaw, tqdm; print('✅ All core packages imported successfully')" | tee -a "$INSTALL_LOG"
        
        # Test GPAW
        print_status "Testing GPAW..." "step"
        python -c "from gpaw import GPAW, PW; print('✅ GPAW working')" | tee -a "$INSTALL_LOG"
        
        # Test MPI
        print_status "Testing MPI..." "step"
        python -c "from mpi4py import MPI; print(f'✅ MPI working: {MPI.COMM_WORLD.Get_size()} processes')" | tee -a "$INSTALL_LOG"
    fi
    
    print_status "Validation completed" "pass"
}

print_completion_summary() {
    print_header "INSTALLATION COMPLETE"
    
    print_status "MLD Smart Chemistry Framework successfully installed!" "pass"
    print_status "Installation log: $INSTALL_LOG" "info"
    
    echo "" | tee -a "$INSTALL_LOG"
    echo "🎉 Next steps:" | tee -a "$INSTALL_LOG"
    echo "1. Restart your terminal or run: source ~/.bashrc" | tee -a "$INSTALL_LOG"
    echo "2. Activate environment: conda activate $CONDA_ENV_NAME" | tee -a "$INSTALL_LOG"
    echo "3. Test installation: python install/validate_installation.py" | tee -a "$INSTALL_LOG"
    echo "4. Run first example: python scripts/step1_create_h2o.py" | tee -a "$INSTALL_LOG"
    echo "" | tee -a "$INSTALL_LOG"
    
    echo "📁 Environment details:" | tee -a "$INSTALL_LOG"
    echo "   Conda environment: $CONDA_ENV_NAME" | tee -a "$INSTALL_LOG"
    echo "   Python version: $PYTHON_VERSION" | tee -a "$INSTALL_LOG"
    echo "   GPAW datasets: $GPAW_SETUP_PATH" | tee -a "$INSTALL_LOG"
    echo "" | tee -a "$INSTALL_LOG"
}

main() {
    print_header "MLD SMART CHEMISTRY FRAMEWORK INSTALLER"
    
    echo "This script will install the complete MLD computational chemistry environment"
    echo "including GPAW, ASE, and all dependencies for WSL2/Ubuntu."
    echo ""
    echo "Installation will take approximately 15-30 minutes depending on your internet speed."
    echo ""
    echo "Continue? (Y/n)"
    read -r response
    if [[ "$response" =~ ^[Nn]$ ]]; then
        echo "Installation cancelled"
        exit 0
    fi
    
    # Run installation steps
    check_prerequisites
    install_system_dependencies
    install_conda
    create_mld_environment
    download_gpaw_datasets
    configure_environment
    run_validation
    print_completion_summary
    
    print_status "Total installation time: $SECONDS seconds" "info"
}

# Trap errors and provide helpful information
trap 'echo "❌ Installation failed at line $LINENO. Check $INSTALL_LOG for details."' ERR

# Run main installation
main "$@"