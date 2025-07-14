#!/bin/bash
# WSL2 Ubuntu 24.04 Environment Setup for MLD
# Handles updated package names and WSL-specific requirements

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
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
    echo -e "${CYAN}"
    echo "==========================================="
    echo "$1"
    echo "==========================================="
    echo -e "${NC}"
}

detect_ubuntu_version() {
    if [ -f /etc/os-release ]; then
        . /etc/os-release
        echo "$VERSION_ID"
    else
        echo "unknown"
    fi
}

install_system_dependencies() {
    print_header "INSTALLING SYSTEM DEPENDENCIES"
    
    # Update package list
    print_status "Updating package list..." "step"
    sudo apt update
    
    # Detect Ubuntu version for correct package names
    ubuntu_version=$(detect_ubuntu_version)
    print_status "Detected Ubuntu version: $ubuntu_version" "info"
    
    # Essential build tools
    print_status "Installing build essentials..." "step"
    sudo apt install -y \
        build-essential \
        gfortran \
        cmake \
        pkg-config \
        curl \
        wget \
        git
    
    # MPI and parallel computing
    print_status "Installing MPI and parallel computing tools..." "step"
    sudo apt install -y \
        libopenmpi-dev \
        openmpi-bin \
        libscalapack-openmpi-dev
    
    # Scientific computing libraries
    print_status "Installing scientific libraries..." "step"
    sudo apt install -y \
        libfftw3-dev \
        liblapack-dev \
        libblas-dev \
        libhdf5-serial-dev
    
    # Graphics and display libraries (Ubuntu 24.04 compatible)
    print_status "Installing graphics libraries..." "step"
    
    if [[ "$ubuntu_version" == "24.04" || "$ubuntu_version" > "22.04" ]]; then
        # Ubuntu 24.04+ package names
        sudo apt install -y \
            libgl1-mesa-dev \
            libglu1-mesa-dev \
            libglx-mesa0 \
            libxrender1 \
            libxkbcommon-x11-0 \
            libxcb-xinerama0 \
            libxrandr2 \
            libxss1 \
            libgconf-2-4 \
            libgtk-3-0 \
            libdrm2 \
            libxcomposite1 \
            libxdamage1 \
            libxfixes3 || {
                print_status "Some graphics packages failed - continuing" "warn"
            }
    else
        # Older Ubuntu versions
        sudo apt install -y \
            libgl1-mesa-glx \
            libglu1-mesa-dev \
            libxrender1 \
            libxkbcommon-x11-0 \
            libxcb-xinerama0 \
            libxrandr2 \
            libxss1 || {
                print_status "Some graphics packages failed - continuing" "warn"
            }
    fi
    
    print_status "System dependencies installed" "pass"
}

configure_wsl_environment() {
    print_header "CONFIGURING WSL2 ENVIRONMENT"
    
    # Check if we're actually in WSL
    if ! grep -q Microsoft /proc/version; then
        print_status "Not running in WSL - skipping WSL configuration" "info"
        return 0
    fi
    
    print_status "Configuring WSL2-specific settings..." "step"
    
    # Create WSL configuration block in bashrc
    if ! grep -q "# WSL MLD Configuration" ~/.bashrc; then
        cat >> ~/.bashrc << 'EOF'

# WSL MLD Configuration
if grep -q Microsoft /proc/version; then
    # Force headless backends for WSL
    export MPLBACKEND=Agg
    export QT_QPA_PLATFORM=offscreen
    export DISPLAY=""
    
    # Disable GUI warnings
    export LIBGL_ALWAYS_INDIRECT=1
    export MESA_GL_VERSION_OVERRIDE=3.3
    
    # MPI settings for WSL
    export OMPI_MCA_btl_vader_single_copy_mechanism=none
fi
EOF
        print_status "Added WSL configuration to ~/.bashrc" "pass"
    else
        print_status "WSL configuration already exists in ~/.bashrc" "info"
    fi
    
    # Set for current session
    export MPLBACKEND=Agg
    export QT_QPA_PLATFORM=offscreen
    export DISPLAY=""
    export LIBGL_ALWAYS_INDIRECT=1
    export MESA_GL_VERSION_OVERRIDE=3.3
    export OMPI_MCA_btl_vader_single_copy_mechanism=none
    
    print_status "WSL2 environment configured" "pass"
}

fix_conda_permissions() {
    print_header "FIXING CONDA PERMISSIONS"
    
    # Sometimes conda has permission issues in WSL
    if [ -d "$HOME/miniconda3" ]; then
        print_status "Fixing conda permissions..." "step"
        chmod -R u+rwx "$HOME/miniconda3" 2>/dev/null || true
        print_status "Conda permissions fixed" "pass"
    fi
    
    # Ensure conda is in PATH
    if ! command -v conda &> /dev/null; then
        if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
            print_status "Sourcing conda environment..." "step"
            source "$HOME/miniconda3/etc/profile.d/conda.sh"
            print_status "Conda environment sourced" "pass"
        else
            print_status "Conda not found - please install miniconda first" "fail"
            exit 1
        fi
    else
        print_status "Conda is available" "pass"
    fi
}

setup_tmp_directory() {
    print_header "SETTING UP TEMPORARY DIRECTORIES"
    
    # Ensure temp directories exist and are writable
    for tmpdir in "/tmp" "$HOME/tmp" "$HOME/.cache"; do
        if [ ! -d "$tmpdir" ]; then
            mkdir -p "$tmpdir"
            print_status "Created directory: $tmpdir" "info"
        fi
        
        if [ ! -w "$tmpdir" ]; then
            chmod 755 "$tmpdir" 2>/dev/null || {
                print_status "Cannot write to $tmpdir - may cause issues" "warn"
            }
        fi
    done
    
    print_status "Temporary directories configured" "pass"
}

test_system_setup() {
    print_header "TESTING SYSTEM SETUP"
    
    # Test essential commands
    commands=("gcc" "gfortran" "mpirun" "cmake")
    
    for cmd in "${commands[@]}"; do
        if command -v "$cmd" &> /dev/null; then
            version=$($cmd --version 2>&1 | head -n1 || echo "unknown")
            print_status "$cmd: available ($version)" "pass"
        else
            print_status "$cmd: not found" "fail"
        fi
    done
    
    # Test MPI functionality
    print_status "Testing MPI..." "step"
    if mpirun -np 2 echo "MPI test" &> /dev/null; then
        print_status "MPI: working" "pass"
    else
        print_status "MPI: may have issues" "warn"
    fi
    
    # Test graphics libraries
    print_status "Testing graphics libraries..." "step"
    if ldconfig -p | grep -q "libGL.so"; then
        print_status "OpenGL libraries: found" "pass"
    else
        print_status "OpenGL libraries: not found (may affect visualization)" "warn"
    fi
    
    print_status "System setup test completed" "pass"
}

main() {
    print_header "WSL2 UBUNTU 24.04 MLD SETUP"
    
    echo "This script will:"
    echo "1. Install system dependencies (Ubuntu 24.04 compatible)"
    echo "2. Configure WSL2 environment variables"
    echo "3. Fix conda permissions"
    echo "4. Set up temporary directories"
    echo "5. Test system setup"
    echo ""
    echo "Continue? (Y/n)"
    read -r response
    if [[ "$response" =~ ^[Nn]$ ]]; then
        echo "Setup cancelled"
        exit 0
    fi
    
    # Run setup steps
    install_system_dependencies
    configure_wsl_environment
    fix_conda_permissions
    setup_tmp_directory
    test_system_setup
    
    print_header "WSL2 SETUP COMPLETE"
    print_status "System dependencies installed and configured" "pass"
    print_status "WSL2 environment optimized for MLD" "pass"
    
    echo ""
    echo -e "${GREEN}🎉 WSL2 setup complete!${NC}"
    echo ""
    echo "Next steps:"
    echo "1. source ~/.bashrc  # Load new environment variables"
    echo "2. ./clean_install_mld.sh  # Install MLD packages"
    echo "3. python validate_clean_install.py  # Test installation"
    echo ""
    echo "Note: You may need to restart your terminal for all changes to take effect."
}

# Run main function
main "$@"