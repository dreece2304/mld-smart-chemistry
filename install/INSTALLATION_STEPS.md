# MLD Installation Guide - Manual Steps Required

This guide separates what **you must do manually** vs what **scripts can automate** for setting up the MLD framework on a new machine.

## Phase 1: Manual Prerequisites (USER ACTIONS REQUIRED)

### 1.1 Windows WSL2 Setup (USER)
```powershell
# Run in Windows PowerShell as Administrator
wsl --install
# Restart computer when prompted
```

### 1.2 Ubuntu System Dependencies (USER)
```bash
# Run these commands manually in Ubuntu WSL2
sudo apt update
sudo apt install -y build-essential gfortran cmake pkg-config curl wget git python3-dev

# MPI and scientific libraries (requires sudo)
sudo apt install -y libopenmpi-dev openmpi-bin libscalapack-openmpi-dev
sudo apt install -y libfftw3-dev liblapack-dev libblas-dev libxc-dev

# Graphics libraries (Ubuntu version-specific)
sudo apt install -y libgl1-mesa-dev libglu1-mesa-dev
```

### 1.3 Clone Repository (USER)
```bash
# Clone the repository
git clone https://github.com/yourusername/mld-smart-chemistry.git
cd mld-smart-chemistry
```

## Phase 2: Automated Installation (SCRIPTS HANDLE THIS)

### 2.1 Run Automated Installer
```bash
# This script handles everything that doesn't need sudo
./install/install_conda_and_packages.sh
```

### 2.2 Validate Installation
```bash
# Test everything works
python install/validate_installation.py
```

## Phase 3: Manual Verification (USER ACTIONS)

### 3.1 Test Basic Functionality
```bash
# Activate environment
conda activate mld_modeling

# Test basic optimization
python scripts/step1_create_h2o.py
python scripts/step2_smart_tma.py
```

---

## Troubleshooting Common Issues

### Issue: Package Installation Fails
**Solution**: Run the system dependency commands manually (Phase 1.2)

### Issue: Permission Denied
**Solution**: Ensure you ran `sudo apt install` commands in Phase 1.2

### Issue: GPAW Datasets Missing
**Solution**: Run `gpaw install-data` manually after conda activation

---

## What Each Script Does

### Manual Commands (require sudo/user action):
- WSL2 installation
- System package installation (`sudo apt install`)
- Repository cloning
- Testing and verification

### Automated Scripts (no user intervention):
- Conda installation
- Python package installation
- Environment configuration
- GPAW dataset download
- Validation testing