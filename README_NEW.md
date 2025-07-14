# MLD Smart Chemistry Framework

A comprehensive computational chemistry framework for **Molecular Layer Deposition (MLD)** simulations using DFT (Density Functional Theory) with intelligent optimization strategies.

## 🚀 Quick Start

```bash
# Clone repository
git clone https://github.com/yourusername/mld-smart-chemistry.git
cd mld-smart-chemistry

# Automated installation (WSL2/Ubuntu)
./install/setup_complete_environment.sh

# Validate installation
python install/validate_installation.py

# Run first example
python scripts/step1_create_h2o.py
```

## 🧪 What This Framework Does

### Smart Molecular Optimization
- **Database → Classical → DFT** hierarchy for efficient optimization
- Automatic lookup in multiple molecular databases (PubChem, NIST, ChEBI)
- Intelligent fallback to classical force fields and full DFT when needed
- **10-100x faster** than traditional DFT-only approaches for known molecules

### MLD Chemistry Simulation
- **TMA + H2O → Al-OH + CH4** gas phase reactions
- **Surface chemistry** on hydroxylated Si(100) surfaces
- **Complete MLD cycles** with realistic temperature/pressure conditions
- **Energy analysis** and reaction pathway visualization

### Progress Tracking & Visualization
- Real-time optimization progress with force convergence monitoring
- Headless visualization compatible with WSL2
- Interactive HTML reports for molecular structures
- Trajectory analysis and convergence plots

## 📁 Project Structure

```
mld-smart-chemistry/
├── src/mld_chemistry/          # Core library modules
├── scripts/                    # Step-by-step workflow scripts
├── install/                    # Automated installation tools
├── structures/                 # Input molecular structures
├── tests/                      # Validation and testing
├── examples/                   # Usage examples
└── docs/                       # Documentation
```

## 🔬 Supported Chemistry

### Molecules
- **TMA (Trimethylaluminum)**: Al(CH₃)₃ - MLD precursor
- **H2O (Water)**: Surface hydroxylation and reaction partner
- **Silicon surfaces**: Si(100) with realistic OH termination

### Reactions
- **Gas phase**: TMA + H2O → Al(CH₃)₂OH + CH4
- **Surface**: TMA + Surface-OH → Surface-O-Al(CH₃)₂ + CH4
- **Complete MLD cycle**: Precursor + co-reactant + surface

## 🛠️ Installation

### Requirements
- **Windows 10/11** with WSL2 enabled
- **Ubuntu 20.04+** in WSL2
- **8GB+ RAM** (16GB recommended for large systems)
- **Internet connection** for initial setup

### Automated Installation
```bash
# Run the master installer
./install/setup_complete_environment.sh
```

This will:
1. Install system dependencies (build tools, MPI, scientific libraries)
2. Create isolated conda environment with Python 3.10
3. Install GPAW, ASE, and all required packages
4. Download GPAW PAW datasets
5. Validate complete installation
6. Run performance benchmarks

## 📖 Usage Examples

### Basic Molecule Optimization
```python
from mld_chemistry import optimize_molecule

# Smart optimization with database lookup
atoms, log = optimize_molecule("trimethylaluminum")
print(f"Optimization path: {log[-1]}")
```

### Step-by-Step Workflow
```bash
# Step 1: Optimize individual molecules
python scripts/step2_smart_tma.py      # Optimize TMA
python scripts/step3_smart_h2o.py      # Optimize H2O

# Step 2: Create surface
python scripts/step4_create_surface.py # Si(100) + OH

# Step 3: Simulate reactions  
python scripts/step5_tma_h2o_reaction.py # Gas phase reaction

# Step 4: Visualize results
python scripts/visualize_results.py
```

## 🧪 Testing

```bash
# Run complete test suite
python -m pytest tests/

# Quick validation
python install/validate_installation.py

# Performance benchmark
python tests/benchmark_performance.py
```

## 📊 Multi-Machine Deployment

### Setup on New Machine
1. **Install WSL2** on target Windows machine
2. **Clone repository**: `git clone <repo-url>`
3. **Run installer**: `./install/setup_complete_environment.sh`
4. **Validate**: `python install/validate_installation.py`
5. **Benchmark**: Compare performance with reference machine

## 🔬 Scientific Background

This framework implements modern computational chemistry best practices:

- **Density Functional Theory (DFT)** using GPAW with PAW method
- **PBE exchange-correlation functional** for accurate energetics  
- **Plane wave basis sets** with appropriate cutoff energies
- **Smart optimization strategies** combining multiple approaches
- **Realistic MLD conditions** including temperature and pressure effects

## 📄 License

MIT License - see [LICENSE](LICENSE) for details.

---

*Built with ❤️ for the computational chemistry community*