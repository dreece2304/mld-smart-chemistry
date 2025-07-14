# MLD Simulation with Progress Tracking - User Guide

## Overview

This enhanced MLD simulation system provides comprehensive progress tracking for all calculations using beautiful progress bars, real-time monitoring, and detailed reporting. You can now run complex DFT calculations with full visibility into their progress.

## 🚀 Quick Start

### 1. Test Everything Works
```bash
python quick_run.py --test
```
This runs a quick H2 molecule optimization to verify GPAW is working.

### 2. Run Your First MLD Calculation
```bash
python quick_run.py --molecule tma --precision fast
```
Optimizes the TMA molecule with progress tracking.

### 3. Run Complete MLD Cycle
```bash
python quick_run.py --full-cycle --precision fast
```
Runs the complete TMA + diol deposition simulation.

## 📊 Progress Tracking Features

### Multi-Level Progress Bars
- **🎯 Convergence Progress**: Shows percentage towards force convergence (logarithmic scale)
- **📈 Step Progress**: Current optimization step with ETA calculation
- **💻 System Resources**: CPU usage and memory consumption
- **🧪 Workflow Progress**: Overall completion for multi-stage calculations

### Real-Time Information
- Energy per step with changes
- Force convergence tracking
- Time per step and ETA estimates
- Memory usage monitoring
- Automatic final summaries

### Visual Examples
```
🎯 TMA_Molecule Convergence: ████████░░ 85% [00:45<00:08]
📈 TMA_Molecule Steps: 15/100 steps | E:-6.234 ΔE:-2.1e-04 F_max:0.0234 ETA:0:02:15
💻 System Resources: CPU 45.2% | Mem: 1250MB (120%)
```

## 🛠️ Available Commands

### Quick Interface (`quick_run.py`)

#### Basic Tests
```bash
python quick_run.py --test                    # Quick functionality test
python quick_run.py --status                  # Check installation status
```

#### Single Components
```bash
python quick_run.py --molecule tma            # Optimize TMA molecule
python quick_run.py --molecule diol           # Optimize diol molecule  
python quick_run.py --surface                 # Optimize surface only
```

#### Complete Workflows
```bash
python quick_run.py --full-cycle              # Complete MLD cycle
python quick_run.py --full-cycle --precision medium  # Higher accuracy
```

#### Analysis
```bash
python quick_run.py --radiation               # Radiation damage test
python quick_run.py --analysis                # Bulk structure analysis
```

### Advanced Interface (`run_mld_with_progress.py`)

#### Precision Levels
```bash
# Fast (testing): PW(200), loose convergence, ~minutes
python run_mld_with_progress.py --precision fast

# Medium (research): PW(300), medium convergence, ~hours  
python run_mld_with_progress.py --precision medium

# Production (publication): PW(400), tight convergence, ~days
python run_mld_with_progress.py --precision production
```

#### Specific Tasks
```bash
python run_mld_with_progress.py --molecule tma --precision fast
python run_mld_with_progress.py --output my_results --precision medium
```

### Direct Optimization (`progress_optimization.py`)

For custom calculations:
```python
from progress_optimization import optimize_structure

# Your GPAW parameters
calc_params = {
    'mode': PW(300),
    'xc': 'PBE',
    'kpts': (2, 2, 1)
}

# Run with progress tracking
atoms, stats = optimize_structure(
    your_atoms, calc_params,
    optimizer='BFGS', fmax=0.05,
    name="My_Calculation"
)
```

## 📈 Real-Time Monitoring

### Live Monitoring (`monitor_calculations.py`)
Monitor running calculations in real-time:

```bash
# Continuous monitoring with 5-second updates
python monitor_calculations.py

# Monitor specific directory
python monitor_calculations.py --dir /path/to/calculations

# Single status check
python monitor_calculations.py --once

# Custom update intervals
python monitor_calculations.py --interval 10 --plot-interval 120
```

### Monitoring Features
- **📊 Active Calculation Detection**: Automatically finds running GPAW jobs
- **💾 File Growth Tracking**: Shows file sizes and growth rates
- **🔧 GPAW Progress**: Parses SCF iterations and energies
- **📈 Optimization Status**: Tracks force convergence from log files
- **📊 Automatic Plots**: Generates convergence plots periodically
- **⏰ Timestamps**: Shows when calculations were last updated

### Monitoring Display Example
```
🔍 MLD CALCULATION MONITOR
⏰ 2025-01-14 15:30:45 | Runtime: 0:15:23
================================================================================

📊 Calculation 1: results/tma_opt_gpaw.out
----------------------------------------
💾 File size: 15.23 MB | Growth: 0.156 MB/s
🔧 GPAW calculation:
   SCF steps: 45
   Last energy: -234.567890 eV
   Last max force: 0.034567 eV/Å
   Status: 🔄 Running
⏱️  Last updated: 3 seconds ago

Press Ctrl+C to stop monitoring | Next update in 5s
```

## 🔬 Enhanced Analysis

### Automated Analysis (`enhanced_analysis.py`)

#### Complete Analysis Suite
```bash
python enhanced_analysis.py structure.xyz
python enhanced_analysis.py structure.xyz --name "My_Structure"
```

#### Specific Analysis Types
```bash
python enhanced_analysis.py structure.xyz --bulk-only
python enhanced_analysis.py structure.xyz --radiation-only
python enhanced_analysis.py structure.xyz --no-radiation
```

### Analysis Features
- **🔍 Bulk Structure Analysis**: Crystallinity, density, mechanical properties
- **☢️ Radiation Damage**: UV photolysis and electron beam damage
- **📊 Progress Tracking**: Visual progress bars for all analysis steps
- **📄 Comprehensive Reports**: Detailed results with statistics
- **🔗 Integrated Workflow**: Seamless integration with optimization results

### Analysis Progress Display
```
🔬 Enhanced Analysis: TMA_Optimized
📁 Structure: tma_optimized.xyz
🕒 Started: 15:45:30
============================================================

📊 Bulk Structure Analysis
----------------------------------------
🔍 Bulk Analysis: ████████░░ 80% | 6/8 steps [00:12<00:03]

☢️  Radiation Damage Analysis  
----------------------------------------
☢️  Radiation Tests: ██████░░░░ 60% | 3/5 tests [00:08<00:05]
```

## 🎯 Precision Levels Guide

### Fast (Development/Testing)
- **Time**: Minutes to hours
- **Settings**: PW(200), loose convergence
- **Use for**: Method development, testing, debugging
- **Accuracy**: ~0.01 eV

### Medium (Research)
- **Time**: Hours to days  
- **Settings**: PW(300), medium convergence
- **Use for**: Research calculations, trends
- **Accuracy**: ~0.005 eV

### Production (Publication)
- **Time**: Days to weeks
- **Settings**: PW(400), tight convergence  
- **Use for**: Final results, publication quality
- **Accuracy**: ~0.0005 eV

## 📊 Understanding Progress Information

### Convergence Progress Bar
Shows logarithmic progress towards force convergence:
- **0%**: Just started
- **50%**: Halfway to target on log scale
- **100%**: Converged to target force

### Step Information
- **E**: Current energy in eV
- **ΔE**: Energy change from previous step
- **F_max**: Maximum force on any atom
- **ETA**: Estimated time to completion

### System Resources
- **CPU**: Current CPU usage percentage
- **Mem**: Current memory usage and change from start

## 🗂️ Output Files

### Structure Files
- `*_optimized.xyz`: Final optimized structures
- `*_opt.traj`: Full optimization trajectory
- `*_opt.log`: Optimization log with forces/energies

### Progress and Logs
- `*_gpaw.out`: GPAW DFT calculation output
- `*_opt.log`: ASE optimization log
- Progress bars appear in terminal (not saved)

### Analysis Reports
- `*_bulk_analysis.txt`: Detailed bulk properties
- `*_radiation_analysis.txt`: Radiation damage results  
- `*_complete_analysis.txt`: Combined analysis summary
- `mld_workflow_report.txt`: Complete workflow summary

### Monitoring
- `convergence_plot.png`: Real-time convergence plots
- Monitoring output appears in terminal

## 🚨 Troubleshooting

### Common Issues

#### Progress bars not appearing
```bash
pip install tqdm  # Install progress bar library
```

#### GPAW datasets missing
```bash
gpaw install-data ~/.local/share/gpaw
export GPAW_SETUP_PATH=~/.local/share/gpaw/gpaw-setups-24.11.0
```

#### Calculations running but no progress
- Check if calculation is actually running: `python monitor_calculations.py --once`
- Look for error messages in `*_gpaw.out` files
- Verify GPAW setup with: `python quick_run.py --test`

#### Memory issues
- Use lower precision: `--precision fast`
- Reduce system size
- Monitor with: `python monitor_calculations.py`

### Getting Help

1. **Check Status**: `python quick_run.py --status`
2. **Run Test**: `python quick_run.py --test`  
3. **Monitor Progress**: `python monitor_calculations.py --once`
4. **Check Logs**: Look at `*_gpaw.out` and `*_opt.log` files

## 🎮 Interactive Usage Tips

### Running in Background
```bash
# Run with output to file
python quick_run.py --full-cycle > calculation.log 2>&1 &

# Monitor from another terminal
python monitor_calculations.py
```

### Parallel Calculations
The system automatically detects MPI and adjusts progress bars:
```bash
mpirun -np 4 python run_mld_with_progress.py --precision medium
```

### Jupyter Notebook
All scripts work in Jupyter with enhanced progress bars:
```python
%run quick_run.py --test
```

### Customization
All progress tracking can be customized by modifying the respective Python files:
- `progress_optimization.py`: Core optimization progress
- `run_mld_with_progress.py`: Workflow progress  
- `monitor_calculations.py`: Monitoring settings
- `enhanced_analysis.py`: Analysis progress

## 📈 Performance Tips

1. **Start with Fast**: Always test with `--precision fast` first
2. **Monitor Resources**: Use the monitoring script to track CPU/memory
3. **Use Appropriate Precision**: Don't use production settings for testing
4. **Check Convergence**: Watch the convergence progress bar behavior
5. **Save Intermediate Results**: Scripts automatically save progress

## 🎉 Example Full Workflow

```bash
# 1. Test installation
python quick_run.py --test

# 2. Check status  
python quick_run.py --status

# 3. Run single molecule test
python quick_run.py --molecule tma --precision fast

# 4. Monitor the calculation (in another terminal)
python monitor_calculations.py

# 5. Run complete MLD cycle
python quick_run.py --full-cycle --precision medium

# 6. Analyze results
python quick_run.py --analysis

# 7. Test radiation effects
python quick_run.py --radiation
```

This gives you a complete MLD simulation with full progress visibility and comprehensive analysis!