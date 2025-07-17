# Usage Guide

## Quick Start Commands

### Immediate Testing
```bash
# Test environment
python tests/test_suite.py --level 1

# Quick functionality test
python examples/quick_run.py --test

# Run first simulation
python examples/quick_run.py --molecule tma --precision fast
```

### Main Workflows
```bash
# Fast MLD simulation
python src/simulation/dft_fast_mld.py

# Production MLD simulation
python src/simulation/dft_small_mld.py

# Surface optimization
python src/simulation/optimize_surface_only.py surface.xyz
```

---

## Progress Tracking Features

### Real-time Monitoring
All simulation scripts include comprehensive progress tracking:

- **🎯 Convergence Progress**: Shows percentage towards force convergence
- **📈 Step Progress**: Current optimization step with ETA
- **💻 System Resources**: CPU usage and memory consumption
- **🧪 Workflow Progress**: Overall completion for multi-stage calculations

### Visual Progress Examples
```
🎯 TMA_Molecule Convergence: ████████░░ 85% [00:45<00:08]
📈 TMA_Molecule Steps: 15/100 | E:-6.234 ΔE:-2.1e-04 F_max:0.0234 ETA:0:02:15
💻 System Resources: CPU 45.2% | Mem: 1250MB (120%)
```

### Monitoring Commands
```bash
# Real-time DFT progress monitoring
python scripts/monitoring/monitor_dft_progress.py output.log

# With custom refresh rate
python scripts/monitoring/monitor_dft_progress.py output.log --refresh 0.5

# Terminal-based result visualization
python src/visualization/terminal_visualizer.py
```

---

## Examples & Use Cases

### 1. Basic Molecule Optimization
```bash
# TMA molecule optimization
python examples/quick_run.py --molecule tma --precision fast

# 2-butyne-1,4-diol optimization
python examples/quick_run.py --molecule diol --precision fast
```

### 2. Surface Modeling
```bash
# Surface optimization (can be run independently)
python src/simulation/optimize_surface_only.py small_si100_2x2_4L.xyz

# With restart capability
python src/simulation/optimize_surface_only.py small_si100_2x2_4L.xyz --restart

# Monitor surface optimization
python scripts/monitoring/monitor_dft_progress.py small_si100_2x2_4L_surf_opt.log
```

### 3. Complete MLD Cycles
```bash
# Fast simulation with pre-optimized structures
python src/simulation/dft_fast_mld.py

# Full production simulation
python src/simulation/dft_small_mld.py

# Complete workflow with enhanced progress tracking
python src/simulation/run_mld_with_progress.py
```

### 4. Analysis & Visualization
```bash
# Terminal-based visualization
python src/visualization/terminal_visualizer.py

# Analysis of bulk structures
python src/analysis/mld_analyzer.py

# Enhanced analysis
python src/analysis/enhanced_analysis.py
```

---

## Configuration & Parameters

### Precision Modes
```bash
# Fast mode (quick testing)
python examples/quick_run.py --full-cycle --precision fast

# Medium mode (balanced)
python examples/quick_run.py --full-cycle --precision medium

# Production mode (high accuracy)
python examples/quick_run.py --full-cycle --precision production
```

### DFT Parameters
The framework automatically adjusts DFT parameters based on system size:

- **Small systems (<50 atoms)**: High precision, tight convergence
- **Medium systems (50-100 atoms)**: Balanced parameters
- **Large systems (100-200 atoms)**: Optimized for efficiency

### Energy Caching
The system automatically caches calculated energies in `config/pre_calculated_energies.json`:
```json
{
  "surface": -654.3515,
  "tma_molecule": -123.4567,
  "h2o_molecule": -45.6789
}
```

---

## Background Execution

### For Long Calculations
```bash
# Background execution with output logging
nohup python src/simulation/dft_small_mld.py > mld_output.log 2>&1 &

# Save process ID for monitoring
echo $! > simulation.pid

# Monitor progress later
python scripts/monitoring/monitor_dft_progress.py mld_output.log
```

### Process Management
```bash
# Check if simulation is running
ps aux | grep python

# Check specific process
ps -p $(cat simulation.pid)

# Stop simulation gracefully
kill -TERM $(cat simulation.pid)
```

---

## Error Handling & Troubleshooting

### Common Issues

#### 1. GPAW Setup Path
```bash
# Error: GPAW setup files not found
export GPAW_SETUP_PATH=~/.local/share/gpaw/gpaw-setups-24.11.0

# Verify setup
python -c "import gpaw; print(gpaw.__file__)"
```

#### 2. Memory Issues
```bash
# Monitor memory usage
free -h

# Use fast simulation for large systems
python src/simulation/dft_fast_mld.py

# Reduce system size
python scripts/create_small_surfaces.py
```

#### 3. Convergence Problems
```bash
# Check convergence progress
python scripts/monitoring/monitor_dft_progress.py output.log

# Use restart functionality
python src/simulation/optimize_surface_only.py surface.xyz --restart
```

### Diagnostic Commands
```bash
# Comprehensive system check
python tests/test_suite.py

# Specific level testing
python tests/test_suite.py --level 1-3

# Installation validation
python scripts/setup/validate_installation.py
```

---

## Output Files & Results

### File Structure
```
data/outputs/
├── structures/          # Optimized structures (.xyz files)
├── plots/              # Generated plots (.png files)
└── dft_logs/           # DFT output logs (.txt files)
```

### Key Output Files
- `*_optimized.xyz`: Final optimized structures
- `*_optimization_progress.json`: Detailed optimization data
- `*_optimization_progress.png`: Convergence plots
- `*_surf_results.json`: Surface optimization results

### Analysis Files
- `config/pre_calculated_energies.json`: Energy cache
- `*_analysis_report.html`: Comprehensive analysis reports
- `*_trajectory.xyz`: Optimization trajectory

---

## Advanced Features

### Restart Capability
```bash
# Surface optimization with restart
python src/simulation/optimize_surface_only.py surface.xyz --restart

# The system will automatically:
# 1. Look for completed optimization
# 2. Check for partial structures
# 3. Load saved GPAW state
# 4. Continue from where it left off
```

### Intelligent Initial Guess
The framework uses a smart hierarchy:
1. **Database lookup**: Check for known structures
2. **Pre-calculated energies**: Use cached values
3. **Classical methods**: Fast initial optimization
4. **Full DFT**: When needed

### Automatic Parameter Optimization
- **System size detection**: Adjusts parameters based on atom count
- **Convergence criteria**: Adaptive thresholds
- **Memory management**: Prevents system overload
- **Error recovery**: Automatic restart with modified parameters

---

## Performance Tips

### For Faster Calculations
```bash
# Use fast simulation mode
python src/simulation/dft_fast_mld.py

# Use pre-optimized structures
# (automatically loaded from config/pre_calculated_energies.json)

# Monitor system resources
htop
```

### For Large Systems
```bash
# Create smaller surface models
python scripts/setup/create_small_surfaces.py

# Use progressive optimization
python src/simulation/optimize_surface_only.py surface.xyz

# Monitor memory usage
python scripts/monitoring/monitor_dft_progress.py --memory output.log
```

### For Remote Execution
```bash
# Use SSH with compression
ssh -C -p 2222 dreece23@172.25.145.180

# Background execution
nohup python src/simulation/dft_small_mld.py > output.log 2>&1 &

# Monitor remotely
python scripts/monitoring/monitor_dft_progress.py output.log
```

---

## Validation & Testing

### Progressive Test Suite
```bash
# Level 1: Basic imports and packages
python tests/test_suite.py --level 1

# Level 2: GPAW functionality
python tests/test_suite.py --level 2

# Level 3: Optimization algorithms
python tests/test_suite.py --level 3

# All levels
python tests/test_suite.py
```

### Individual Component Testing
```bash
# Test visualization
python src/visualization/terminal_visualizer.py

# Test analysis
python src/analysis/mld_analyzer.py

# Test monitoring
python scripts/monitoring/monitor_dft_progress.py test_output.log
```

---

## Integration with Experimental Work

### Output for Experimental Comparison
```python
# Example: Extract growth rate data
from src.analysis.mld_analyzer import MLDAnalyzer

analyzer = MLDAnalyzer('data/outputs/structures/')
growth_rate = analyzer.calculate_growth_rate()
thickness = analyzer.measure_thickness()
```

### Property Calculations
```python
# Example: Calculate bulk properties
from src.analysis.enhanced_analysis import EnhancedAnalyzer

analyzer = EnhancedAnalyzer('final_structure.xyz')
density = analyzer.calculate_density()
mechanical_properties = analyzer.estimate_mechanical_properties()
```

---

## Best Practices

### 1. Always Test First
```bash
# Start with quick test
python examples/quick_run.py --test

# Then try small system
python examples/quick_run.py --molecule tma --precision fast
```

### 2. Monitor Progress
```bash
# Use real-time monitoring
python scripts/monitoring/monitor_dft_progress.py output.log

# Check system resources
htop
free -h
```

### 3. Use Restart Capability
```bash
# For long calculations
python src/simulation/optimize_surface_only.py surface.xyz --restart

# Save progress regularly (automatic every 5 steps)
```

### 4. Keep Results Organized
```bash
# Results are automatically organized in data/outputs/
ls data/outputs/structures/  # Optimized structures
ls data/outputs/plots/       # Generated plots
ls data/outputs/dft_logs/    # DFT output logs
```

---

## Quick Reference

### Essential Commands
```bash
# Test environment
python tests/test_suite.py --level 1

# Quick simulation
python examples/quick_run.py --test

# Main simulation
python src/simulation/dft_fast_mld.py

# Monitor progress
python scripts/monitoring/monitor_dft_progress.py output.log

# Visualize results
python src/visualization/terminal_visualizer.py
```

### File Locations
- **Simulations**: `src/simulation/`
- **Analysis**: `src/analysis/`
- **Visualization**: `src/visualization/`
- **Monitoring**: `scripts/monitoring/`
- **Results**: `data/outputs/`
- **Configuration**: `config/`

---

**Last Updated**: July 2025