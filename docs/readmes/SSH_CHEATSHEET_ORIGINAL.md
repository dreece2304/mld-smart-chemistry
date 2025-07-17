# MLD SSH & Simulation Cheat Sheet
*Complete guide for remote MLD simulations on lab desktop*

---

## 🔗 Connecting to Lab Desktop

### Basic SSH Connection
```bash
ssh -p 2222 dreece23@172.25.145.180
```

### Quick Connection (if SSH config set up)
```bash
ssh lab-desktop
```

### Connection Troubleshooting
```bash
# Test connection diagnostics
./scripts/ssh_setup/04_test_connection.sh

# If connection fails, restart SSH on lab desktop (via Remote Desktop):
~/restart_ssh.sh
~/check_ssh_status.sh

# Check Windows port forwarding (lab desktop PowerShell as Admin):
netsh interface portproxy show all
```

### VPN Requirements
- Must be connected to **UW VPN (Big-IP Edge)** for access
- Lab Desktop IP: `172.25.145.180`
- SSH Port: `2222`

---

## 📊 Checking Simulation Progress

### Check Running Processes
```bash
# See all Python/DFT processes
ps aux | grep python
ps aux | grep gpaw

# See running background jobs
jobs

# System resource monitoring
htop        # Interactive process viewer
top         # Simple process viewer
```

### View Simulation Logs
```bash
# Navigate to project directory
cd ~/mld-smart-chemistry

# List all log files
ls -la *.log *.txt *.out

# Monitor real-time output
tail -f output.log                    # General output
tail -f *_progress.txt                # Progress logs
tail -f *_gpaw.out                    # GPAW DFT output

# View last 50 lines of logs
tail -50 dft_output.log
tail -50 mld_simulation.log
```

### Check DFT Calculation Status
```bash
# Look for GPAW output files
ls -la *_gpaw.out *_test.txt *_dft.txt

# Check convergence status
grep -i "converged" *.txt *.out
grep -i "scf.*converged" *.out

# Check for calculation errors
grep -i "error\|failed\|crash\|exception" *.txt *.out

# Monitor SCF iterations
grep -i "scf\|iteration" *.out | tail -20
```

### Progress Tracking
```bash
# Check smart caching status
ls -la geometry_cache/ dft_cache/

# View progress summary
grep -i "progress\|cycle\|step" *.log
```

---

## 🚀 Simulation Commands

### Test Suite (Run First)
```bash
# Test DFT setup with 4 progressive tests
python scripts/dft_simple_test.py

# Expected runtime: 5-15 minutes
# Tests: H2O (3 atoms) → TMA (13 atoms) → Surface cluster (30 atoms) → H2O+Surface (40 atoms)
```

### Surface Generation
```bash
# Create literature-validated small surfaces
python scripts/create_small_surfaces.py

# Expected runtime: 1-2 minutes
# Creates: Si(100) 2×2, 3×3 slabs with proper hydroxylation
```

### Full MLD Simulation
```bash
# Run complete MLD simulation with DFT
python scripts/dft_small_mld.py

# Expected runtime: 2-6 hours
# Features: TMA/H2O cycles, smart caching, progress tracking
```

### Smart Features Testing
```bash
# Test intelligent initial guess system
python scripts/smart_initial_guess.py

# Test progress tracking demos
python scripts/dft_with_progress.py

# Add progress bars to existing scripts
python scripts/add_progress_to_dft.py
```

### Background Execution (Recommended for Long Jobs)
```bash
# Run simulation in background with logging
nohup python scripts/dft_small_mld.py > mld_output.log 2>&1 &

# Check background job status
jobs

# View background job output
tail -f mld_output.log
```

---

## 🔧 Troubleshooting & Diagnostics

### DFT Convergence Issues
```bash
# Fix surface SCF convergence problems
python scripts/fix_surface_convergence.py

# Progressive DFT parameter testing
python scripts/dft_troubleshoot.py

# Validate molecular geometries (prevents "atoms too close" errors)
python scripts/geometry_validator.py
```

### System Health Checks
```bash
# Memory usage
free -h
df -h                    # Disk space

# CPU utilization
top
htop

# Network connectivity
ping 8.8.8.8            # Internet
ping 172.25.145.180     # Lab desktop from laptop
```

### Common Error Solutions
```bash
# "Atoms too close" error:
python scripts/geometry_validator.py

# SCF convergence failure:
python scripts/fix_surface_convergence.py

# Memory issues:
# Reduce system size or use smaller cutoff energies

# SSH connection drops:
# Check VPN connection, restart SSH service
```

---

## 📋 Development Roadmap

### ✅ **Completed Features**
- **Smart Initial Guess System**
  - Geometry similarity matching
  - Material-specific DFT parameters
  - Automatic caching and reuse
  - 5-10x speedup for similar systems

- **Real-time Progress Tracking**
  - SCF iteration monitoring
  - Optimization step progress
  - ASCII convergence plots
  - Time estimation and ETA

- **SSH Remote Development**
  - Laptop ↔ lab desktop connection
  - Key-based authentication
  - File transfer automation
  - Background job management

- **Literature-Validated Models**
  - Small surface generation (<200 atoms)
  - Geometry validation
  - Conservative DFT parameters
  - Troubleshooting tools

### 🔄 **In Progress**
- DFT test suite validation on lab hardware
- Surface convergence optimization
- Performance benchmarking

### 📝 **Next Implementation Priorities**

#### **High Priority**
1. **Advanced GUI Visualization**
   - QT-based molecular viewer with real-time updates
   - Interactive structure editing
   - Progress monitoring dashboard
   - Publication-quality rendering

2. **Enhanced DFT Workflow**
   - Automatic parameter optimization based on system type
   - Parallel calculation management
   - Result validation and quality checks
   - Error recovery and restart capabilities

3. **MLD-Specific Analysis Tools**
   - Growth rate calculation and visualization
   - Surface coverage analysis over cycles
   - Reaction pathway identification
   - Thermodynamic property calculation

#### **Medium Priority**
1. **Data Management System**
   - Automated result archiving with metadata
   - Calculation database with search capabilities
   - Experiment tracking and comparison tools
   - Integration with electronic lab notebooks

2. **Performance Optimization**
   - GPU acceleration (CUDA/OpenCL if available)
   - Memory management optimization
   - Intelligent job queueing system
   - Cluster computing integration

3. **Advanced MLD Features**
   - Multi-precursor cycle optimization
   - Temperature-dependent kinetics
   - Substrate effect modeling
   - Defect formation and healing

#### **Future Enhancements**
1. **Machine Learning Integration**
   - Predictive models for optimal conditions
   - Automated parameter tuning
   - Pattern recognition in growth behavior

2. **Experimental Integration**
   - Real-time comparison with experimental data
   - QCM/ellipsometry data correlation
   - Process control optimization

---

## 🎨 Visualization Solutions

### **Current Terminal-Based Tools** (Available Now)
```bash
# ASCII structure viewer
python ascii_structure_viewer.py

# Terminal molecular display with bond analysis
python terminal_visualizer.py

# Surface chemistry visualization
python surface_visualizer.py

# MLD cycle progress analysis
python mld_analyzer.py
```

### **Recommended GUI Visualization Tools**

#### **1. VESTA** ⭐ *Best for Crystal/Surface Structures*
- **Download**: https://jp-minerals.org/vesta/
- **Features**: Professional crystal structure visualization, excellent for surfaces
- **Use Case**: Publication-quality figures, crystal analysis
- **Installation**: Download and extract, runs on Windows/Linux

#### **2. Avogadro** ⭐ *Best for Molecular Systems*
```bash
# Linux installation
sudo apt install avogadro

# Windows: Download from https://avogadro.cc/
```
- **Features**: Molecular editor, DFT result visualization, animation
- **Use Case**: Molecular optimization, reaction pathways

#### **3. PyMOL** ⭐ *Professional Molecular Graphics*
- **Download**: https://pymol.org/
- **Features**: Industry standard, excellent rendering, scripting
- **Use Case**: Publication figures, complex molecular systems
- **Note**: Free educational version available

#### **4. ParaView** ⭐ *Advanced Scientific Visualization*
- **Download**: https://www.paraview.org/
- **Features**: Large dataset handling, custom visualizations
- **Use Case**: Complex analysis, custom plotting

#### **5. ASE GUI** *Simple, Built-in*
```python
from ase.gui import gui
gui.view(atoms)  # Quick structure viewing
```

#### **6. Custom QT Application** *Future Development*
- **Plan**: Create dedicated MLD visualization tool
- **Features**: Real-time progress, integrated workflow, surface analysis
- **Timeline**: High priority for next development phase

### **Visualization Workflow Recommendations**
1. **Quick checks**: ASE GUI, terminal visualizers
2. **Analysis**: VESTA (surfaces), Avogadro (molecules)
3. **Publications**: PyMOL, VESTA
4. **Data analysis**: ParaView, custom scripts

---

## 🔄 File Transfer & Synchronization

### **Upload Files to Lab Desktop**
```bash
# Single file
scp file.txt lab-desktop:~/mld-smart-chemistry/

# Multiple files
scp *.py lab-desktop:~/mld-smart-chemistry/scripts/

# Entire directory
scp -r local_folder/ lab-desktop:~/mld-smart-chemistry/
```

### **Download Results from Lab Desktop**
```bash
# Single file
scp lab-desktop:~/mld-smart-chemistry/results.xyz ./

# All results
scp lab-desktop:~/mld-smart-chemistry/*.xyz ./results/

# Entire results directory
scp -r lab-desktop:~/mld-smart-chemistry/results/ ./
```

### **Synchronization (Recommended)**
```bash
# Download all new results
rsync -avz lab-desktop:~/mld-smart-chemistry/results/ ./results/

# Upload modified scripts
rsync -avz ./scripts/ lab-desktop:~/mld-smart-chemistry/scripts/

# Two-way sync (careful!)
rsync -avz --delete ./scripts/ lab-desktop:~/mld-smart-chemistry/scripts/
```

### **Git-Based Workflow** (Recommended)
```bash
# On laptop: commit and push changes
git add .
git commit -m "Update analysis scripts"
git push origin main

# On lab desktop: pull latest changes
git pull origin main
```

---

## 💡 Advanced Tips & Tricks

### **Persistent Sessions**
```bash
# Using screen (recommended)
screen -S mld_simulation
# Run your simulation
# Detach: Ctrl+A, then D
# Reattach: screen -r mld_simulation

# Using tmux (alternative)
tmux new -s mld_simulation
# Detach: Ctrl+B, then D
# Reattach: tmux attach -t mld_simulation
```

### **Monitoring & Alerts**
```bash
# Email notification when job completes
echo "Simulation completed" | mail -s "MLD Job Done" your.email@uw.edu

# System status in one line
echo "Load: $(uptime | cut -d: -f5) | Memory: $(free -h | grep Mem | awk '{print $3"/"$2}') | Disk: $(df -h / | tail -1 | awk '{print $5}')"

# Auto-check simulation every hour
watch -n 3600 'ps aux | grep python && tail -5 *.log'
```

### **Performance Optimization**
```bash
# Check CPU usage by process
ps aux --sort=-%cpu | head -10

# Monitor memory usage
watch -n 5 'free -h && ps aux --sort=-%mem | head -5'

# Check I/O usage
iotop

# System temperature (if available)
sensors
```

### **Backup & Safety**
```bash
# Automatic backup of results
rsync -avz ~/mld-smart-chemistry/results/ /backup/mld_results_$(date +%Y%m%d)/

# Create checkpoint saves
cp important_structure.xyz important_structure_$(date +%Y%m%d_%H%M).xyz.bak

# Git auto-commit results
git add results/ && git commit -m "Auto-save results $(date)" && git push
```

---

## 📞 Emergency Procedures

### **Kill Runaway Processes**
```bash
# Find process ID
ps aux | grep python
ps aux | grep gpaw

# Kill specific process
kill -9 PID

# Kill all Python processes (use with caution!)
pkill -f python

# Kill specific simulation
pkill -f "dft_small_mld"
```

### **Free Up System Resources**
```bash
# Clear system cache
sync && echo 3 | sudo tee /proc/sys/vm/drop_caches

# Find memory-hungry processes
ps aux --sort=-%mem | head -10

# Find large files
find . -size +1G -ls

# Clean temporary files
rm -rf /tmp/* ~/.cache/*
```

### **SSH Connection Recovery**
```bash
# If SSH connection is lost, try:
ssh -p 2222 dreece23@172.25.145.180

# If connection refused, restart SSH (via Remote Desktop):
~/restart_ssh.sh

# Check SSH status
~/check_ssh_status.sh

# Reset Windows port forwarding (PowerShell as Admin):
netsh interface portproxy delete v4tov4 listenport=2222 listenaddress=0.0.0.0
netsh interface portproxy add v4tov4 listenport=2222 connectaddress=172.29.234.10 connectport=2222 listenaddress=0.0.0.0
```

### **Data Recovery**
```bash
# Check for auto-saved files
ls -la *_backup.* *_checkpoint.* *.bak

# Find recent modifications
find . -mtime -1 -name "*.xyz" -o -name "*.log"

# Git recovery
git log --oneline | head -10
git checkout HEAD~1 -- important_file.py
```

---

## 📊 Expected Performance Metrics

### **Typical Runtimes**
- **H2O molecule test**: 30-60 seconds
- **TMA molecule test**: 2-5 minutes  
- **Surface cluster test**: 10-30 minutes
- **H2O+surface test**: 20-60 minutes
- **Full MLD cycle**: 1-3 hours
- **Complete MLD simulation**: 4-12 hours

### **System Requirements Met**
- **Atom limits**: <200 atoms per calculation
- **Memory usage**: 2-8 GB typical
- **Disk space**: 1-5 GB per simulation
- **Network**: Stable VPN connection required

### **Success Indicators**
- ✅ All DFT tests pass
- ✅ SCF convergence achieved
- ✅ Smart caching reduces repeat calculations
- ✅ Progress bars update in real-time
- ✅ SSH connection stable for hours

---

## 📈 Performance Monitoring

### **Real-time Monitoring Commands**
```bash
# CPU, memory, and process monitoring
htop

# Network activity
nethogs

# Disk I/O
iotop

# GPU usage (if applicable)
nvidia-smi

# Temperature monitoring
watch -n 2 sensors
```

### **Log Analysis**
```bash
# Find calculation bottlenecks
grep -i "time\|duration\|elapsed" *.log

# Monitor convergence trends
grep -A 2 -B 2 "scf.*converged" *.out

# Check memory usage trends
grep -i "memory\|mem" *.log
```

---

## 🎯 Quick Reference Card

### **Essential Commands**
```bash
# Connect to lab desktop
ssh -p 2222 dreece23@172.25.145.180

# Check running simulations
ps aux | grep python && jobs

# Monitor progress
tail -f *.log

# Run DFT test
python scripts/dft_simple_test.py

# Start MLD simulation
nohup python scripts/dft_small_mld.py > output.log 2>&1 &
```

### **Emergency Contacts & Resources**
- **UW IT Support**: For VPN and network issues
- **Lab Manager**: For hardware and access issues
- **GitHub Repository**: https://github.com/your-username/mld-smart-chemistry
- **Documentation**: `scripts/README_DFT_WORKFLOW.md`

---

**Last Updated**: July 16, 2025  
**Version**: 2.0  
**Author**: MLD Smart Chemistry Project

---

*Save this cheat sheet for quick reference during remote simulations!*