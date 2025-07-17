# SSH & Remote Execution Guide

## Overview

This guide covers remote execution of MLD simulations on the lab desktop computer via SSH, including setup, monitoring, and troubleshooting.

---

## Connection Setup

### Network Requirements
- **VPN**: Must be connected to UW VPN (Big-IP Edge)
- **Lab Desktop IP**: `172.25.145.180`
- **SSH Port**: `2222`

### Basic SSH Connection
```bash
# Direct connection
ssh -p 2222 dreece23@172.25.145.180

# With X11 forwarding (for GUI applications)
ssh -X -p 2222 dreece23@172.25.145.180
```

### SSH Config Setup
Create `~/.ssh/config` on your local machine:
```bash
Host lab-desktop
    HostName 172.25.145.180
    Port 2222
    User dreece23
    ForwardX11 yes
    ServerAliveInterval 60
    ServerAliveCountMax 3
```

Then connect with:
```bash
ssh lab-desktop
```

---

## Remote Execution Workflow

### 1. Environment Setup (First Time)
```bash
# Connect to lab desktop
ssh lab-desktop

# Navigate to project
cd ~/mld-smart-chemistry

# Activate environment
conda activate mld_modeling

# Verify setup
python tests/test_suite.py --level 1
```

### 2. Background Execution
```bash
# Start long-running calculation in background
nohup python src/simulation/dft_small_mld.py > dft_output.log 2>&1 &

# Get process ID
echo $! > simulation.pid

# Logout safely (calculation continues)
exit
```

### 3. Progress Monitoring
```bash
# Reconnect later
ssh lab-desktop

# Check if process is running
ps aux | grep python

# Or check specific PID
kill -0 $(cat simulation.pid) && echo "Running" || echo "Stopped"

# Monitor progress
python scripts/monitoring/monitor_dft_progress.py dft_output.log

# Or simple log monitoring
tail -f dft_output.log
```

---

## Monitoring & Management

### Check Running Processes
```bash
# All Python processes
ps aux | grep python

# All GPAW processes
ps aux | grep gpaw

# Background jobs in current session
jobs

# System resources
htop
```

### Real-time Monitoring
```bash
# Live progress monitoring
python scripts/monitoring/monitor_dft_progress.py dft_output.log

# Monitor with refresh rate
python scripts/monitoring/monitor_dft_progress.py dft_output.log --refresh 0.5

# System resource monitoring
watch -n 1 'free -h && echo "CPU:" && uptime'
```

### Process Management
```bash
# Check specific process
ps -p $(cat simulation.pid)

# Stop process gracefully
kill -TERM $(cat simulation.pid)

# Force stop (last resort)
kill -9 $(cat simulation.pid)

# Check exit status
echo $?
```

---

## File Management

### Transfer Files
```bash
# Copy results to local machine
scp -P 2222 dreece23@172.25.145.180:~/mld-smart-chemistry/data/outputs/structures/* ./local_results/

# Copy entire results directory
scp -r -P 2222 dreece23@172.25.145.180:~/mld-smart-chemistry/data/outputs/ ./

# Upload input files
scp -P 2222 ./local_input.xyz dreece23@172.25.145.180:~/mld-smart-chemistry/data/structures/
```

### Synchronization
```bash
# Sync project directory (from local machine)
rsync -avz -e "ssh -p 2222" ./mld-smart-chemistry/ dreece23@172.25.145.180:~/mld-smart-chemistry/

# Sync results back (from local machine)
rsync -avz -e "ssh -p 2222" dreece23@172.25.145.180:~/mld-smart-chemistry/data/outputs/ ./results/
```

---

## Typical Workflow Examples

### Quick Test Run
```bash
# Connect and test
ssh lab-desktop
cd ~/mld-smart-chemistry
conda activate mld_modeling
python examples/quick_run.py --test
```

### Full MLD Simulation
```bash
# Start background simulation
ssh lab-desktop
cd ~/mld-smart-chemistry
conda activate mld_modeling
nohup python src/simulation/dft_small_mld.py > mld_$(date +%Y%m%d_%H%M%S).log 2>&1 &
echo $! > mld_simulation.pid

# Monitor progress (can disconnect/reconnect)
python scripts/monitoring/monitor_dft_progress.py mld_$(date +%Y%m%d_%H%M%S).log
```

### Surface Optimization
```bash
# Optimize surface with restart capability
ssh lab-desktop
cd ~/mld-smart-chemistry
conda activate mld_modeling
nohup python src/simulation/optimize_surface_only.py small_si100_2x2_4L.xyz --restart > surface_opt.log 2>&1 &
```

---

## Connection Troubleshooting

### Test Connection Diagnostics
```bash
# Basic connectivity test
ping 172.25.145.180

# Test SSH port
telnet 172.25.145.180 2222

# Verbose SSH connection
ssh -v -p 2222 dreece23@172.25.145.180

# Test with timeout
timeout 10 ssh -p 2222 dreece23@172.25.145.180 echo "Connected"
```

### Common Issues & Solutions

#### Connection Refused
```bash
# Check if lab desktop is on
ping 172.25.145.180

# Check Windows port forwarding (on lab desktop, PowerShell as Admin)
netsh interface portproxy show all

# Expected output should show:
# 172.25.145.180:2222 -> 127.0.0.1:22

# Restart SSH service on lab desktop
# Run these on lab desktop via Remote Desktop:
~/restart_ssh.sh
~/check_ssh_status.sh
```

#### VPN Issues
```bash
# Check VPN connection
ping 172.25.145.180

# If fails, reconnect to UW VPN (Big-IP Edge)
# Verify IP is in UW range
curl ifconfig.me
```

#### Slow Connection
```bash
# Test connection speed
ssh -p 2222 dreece23@172.25.145.180 "time dd if=/dev/zero of=/tmp/test bs=1M count=100"

# Use compression
ssh -C -p 2222 dreece23@172.25.145.180

# Add to ~/.ssh/config
Compression yes
```

---

## Performance Optimization

### For Large Calculations
```bash
# Use screen/tmux for persistent sessions
ssh lab-desktop
screen -S mld_simulation
# Or
tmux new-session -s mld_simulation

# Run calculation
python src/simulation/dft_small_mld.py

# Detach: Ctrl+A, D (screen) or Ctrl+B, D (tmux)
# Reattach later: screen -r mld_simulation or tmux attach -t mld_simulation
```

### System Resource Management
```bash
# Check available resources
free -h
df -h
nproc  # Number of CPU cores

# Monitor during calculation
htop

# Nice priority for long calculations
nice -n 10 python src/simulation/dft_small_mld.py
```

---

## Automation Scripts

### Connection Test Script
```bash
#!/bin/bash
# test_lab_connection.sh
echo "Testing lab desktop connection..."
if timeout 5 ssh -p 2222 dreece23@172.25.145.180 echo "Connection successful"; then
    echo "✅ Lab desktop accessible"
else
    echo "❌ Lab desktop not accessible"
    echo "Check VPN connection and try again"
fi
```

### Simulation Status Script
```bash
#!/bin/bash
# check_simulation_status.sh
PID_FILE="simulation.pid"
if [ -f "$PID_FILE" ]; then
    PID=$(cat $PID_FILE)
    if kill -0 $PID 2>/dev/null; then
        echo "✅ Simulation running (PID: $PID)"
        ps -p $PID -o pid,etime,cmd
    else
        echo "❌ Simulation not running"
    fi
else
    echo "❌ No PID file found"
fi
```

---

## Security Best Practices

### SSH Key Setup
```bash
# Generate SSH key pair (on local machine)
ssh-keygen -t rsa -b 4096 -f ~/.ssh/lab_desktop_key

# Copy public key to lab desktop
ssh-copy-id -i ~/.ssh/lab_desktop_key.pub -p 2222 dreece23@172.25.145.180

# Add to SSH config
Host lab-desktop
    HostName 172.25.145.180
    Port 2222
    User dreece23
    IdentityFile ~/.ssh/lab_desktop_key
```

### Connection Security
```bash
# Use SSH config for consistent settings
# Add to ~/.ssh/config
Host lab-desktop
    HostName 172.25.145.180
    Port 2222
    User dreece23
    ForwardX11 yes
    ServerAliveInterval 60
    ServerAliveCountMax 3
    PreferredAuthentications publickey
```

---

## Emergency Procedures

### If Connection is Lost
```bash
# Check if processes are still running
ssh lab-desktop "ps aux | grep python"

# Kill runaway processes if needed
ssh lab-desktop "pkill -f 'python.*dft'"

# Check system status
ssh lab-desktop "uptime && free -h"
```

### System Recovery
```bash
# If lab desktop is unresponsive
# Contact lab admin for remote desktop access
# Or use physical access to restart services

# On lab desktop (via Remote Desktop):
~/restart_ssh.sh
~/check_ssh_status.sh
```

---

## Quick Reference

### Essential Commands
```bash
# Connect
ssh lab-desktop

# Background execution
nohup python script.py > output.log 2>&1 &

# Monitor progress
python scripts/monitoring/monitor_dft_progress.py output.log

# Check processes
ps aux | grep python

# Transfer files
scp -P 2222 file.xyz dreece23@172.25.145.180:~/
```

### Monitoring Commands
```bash
# System resources
htop
free -h
df -h

# Process monitoring
ps aux | grep python
jobs
pgrep -f gpaw

# Log monitoring
tail -f output.log
grep -i error output.log
```

---

**Last Updated**: July 2025