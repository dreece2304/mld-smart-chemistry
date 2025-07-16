# SSH Setup Guide for MLD Project

This directory contains scripts to set up SSH connection between your laptop and lab desktop for running MLD simulations remotely.

## Architecture

```
Your Laptop (Client)          Lab Desktop (Server)
┌─────────────────┐          ┌──────────────────┐
│                 │   SSH    │  Windows         │
│  Windows/WSL    │ -------> │  └── WSL/Ubuntu  │
│  (You type here)│   2222   │      (Runs code) │
└─────────────────┘          └──────────────────┘
```

## Quick Setup

### 1. Check System Information (Both Machines)
```bash
./01_check_system_info.sh
```
Run this on both laptop and lab desktop to gather network and SSH information.

### 2. Setup SSH Server (Lab Desktop Only)
```bash
./02_setup_ssh_server.sh
```
Run this on the lab desktop to install and configure SSH server.

### 3. Setup SSH Client (Laptop Only)
```bash
./03_setup_ssh_client.sh
```
Run this on your laptop to configure SSH client and generate keys.

### 4. Test Connection
```bash
./04_test_connection.sh
```
Test the SSH connection and diagnose any issues.

## Files Description

- **01_check_system_info.sh** - Gathers system and network information
- **02_setup_ssh_server.sh** - Sets up SSH server on lab desktop
- **03_setup_ssh_client.sh** - Configures SSH client on laptop
- **04_test_connection.sh** - Tests connection and provides troubleshooting
- **README.md** - This file

## Expected Network Configuration

Based on your setup:
- **Lab Desktop Windows IP**: 172.25.145.180
- **Lab Desktop WSL IP**: 172.29.234.10  
- **Laptop VPN IP**: 10.102.26.144
- **SSH Port**: 2222 (to avoid conflicts)

## Connection Commands

After setup, connect with:
```bash
# Using SSH config alias
ssh lab-desktop

# Manual connection
ssh -p 2222 dreece23@172.25.145.180
```

## File Transfer

```bash
# Upload files to lab desktop
scp file.txt lab-desktop:~/MLD_ASE_GPAW/

# Download results from lab desktop
scp lab-desktop:~/MLD_ASE_GPAW/results.txt ./
```

## MLD Workflow

1. **Code on laptop** → Edit files locally
2. **Upload to lab desktop** → Transfer files via SSH/git
3. **Run simulations** → SSH to lab desktop, start calculations
4. **Disconnect safely** → Simulations continue running
5. **Reconnect later** → Check results and analyze data

### Running Long Simulations

```bash
# Connect to lab desktop
ssh lab-desktop

# Start simulation in background
nohup python scripts/dft_small_mld.py > output.log 2>&1 &

# Disconnect (simulation continues)
exit

# Later, reconnect and check progress
ssh lab-desktop
tail -f output.log
```

## Troubleshooting

If connection fails:

1. **Check SSH server status** (lab desktop):
   ```bash
   sudo systemctl status ssh
   ss -tln | grep :2222
   ```

2. **Check Windows port forwarding** (lab desktop PowerShell as Admin):
   ```powershell
   netsh interface portproxy show all
   ```

3. **Test network connectivity** (laptop):
   ```bash
   ping 172.25.145.180
   telnet 172.25.145.180 2222
   ```

4. **Check firewall rules** (both machines):
   - Windows Defender Firewall
   - Corporate/university firewall
   - VPN restrictions

## Security Notes

- SSH server configured on non-standard port (2222)
- Password and key-based authentication enabled
- Root login disabled
- X11 forwarding enabled for GUI applications
- Connection keep-alive configured

## VPN Considerations

Your setup uses UW VPN (Big-IP Edge):
- Both machines must be connected to UW network/VPN
- SSH traffic should be allowed through VPN
- Check with IT if SSH connections are blocked

## Support

If you encounter issues:
1. Run the test script: `./04_test_connection.sh`
2. Check the troubleshooting section above
3. Verify all setup steps were completed
4. Contact IT support for network-related issues