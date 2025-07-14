# Updating WSL2 and Ubuntu on Target Machine

## 1. Update WSL2 (Windows Side)

### In Windows PowerShell (as Administrator):
```powershell
# Check current WSL version
wsl --version

# Update WSL to latest version
wsl --update

# Shutdown WSL to apply updates
wsl --shutdown

# List installed distributions
wsl --list --verbose
```

## 2. Install/Update Ubuntu Distribution

### If Ubuntu not installed:
```powershell
# Install latest Ubuntu LTS (24.04)
wsl --install -d Ubuntu-24.04

# Or install Ubuntu 22.04 if preferred
wsl --install -d Ubuntu-22.04
```

### If Ubuntu already installed:
```powershell
# Set default WSL version to 2
wsl --set-default-version 2

# Convert existing Ubuntu to WSL2 if needed
wsl --set-version Ubuntu 2
```

## 3. Update Ubuntu System (Inside WSL)

### First time setup:
```bash
# Update package lists
sudo apt update

# Upgrade all packages
sudo apt upgrade -y

# Upgrade distribution (if major version update available)
sudo apt dist-upgrade -y

# Remove unnecessary packages
sudo apt autoremove -y

# Check Ubuntu version
lsb_release -a
```

### Enable systemd (recommended for WSL2):
```bash
# Edit wsl.conf
sudo nano /etc/wsl.conf

# Add these lines:
[boot]
systemd=true

# Save and exit, then in PowerShell:
wsl --shutdown
# Restart WSL
```

## 4. Optimize WSL2 Performance

### Create .wslconfig in Windows (C:\Users\YourUsername\.wslconfig):
```ini
[wsl2]
memory=8GB              # Limit memory (adjust based on system)
processors=4            # Number of processors
localhostForwarding=true
kernelCommandLine=vsyscall=emulate
swap=4GB
swapFile=C:\\temp\\wsl-swap.vhdx
```

### In Ubuntu, optimize for computational work:
```bash
# Increase file watchers
echo fs.inotify.max_user_watches=524288 | sudo tee -a /etc/sysctl.conf

# Apply changes
sudo sysctl -p
```

## 5. Git and SSH Setup

### Configure Git:
```bash
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"
```

### Setup SSH for GitHub:
```bash
# Generate SSH key
ssh-keygen -t ed25519 -C "your.email@example.com"

# Start ssh-agent
eval "$(ssh-agent -s)"

# Add key to agent
ssh-add ~/.ssh/id_ed25519

# Copy public key to clipboard
cat ~/.ssh/id_ed25519.pub
# Add this key to GitHub: Settings → SSH and GPG keys → New SSH key

# Test connection
ssh -T git@github.com
```

## 6. Clone and Setup MLD Repository

```bash
# Clone via SSH
git clone git@github.com:yourusername/mld-smart-chemistry.git
cd mld-smart-chemistry

# Run installation
./install/INSTALLATION_STEPS.md  # Follow manual steps first
./install/install_conda_and_packages.sh  # Then automated installation
```

## Troubleshooting

### WSL2 Network Issues:
```bash
# Reset network
sudo rm /etc/resolv.conf
sudo bash -c 'echo "nameserver 8.8.8.8" > /etc/resolv.conf'
sudo bash -c 'echo "[network]" > /etc/wsl.conf'
sudo bash -c 'echo "generateResolvConf = false" >> /etc/wsl.conf'
```

### Memory Issues:
Adjust memory in .wslconfig file based on available RAM

### Performance Issues:
- Move large files to Linux filesystem (/home/username/)
- Avoid accessing Windows files from WSL frequently
- Use WSL2 (not WSL1) for better performance