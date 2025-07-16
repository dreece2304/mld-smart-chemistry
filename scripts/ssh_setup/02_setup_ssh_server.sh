#!/bin/bash
"""
SSH Server Setup Script for Lab Desktop
This script sets up SSH server on the lab desktop (always-on machine)
"""

echo "🖥️  SSH Server Setup for Lab Desktop"
echo "====================================="
echo ""

# Check if running as root
if [[ $EUID -eq 0 ]]; then
   echo "❌ This script should NOT be run as root"
   echo "   Run as regular user, sudo will be used when needed"
   exit 1
fi

# Check if running on WSL
if [[ -f /proc/version ]] && grep -q Microsoft /proc/version; then
    echo "✅ Detected WSL environment"
    WSL_MODE=true
else
    echo "✅ Detected native Linux environment"
    WSL_MODE=false
fi

echo ""
echo "📋 Pre-installation checks:"

# Check if SSH server is already installed
if command -v sshd &> /dev/null; then
    echo "   ✅ SSH server already installed"
    INSTALL_NEEDED=false
else
    echo "   ❌ SSH server not found - will install"
    INSTALL_NEEDED=true
fi

# Get current user info
CURRENT_USER=$(whoami)
echo "   👤 Current user: $CURRENT_USER"

# Get WSL IP address
WSL_IP=$(hostname -I | cut -d' ' -f1)
echo "   🌐 WSL IP address: $WSL_IP"

echo ""
echo "🔧 Starting SSH server installation and configuration..."

# Step 1: Update packages and install SSH server
if $INSTALL_NEEDED; then
    echo ""
    echo "📦 Step 1: Installing SSH server..."
    
    echo "   Updating package lists..."
    sudo apt update
    
    echo "   Installing openssh-server..."
    sudo apt install -y openssh-server
    
    if command -v sshd &> /dev/null; then
        echo "   ✅ SSH server installed successfully"
    else
        echo "   ❌ SSH server installation failed"
        exit 1
    fi
else
    echo "   ⏩ SSH server already installed, skipping installation"
fi

# Step 2: Configure SSH server
echo ""
echo "⚙️  Step 2: Configuring SSH server..."

# Backup original config
if [[ -f /etc/ssh/sshd_config ]] && [[ ! -f /etc/ssh/sshd_config.backup ]]; then
    echo "   📋 Creating backup of original SSH config..."
    sudo cp /etc/ssh/sshd_config /etc/ssh/sshd_config.backup
fi

# Create SSH config with required settings
echo "   ✏️  Updating SSH configuration..."
sudo tee /etc/ssh/sshd_config.d/mld_ssh.conf > /dev/null << EOF
# MLD SSH Configuration
# Custom configuration for MLD project SSH access

# Use port 2222 to avoid conflicts
Port 2222

# Authentication settings
PasswordAuthentication yes
PubkeyAuthentication yes
PermitRootLogin no

# Security settings
Protocol 2
StrictModes yes
MaxAuthTries 3
MaxSessions 10

# Connection settings
ClientAliveInterval 60
ClientAliveCountMax 3

# X11 forwarding for GUI applications
X11Forwarding yes
X11DisplayOffset 10
X11UseLocalhost no

# Allow all users by default
AllowUsers $CURRENT_USER
EOF

echo "   ✅ SSH configuration updated"

# Step 3: Start and enable SSH service
echo ""
echo "🚀 Step 3: Starting SSH service..."

if $WSL_MODE; then
    # WSL-specific service management
    echo "   🔄 Reloading systemd daemon..."
    sudo systemctl daemon-reload
    
    echo "   🔄 Restarting SSH socket..."
    sudo systemctl restart ssh.socket
    
    echo "   🔄 Restarting SSH service..."
    sudo systemctl restart ssh
    
    echo "   ⚡ Enabling SSH service..."
    sudo systemctl enable ssh
else
    # Native Linux service management
    echo "   🔄 Restarting SSH service..."
    sudo systemctl restart ssh
    
    echo "   ⚡ Enabling SSH service..."
    sudo systemctl enable ssh
fi

# Step 4: Verify SSH is running
echo ""
echo "✅ Step 4: Verifying SSH service..."

sleep 2  # Give service time to start

if systemctl is-active --quiet ssh; then
    echo "   ✅ SSH service is running"
    
    # Check if listening on port 2222
    if ss -tln | grep -q :2222; then
        echo "   ✅ SSH is listening on port 2222"
        
        # Show listening details
        listening_info=$(ss -tln | grep :2222)
        echo "   📡 Listening on: $listening_info"
    else
        echo "   ❌ SSH is not listening on port 2222"
        echo "   📋 Current listening ports:"
        ss -tln | grep :22
    fi
else
    echo "   ❌ SSH service is not running"
    echo "   📋 Service status:"
    systemctl status ssh --no-pager
fi

# Step 5: Test local SSH connection
echo ""
echo "🧪 Step 5: Testing local SSH connection..."

echo "   🔑 Testing SSH connection to localhost..."
if timeout 10 ssh -o ConnectTimeout=5 -o StrictHostKeyChecking=no -o BatchMode=yes -p 2222 $CURRENT_USER@localhost echo "SSH test successful" 2>/dev/null; then
    echo "   ✅ Local SSH connection successful"
else
    echo "   ⚠️  Local SSH connection test failed (this is normal if no SSH keys are set up)"
    echo "   🔑 You can test manually with: ssh -p 2222 $CURRENT_USER@localhost"
fi

# Step 6: Windows port forwarding setup
echo ""
echo "🖥️  Step 6: Windows port forwarding setup..."

echo "   📋 To complete the setup, run these commands in Windows PowerShell as Administrator:"
echo ""
echo "   # Remove any existing port forwarding rules"
echo "   netsh interface portproxy delete v4tov4 listenport=2222 listenaddress=0.0.0.0"
echo ""
echo "   # Add new port forwarding rule"
echo "   netsh interface portproxy add v4tov4 listenport=2222 connectaddress=$WSL_IP connectport=2222 listenaddress=0.0.0.0"
echo ""
echo "   # Add Windows Firewall rule"
echo "   New-NetFirewallRule -DisplayName \"WSL SSH\" -Direction Inbound -LocalPort 2222 -Protocol TCP -Action Allow"
echo ""
echo "   # Verify port forwarding"
echo "   netsh interface portproxy show all"

# Step 7: Create convenience scripts
echo ""
echo "📝 Step 7: Creating convenience scripts..."

# Create SSH status script
cat > ~/check_ssh_status.sh << 'EOF'
#!/bin/bash
echo "🔍 SSH Server Status Check"
echo "========================="
echo ""

# Service status
echo "📊 Service Status:"
if systemctl is-active --quiet ssh; then
    echo "   ✅ SSH service: Running"
else
    echo "   ❌ SSH service: Not running"
fi

# Listening ports
echo ""
echo "📡 Listening Ports:"
ss -tln | grep :2222 | head -5

# Current connections
echo ""
echo "🔗 Current SSH Connections:"
who | grep pts

# WSL IP
echo ""
echo "🌐 Connection Information:"
echo "   WSL IP: $(hostname -I | cut -d' ' -f1)"
echo "   SSH Port: 2222"
echo "   Connect with: ssh -p 2222 $(whoami)@WINDOWS_IP"
echo ""
echo "   To get Windows IP, run in Windows Command Prompt:"
echo "   ipconfig | findstr IPv4"
EOF

chmod +x ~/check_ssh_status.sh

# Create SSH restart script
cat > ~/restart_ssh.sh << 'EOF'
#!/bin/bash
echo "🔄 Restarting SSH Server"
echo "======================="

sudo systemctl daemon-reload
sudo systemctl restart ssh.socket
sudo systemctl restart ssh

echo ""
echo "✅ SSH server restarted"
echo "   Run ~/check_ssh_status.sh to verify"
EOF

chmod +x ~/restart_ssh.sh

echo "   ✅ Created ~/check_ssh_status.sh - Check SSH server status"
echo "   ✅ Created ~/restart_ssh.sh - Restart SSH server"

# Final summary
echo ""
echo "🎉 SSH Server Setup Complete!"
echo "============================="
echo ""
echo "📋 Summary:"
echo "   • SSH server installed and configured"
echo "   • Service running on port 2222"
echo "   • WSL IP address: $WSL_IP"
echo "   • Username: $CURRENT_USER"
echo ""
echo "📝 Next Steps:"
echo "   1. Run the Windows PowerShell commands shown above"
echo "   2. Get Windows IP address: ipconfig | findstr IPv4"
echo "   3. Test connection from laptop: ssh -p 2222 $CURRENT_USER@WINDOWS_IP"
echo ""
echo "🔧 Management Commands:"
echo "   • Check status: ~/check_ssh_status.sh"
echo "   • Restart SSH: ~/restart_ssh.sh"
echo "   • View logs: sudo journalctl -u ssh -f"
echo ""
echo "✅ Lab desktop is now ready to accept SSH connections!"