#!/bin/bash
"""
System Information Checker for SSH Setup
Run this script on both laptop and lab desktop to gather connection info
"""

echo "🔍 SSH Setup - System Information Checker"
echo "========================================"
echo ""

# Check if running on Windows or WSL
if [[ -f /proc/version ]] && grep -q Microsoft /proc/version; then
    echo "📋 Environment: WSL (Windows Subsystem for Linux)"
    WSL_MODE=true
else
    echo "📋 Environment: Native Linux"
    WSL_MODE=false
fi

echo ""
echo "🖥️  System Information:"
echo "   Hostname: $(hostname)"
echo "   User: $(whoami)"
echo "   OS: $(uname -a)"

# Get IP addresses
echo ""
echo "🌐 Network Information:"
if $WSL_MODE; then
    echo "   WSL IP Address: $(hostname -I | cut -d' ' -f1)"
    echo "   Windows IP: Run 'ipconfig' in Windows Command Prompt"
else
    echo "   IP Address: $(hostname -I | cut -d' ' -f1)"
fi

# Check SSH availability
echo ""
echo "🔑 SSH Status:"

# Check if SSH client is available
if command -v ssh &> /dev/null; then
    echo "   ✅ SSH Client: Available"
    ssh -V
else
    echo "   ❌ SSH Client: Not available"
    echo "      Install with: sudo apt install openssh-client"
fi

# Check if SSH server is available
if command -v sshd &> /dev/null; then
    echo "   ✅ SSH Server: Available"
    
    # Check if SSH server is running
    if systemctl is-active --quiet ssh; then
        echo "   ✅ SSH Service: Running"
        
        # Check what port SSH is listening on
        listening_ports=$(ss -tln | grep -E ":22|:2222" | head -5)
        if [[ -n "$listening_ports" ]]; then
            echo "   📡 SSH Listening Ports:"
            echo "$listening_ports" | while read line; do
                echo "      $line"
            done
        else
            echo "   ⚠️  SSH Service: Not listening on standard ports"
        fi
    else
        echo "   ❌ SSH Service: Not running"
    fi
else
    echo "   ❌ SSH Server: Not available"
    echo "      Install with: sudo apt install openssh-server"
fi

# Check for existing SSH keys
echo ""
echo "🔐 SSH Key Information:"
if [[ -d ~/.ssh ]]; then
    echo "   📁 SSH Directory: ~/.ssh exists"
    
    # Check for common key types
    for key_type in id_rsa id_ed25519 id_ecdsa; do
        if [[ -f ~/.ssh/$key_type ]]; then
            echo "   🔑 Private Key: ~/.ssh/$key_type exists"
        fi
        if [[ -f ~/.ssh/$key_type.pub ]]; then
            echo "   🔓 Public Key: ~/.ssh/$key_type.pub exists"
        fi
    done
    
    # Check for SSH config
    if [[ -f ~/.ssh/config ]]; then
        echo "   ⚙️  SSH Config: ~/.ssh/config exists"
        echo "      Configured hosts:"
        grep -E "^Host " ~/.ssh/config | sed 's/Host /      - /'
    else
        echo "   ⚙️  SSH Config: No custom config found"
    fi
else
    echo "   📁 SSH Directory: ~/.ssh does not exist"
fi

# Check for Windows SSH (if in WSL)
if $WSL_MODE; then
    echo ""
    echo "🖥️  Windows SSH Information:"
    echo "   To check Windows SSH client:"
    echo "   - Open Windows PowerShell"
    echo "   - Run: Get-WindowsCapability -Online | Where-Object Name -like 'OpenSSH.Client*'"
    echo "   - Run: ssh -V"
    echo ""
    echo "   To check Windows OpenSSH server:"
    echo "   - Run: Get-WindowsCapability -Online | Where-Object Name -like 'OpenSSH.Server*'"
fi

# Check firewall status (if available)
echo ""
echo "🛡️  Firewall Information:"
if command -v ufw &> /dev/null; then
    echo "   UFW Status: $(sudo ufw status | head -1)"
elif command -v iptables &> /dev/null; then
    echo "   iptables: Available (check manually: sudo iptables -L)"
else
    echo "   Firewall: No common firewall tools found"
fi

if $WSL_MODE; then
    echo "   Windows Firewall: Check Windows Defender Firewall"
fi

# Check for MongoDB (potential conflict)
echo ""
echo "🗄️  Database Services (Potential Conflicts):"
if command -v mongod &> /dev/null; then
    echo "   MongoDB: Available"
    if pgrep mongod &> /dev/null; then
        echo "      Status: Running"
        echo "      Ports: $(ss -tln | grep :27017 || echo 'Not listening on 27017')"
    else
        echo "      Status: Not running"
    fi
else
    echo "   MongoDB: Not installed"
fi

# Summary and recommendations
echo ""
echo "📝 Summary:"
echo "   Run this script on both machines to compare configurations"
echo "   Note down the IP addresses for SSH connection setup"
echo ""
echo "💡 Next Steps:"
echo "   1. Run 'ipconfig' in Windows to get the actual IP address"
echo "   2. If setting up SSH server, run: ./02_setup_ssh_server.sh"
echo "   3. If setting up SSH client, run: ./03_setup_ssh_client.sh"
echo ""
echo "✅ System information check complete!"
echo "   Save this output for reference during SSH setup"