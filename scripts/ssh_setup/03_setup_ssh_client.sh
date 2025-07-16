#!/bin/bash
"""
SSH Client Setup Script for Laptop
This script sets up SSH client configuration on your laptop to connect to lab desktop
"""

echo "💻 SSH Client Setup for Laptop"
echo "=============================="
echo ""

# Get current user info
CURRENT_USER=$(whoami)
echo "👤 Current user: $CURRENT_USER"

# Check if SSH client is available
if command -v ssh &> /dev/null; then
    echo "✅ SSH client is available"
    echo "   Version: $(ssh -V 2>&1)"
else
    echo "❌ SSH client not found"
    echo "   Installing SSH client..."
    sudo apt update
    sudo apt install -y openssh-client
fi

echo ""
echo "🔑 SSH Key Setup"

# Create SSH directory
mkdir -p ~/.ssh
chmod 700 ~/.ssh

# Check for existing SSH keys
echo "   Checking for existing SSH keys..."
has_keys=false

for key_type in id_ed25519 id_rsa; do
    if [[ -f ~/.ssh/$key_type ]]; then
        echo "   ✅ Found existing key: ~/.ssh/$key_type"
        has_keys=true
    fi
done

if ! $has_keys; then
    echo "   🔧 Creating new SSH key pair..."
    echo "   Generating ED25519 key (recommended)..."
    
    # Generate SSH key
    ssh-keygen -t ed25519 -C "$CURRENT_USER@laptop" -f ~/.ssh/id_ed25519 -N ""
    
    echo "   ✅ SSH key generated: ~/.ssh/id_ed25519"
fi

echo ""
echo "⚙️  SSH Configuration Setup"

# Get lab desktop IP from user
echo "   Enter lab desktop Windows IP address (e.g., 172.25.145.180):"
read -p "   Lab Desktop IP: " LAB_IP

echo "   Enter lab desktop WSL username (e.g., dreece23):"
read -p "   Lab Username: " LAB_USER

# Create SSH config
cat > ~/.ssh/config << EOF
# MLD Lab Desktop Connection
Host lab-desktop
    HostName $LAB_IP
    Port 2222
    User $LAB_USER
    IdentityFile ~/.ssh/id_ed25519
    ServerAliveInterval 60
    ServerAliveCountMax 3
    ForwardX11 yes
    ForwardX11Trusted yes

# Backup connection using standard port
Host lab-desktop-standard
    HostName $LAB_IP
    Port 22
    User $LAB_USER
    IdentityFile ~/.ssh/id_ed25519
EOF

chmod 600 ~/.ssh/config

echo "   ✅ SSH config created: ~/.ssh/config"

echo ""
echo "🔐 SSH Key Copy Instructions"
echo ""
echo "   To enable key-based authentication, copy your public key to the lab desktop:"
echo ""
echo "   1. Display your public key:"
echo "      cat ~/.ssh/id_ed25519.pub"
echo ""
echo "   2. Copy the entire output"
echo ""
echo "   3. On lab desktop, run:"
echo "      mkdir -p ~/.ssh"
echo "      echo 'PASTE_YOUR_PUBLIC_KEY_HERE' >> ~/.ssh/authorized_keys"
echo "      chmod 600 ~/.ssh/authorized_keys"
echo ""

echo "📋 Connection Commands"
echo ""
echo "   Quick connect to lab desktop:"
echo "      ssh lab-desktop"
echo ""
echo "   Manual connection:"
echo "      ssh -p 2222 $LAB_USER@$LAB_IP"
echo ""
echo "   Copy files to lab desktop:"
echo "      scp file.txt lab-desktop:~/MLD_ASE_GPAW/"
echo ""
echo "   Copy files from lab desktop:"
echo "      scp lab-desktop:~/MLD_ASE_GPAW/results.txt ./"
echo ""

echo "🧪 Connection Test"
echo ""
echo "   After setting up SSH server on lab desktop, test with:"
echo "      ssh lab-desktop"
echo ""
echo "   Or manually:"
echo "      ssh -p 2222 $LAB_USER@$LAB_IP"
echo ""

echo "✅ SSH client setup complete!"
echo ""
echo "💡 Next Steps:"
echo "   1. Set up SSH server on lab desktop: ./02_setup_ssh_server.sh"
echo "   2. Copy SSH public key to lab desktop (see instructions above)"
echo "   3. Test connection: ssh lab-desktop"
echo ""

# Show public key for easy copying
echo "🔑 Your SSH Public Key (copy this to lab desktop):"
echo "=================================================="
if [[ -f ~/.ssh/id_ed25519.pub ]]; then
    cat ~/.ssh/id_ed25519.pub
elif [[ -f ~/.ssh/id_rsa.pub ]]; then
    cat ~/.ssh/id_rsa.pub
else
    echo "   No public key found"
fi
echo "=================================================="