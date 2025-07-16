#!/bin/bash
"""
SSH Connection Tester
Tests SSH connection and diagnoses issues
"""

echo "🧪 SSH Connection Tester"
echo "========================"
echo ""

# Get connection details
echo "📋 Connection Details"
echo "   Enter lab desktop Windows IP address:"
read -p "   Lab Desktop IP: " LAB_IP

echo "   Enter lab desktop WSL username:"
read -p "   Lab Username: " LAB_USER

echo "   Enter SSH port (default 2222):"
read -p "   SSH Port [2222]: " SSH_PORT
SSH_PORT=${SSH_PORT:-2222}

echo ""
echo "🔍 Testing Connection to $LAB_USER@$LAB_IP:$SSH_PORT"
echo "================================================="

# Test 1: Basic network connectivity
echo ""
echo "Test 1: Network Connectivity"
echo "   Testing ping to $LAB_IP..."
if ping -c 3 $LAB_IP &> /dev/null; then
    echo "   ✅ Ping successful"
else
    echo "   ❌ Ping failed (this may be normal if ping is blocked)"
fi

# Test 2: Port accessibility
echo ""
echo "Test 2: Port Accessibility"
echo "   Testing if port $SSH_PORT is open..."
if timeout 10 bash -c "cat < /dev/null > /dev/tcp/$LAB_IP/$SSH_PORT" 2>/dev/null; then
    echo "   ✅ Port $SSH_PORT is accessible"
else
    echo "   ❌ Port $SSH_PORT is not accessible"
    echo "      This usually means:"
    echo "      - SSH server is not running"
    echo "      - Port forwarding not set up"
    echo "      - Firewall blocking connection"
fi

# Test 3: SSH connection with timeout
echo ""
echo "Test 3: SSH Connection"
echo "   Testing SSH connection with timeout..."
if timeout 10 ssh -o ConnectTimeout=5 -o BatchMode=yes -o StrictHostKeyChecking=no -p $SSH_PORT $LAB_USER@$LAB_IP echo "SSH connection successful" 2>/dev/null; then
    echo "   ✅ SSH connection successful"
else
    echo "   ❌ SSH connection failed"
    echo "      Trying verbose connection for diagnosis..."
    
    # Test with verbose output
    echo ""
    echo "   Verbose SSH output:"
    timeout 10 ssh -vvv -o ConnectTimeout=5 -o BatchMode=yes -o StrictHostKeyChecking=no -p $SSH_PORT $LAB_USER@$LAB_IP echo "test" 2>&1 | head -20
fi

# Test 4: SSH with password prompt
echo ""
echo "Test 4: Interactive SSH Connection"
echo "   Attempting interactive SSH connection..."
echo "   (This will prompt for password if connection works)"
echo ""

ssh -o ConnectTimeout=10 -o StrictHostKeyChecking=no -p $SSH_PORT $LAB_USER@$LAB_IP

# Troubleshooting tips
echo ""
echo "🔧 Troubleshooting Tips"
echo "======================"
echo ""
echo "If connection failed, check:"
echo ""
echo "1. SSH Server on Lab Desktop:"
echo "   - Run: ./02_setup_ssh_server.sh"
echo "   - Check: sudo systemctl status ssh"
echo "   - Verify: ss -tln | grep :$SSH_PORT"
echo ""
echo "2. Windows Port Forwarding:"
echo "   - Run in PowerShell as Admin:"
echo "   - netsh interface portproxy show all"
echo "   - netsh interface portproxy add v4tov4 listenport=$SSH_PORT connectaddress=WSL_IP connectport=$SSH_PORT listenaddress=0.0.0.0"
echo ""
echo "3. Firewall Rules:"
echo "   - Windows: New-NetFirewallRule -DisplayName \"WSL SSH\" -Direction Inbound -LocalPort $SSH_PORT -Protocol TCP -Action Allow"
echo "   - Check if corporate firewall blocks SSH"
echo ""
echo "4. VPN Connection:"
echo "   - Ensure both machines are on same network/VPN"
echo "   - Check routing: ip route"
echo ""
echo "5. Network Diagnostics:"
echo "   - Test different ports: telnet $LAB_IP $SSH_PORT"
echo "   - Check DNS resolution: nslookup $LAB_IP"
echo ""
echo "✅ Connection test complete!"