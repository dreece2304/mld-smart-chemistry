#!/bin/bash
# Setup GitHub SSH access for new machine

set -e

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}GitHub SSH Setup for MLD Repository${NC}"
echo "===================================="

# Check if SSH key already exists
if [ -f ~/.ssh/id_ed25519 ]; then
    echo -e "${YELLOW}SSH key already exists. Using existing key.${NC}"
else
    echo -e "${GREEN}Generating new SSH key...${NC}"
    read -p "Enter your GitHub email: " github_email
    ssh-keygen -t ed25519 -C "$github_email" -f ~/.ssh/id_ed25519
fi

# Start SSH agent
echo -e "${GREEN}Starting SSH agent...${NC}"
eval "$(ssh-agent -s)"

# Add key to agent
ssh-add ~/.ssh/id_ed25519

# Display public key
echo -e "${BLUE}
Your SSH public key:
====================
${NC}"
cat ~/.ssh/id_ed25519.pub

echo -e "${YELLOW}
Next steps:
1. Copy the SSH key above
2. Go to GitHub.com → Settings → SSH and GPG keys → New SSH key
3. Paste the key and save
4. Run: ssh -T git@github.com
5. Clone the repository: git clone git@github.com:yourusername/mld-smart-chemistry.git
${NC}"

# Create SSH config for GitHub
if ! grep -q "Host github.com" ~/.ssh/config 2>/dev/null; then
    echo -e "${GREEN}Creating SSH config for GitHub...${NC}"
    mkdir -p ~/.ssh
    cat >> ~/.ssh/config << 'EOF'

Host github.com
    HostName github.com
    User git
    IdentityFile ~/.ssh/id_ed25519
    AddKeysToAgent yes
EOF
    chmod 600 ~/.ssh/config
fi

echo -e "${GREEN}SSH setup complete!${NC}"