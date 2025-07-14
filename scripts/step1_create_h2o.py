#!/usr/bin/env python3
"""
Step 1: Create a simple H2O molecule structure
"""

from ase import Atoms
from ase.io import write

# Create H2O molecule with realistic bond lengths
h2o = Atoms('H2O', positions=[
    [0.0, 0.0, 0.0],      # O
    [0.96, 0.0, 0.0],     # H1  
    [-0.24, 0.93, 0.0]    # H2 (104.5° angle)
])

# Center with vacuum
h2o.center(vacuum=8.0)

# Save to file
write('structures/h2o_molecule.xyz', h2o)

print("✅ Created H2O molecule")
print(f"   Atoms: {len(h2o)}")
print(f"   Saved to: structures/h2o_molecule.xyz")