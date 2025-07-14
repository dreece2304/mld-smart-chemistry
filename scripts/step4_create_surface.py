#!/usr/bin/env python3
"""
Step 4: Create simple silicon surface for MLD
Keep it simple - just basic Si(100) with some OH groups
"""

from ase import Atoms
from ase.build import surface, add_adsorbate
from ase.io import write
import numpy as np

def create_simple_surface():
    """Create simple Si(100) surface with hydroxyl groups"""
    
    print("🏔️ Step 4: Creating Simple Silicon Surface")
    print("=" * 50)
    
    # Create Si(100) surface - small and simple
    print("📐 Creating Si(100) surface...")
    si_surface = surface('Si', (1, 0, 0), 3, vacuum=10.0)  # 3 layers, 10Å vacuum
    
    print(f"   Surface atoms: {len(si_surface)}")
    print(f"   Surface size: {si_surface.cell[0,0]:.1f} x {si_surface.cell[1,1]:.1f} Å")
    
    # Add some hydroxyl groups (OH) to make it realistic for MLD
    print("🧪 Adding hydroxyl groups (OH)...")
    
    # Find top layer Si atoms
    positions = si_surface.get_positions()
    z_coords = positions[:, 2]
    z_max = z_coords.max()
    
    # Get top Si atoms (within 1Å of top)
    top_si_indices = []
    for i, pos in enumerate(positions):
        if pos[2] > z_max - 1.0:
            top_si_indices.append(i)
    
    print(f"   Top Si atoms: {len(top_si_indices)}")
    
    # Add OH groups to half the top Si atoms (50% coverage)
    np.random.seed(42)  # Reproducible
    oh_sites = np.random.choice(top_si_indices, size=len(top_si_indices)//2, replace=False)
    
    oh_positions = []
    oh_symbols = []
    
    for si_idx in oh_sites:
        si_pos = positions[si_idx]
        # Add O above Si (1.6 Å Si-O bond)
        o_pos = si_pos + [0, 0, 1.6]
        # Add H above O (0.96 Å O-H bond)  
        h_pos = o_pos + [0, 0, 0.96]
        
        oh_positions.extend([o_pos, h_pos])
        oh_symbols.extend(['O', 'H'])
    
    print(f"   Adding {len(oh_sites)} OH groups...")
    
    # Create complete surface with OH groups
    all_positions = np.vstack([positions, oh_positions])
    all_symbols = list(si_surface.get_chemical_symbols()) + oh_symbols
    
    # Create final surface
    surface_with_oh = Atoms(
        symbols=all_symbols,
        positions=all_positions,
        cell=si_surface.get_cell(),
        pbc=si_surface.get_pbc()
    )
    
    # Save surface
    write('simple_surface.xyz', surface_with_oh)
    
    print(f"✅ Surface created successfully!")
    print(f"   Total atoms: {len(surface_with_oh)}")
    print(f"   Si atoms: {len([s for s in all_symbols if s == 'Si'])}")
    print(f"   OH groups: {len(oh_sites)}")
    print(f"   Saved: simple_surface.xyz")
    
    return surface_with_oh

def main():
    """Create simple silicon surface"""
    surface_atoms = create_simple_surface()
    
    print(f"\n🎉 Simple silicon surface ready for MLD!")
    print(f"Next: Visualize with python visualize_molecules.py")
    
    return 0

if __name__ == "__main__":
    exit(main())