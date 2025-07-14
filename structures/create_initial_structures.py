#!/usr/bin/env python3
"""
Create initial structures for MLD simulation
"""

import numpy as np
from ase import Atoms
from ase.io import write
from ase.build import surface, molecule, add_adsorbate
from ase.visualize import view

def create_tma_molecule():
    """Create trimethylaluminum (TMA) molecule"""
    # Start with aluminum at origin
    tma = Atoms('Al', positions=[[0, 0, 0]])
    
    # Add three methyl groups in tetrahedral arrangement
    # Al-C bond length ~2.0 Å
    methyl_positions = [
        [0, 0, 2.0],           # Above
        [1.73, 0, -1.0],       # 120° separation
        [-1.73, 0, -1.0]       # 120° separation
    ]
    
    for i, pos in enumerate(methyl_positions):
        # Add carbon
        tma.append(Atoms('C', positions=[pos]))
        
        # Add three hydrogens around each carbon
        c_pos = np.array(pos)
        c_direction = c_pos / np.linalg.norm(c_pos)  # Direction from Al to C
        
        # Create orthogonal vectors for H placement
        if abs(c_direction[2]) < 0.9:
            v1 = np.cross(c_direction, [0, 0, 1])
        else:
            v1 = np.cross(c_direction, [1, 0, 0])
        v1 = v1 / np.linalg.norm(v1)
        v2 = np.cross(c_direction, v1)
        v2 = v2 / np.linalg.norm(v2)
        
        # C-H bond length ~1.1 Å
        h_distance = 1.1
        for j in range(3):
            angle = j * 2 * np.pi / 3
            h_offset = h_distance * (np.cos(angle) * v1 + np.sin(angle) * v2)
            h_pos = c_pos + h_offset + 0.3 * c_direction  # Slight outward
            tma.append(Atoms('H', positions=[h_pos]))
    
    # Center molecule
    tma.center()
    return tma

def create_butyne_diol():
    """Create 2-butyne-1,4-diol molecule HO-CH2-C≡C-CH2-OH"""
    # Start with acetylene backbone
    positions = []
    symbols = []
    
    # C≡C triple bond, length ~1.2 Å
    c1_pos = np.array([-0.6, 0, 0])
    c2_pos = np.array([0.6, 0, 0])
    
    positions.extend([c1_pos, c2_pos])
    symbols.extend(['C', 'C'])
    
    # CH2 groups, C-C bond ~1.5 Å
    ch2_left = c1_pos + np.array([-1.5, 0, 0])
    ch2_right = c2_pos + np.array([1.5, 0, 0])
    
    positions.extend([ch2_left, ch2_right])
    symbols.extend(['C', 'C'])
    
    # OH groups, C-O bond ~1.4 Å
    oh_left = ch2_left + np.array([-1.4, 0, 0])
    oh_right = ch2_right + np.array([1.4, 0, 0])
    
    positions.extend([oh_left, oh_right])
    symbols.extend(['O', 'O'])
    
    # H on OH groups, O-H bond ~0.96 Å
    h_oh_left = oh_left + np.array([-0.96, 0, 0])
    h_oh_right = oh_right + np.array([0.96, 0, 0])
    
    positions.extend([h_oh_left, h_oh_right])
    symbols.extend(['H', 'H'])
    
    # H on CH2 groups
    for ch2_pos in [ch2_left, ch2_right]:
        h1 = ch2_pos + np.array([0, 1.1, 0])
        h2 = ch2_pos + np.array([0, -1.1, 0])
        positions.extend([h1, h2])
        symbols.extend(['H', 'H'])
    
    diol = Atoms(symbols=symbols, positions=positions)
    diol.center()
    return diol

def create_hydroxylated_silicon_surface():
    """Create hydroxylated Si(100) surface"""
    # Create Si(100) surface
    si_surf = surface('Si', (1, 0, 0), 4)  # 4 layers
    si_surf = si_surf.repeat((4, 4, 1))    # 4x4 supercell
    
    # Add vacuum
    si_surf.center(vacuum=12.0, axis=2)
    
    # Find top layer Si atoms
    positions = si_surf.get_positions()
    z_max = positions[:, 2].max()
    top_si_indices = []
    
    for i, (pos, symbol) in enumerate(zip(positions, si_surf.get_chemical_symbols())):
        if symbol == 'Si' and abs(pos[2] - z_max) < 0.5:
            top_si_indices.append(i)
    
    # Add OH groups to top Si atoms
    new_atoms = []
    for si_idx in top_si_indices:
        si_pos = positions[si_idx]
        
        # Add oxygen 1.6 Å above Si
        o_pos = si_pos + [0, 0, 1.6]
        new_atoms.append(Atoms('O', positions=[o_pos]))
        
        # Add hydrogen 0.96 Å above oxygen
        h_pos = o_pos + [0, 0, 0.96]
        new_atoms.append(Atoms('H', positions=[h_pos]))
    
    # Combine surface with OH groups
    for atoms in new_atoms:
        si_surf.extend(atoms)
    
    return si_surf

def create_test_structures():
    """Create all initial structures for testing"""
    
    print("Creating initial structures...")
    
    # Create molecules
    print("  Creating TMA molecule...")
    tma = create_tma_molecule()
    write('tma_molecule.xyz', tma)
    print(f"    TMA: {len(tma)} atoms")
    
    print("  Creating 2-butyne-1,4-diol...")
    diol = create_butyne_diol()
    write('butyne_diol_molecule.xyz', diol)
    print(f"    Diol: {len(diol)} atoms")
    
    print("  Creating hydroxylated Si surface...")
    surface_oh = create_hydroxylated_silicon_surface()
    write('hydroxylated_si_surface.xyz', surface_oh)
    print(f"    Surface: {len(surface_oh)} atoms")
    
    # Create combined structure for initial MLD setup
    print("  Creating initial MLD system...")
    
    # Position TMA above surface
    mld_system = surface_oh.copy()
    tma_positioned = tma.copy()
    
    # Find surface center and height
    surf_positions = surface_oh.get_positions()
    surf_center = [surf_positions[:, 0].mean(), surf_positions[:, 1].mean()]
    surf_top = surf_positions[:, 2].max()
    
    # Position TMA 4 Å above surface center
    tma_positioned.translate([surf_center[0], surf_center[1], surf_top + 4.0] - 
                           tma_positioned.get_center_of_mass())
    
    # Combine
    mld_system.extend(tma_positioned)
    write('initial_mld_system.xyz', mld_system)
    print(f"    Initial MLD system: {len(mld_system)} atoms")
    
    print("\nStructures created successfully!")
    print("Files generated:")
    print("  - tma_molecule.xyz")
    print("  - butyne_diol_molecule.xyz") 
    print("  - hydroxylated_si_surface.xyz")
    print("  - initial_mld_system.xyz")
    
    return {
        'tma': tma,
        'diol': diol,
        'surface': surface_oh,
        'mld_system': mld_system
    }

if __name__ == "__main__":
    structures = create_test_structures()