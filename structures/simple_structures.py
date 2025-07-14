#!/usr/bin/env python3
"""
Create simple initial structures without ASE dependencies
Generate XYZ files directly for initial testing
"""

import numpy as np

def write_xyz(filename, symbols, positions, comment=""):
    """Write XYZ file format"""
    with open(filename, 'w') as f:
        f.write(f"{len(symbols)}\n")
        f.write(f"{comment}\n")
        for symbol, pos in zip(symbols, positions):
            f.write(f"{symbol}  {pos[0]:12.6f}  {pos[1]:12.6f}  {pos[2]:12.6f}\n")

def create_tma_xyz():
    """Create TMA molecule XYZ"""
    symbols = ['Al']
    positions = [[0.0, 0.0, 0.0]]
    
    # Three methyl groups
    methyl_positions = [
        [0.0, 0.0, 2.0],
        [1.73, 0.0, -1.0],
        [-1.73, 0.0, -1.0]
    ]
    
    for pos in methyl_positions:
        symbols.append('C')
        positions.append(pos)
        
        # Add hydrogens
        c_pos = np.array(pos)
        for i in range(3):
            angle = i * 2 * np.pi / 3
            h_pos = c_pos + 1.1 * np.array([
                0.5 * np.cos(angle), 
                0.5 * np.sin(angle), 
                0.3
            ])
            symbols.append('H')
            positions.append(h_pos.tolist())
    
    write_xyz('tma_molecule.xyz', symbols, positions, 
              "Trimethylaluminum (TMA) - Al(CH3)3")
    return symbols, positions

def create_diol_xyz():
    """Create 2-butyne-1,4-diol XYZ"""
    symbols = []
    positions = []
    
    # Acetylene backbone C≡C
    symbols.extend(['C', 'C'])
    positions.extend([[-0.6, 0.0, 0.0], [0.6, 0.0, 0.0]])
    
    # CH2 carbons
    symbols.extend(['C', 'C'])
    positions.extend([[-2.1, 0.0, 0.0], [2.1, 0.0, 0.0]])
    
    # Oxygen atoms
    symbols.extend(['O', 'O'])
    positions.extend([[-3.5, 0.0, 0.0], [3.5, 0.0, 0.0]])
    
    # OH hydrogens
    symbols.extend(['H', 'H'])
    positions.extend([[-4.46, 0.0, 0.0], [4.46, 0.0, 0.0]])
    
    # CH2 hydrogens
    ch2_h_positions = [
        [-2.1, 1.1, 0.0], [-2.1, -1.1, 0.0],  # Left CH2
        [2.1, 1.1, 0.0], [2.1, -1.1, 0.0]     # Right CH2
    ]
    symbols.extend(['H'] * 4)
    positions.extend(ch2_h_positions)
    
    write_xyz('butyne_diol_molecule.xyz', symbols, positions,
              "2-butyne-1,4-diol - HOCH2-C≡C-CH2OH")
    return symbols, positions

def create_simple_surface():
    """Create simple Si surface with OH groups"""
    symbols = []
    positions = []
    
    # 3x3 Si surface grid
    a = 5.43  # Si lattice parameter
    
    # Bottom Si layer
    for i in range(3):
        for j in range(3):
            x = i * a
            y = j * a
            z = 0.0
            symbols.append('Si')
            positions.append([x, y, z])
    
    # Top Si layer (offset)
    for i in range(3):
        for j in range(3):
            x = i * a + a/2
            y = j * a + a/2
            z = a/4
            symbols.append('Si')
            positions.append([x, y, z])
    
    # Add OH groups on top layer
    for i in range(3):
        for j in range(3):
            x = i * a + a/2
            y = j * a + a/2
            z = a/4
            
            # Oxygen above Si
            symbols.append('O')
            positions.append([x, y, z + 1.6])
            
            # Hydrogen above oxygen
            symbols.append('H')
            positions.append([x, y, z + 2.56])
    
    write_xyz('simple_si_surface.xyz', symbols, positions,
              "Simple hydroxylated Si surface")
    return symbols, positions

def create_initial_mld_system():
    """Create initial MLD system with TMA above surface"""
    # Start with surface
    surf_symbols, surf_positions = create_simple_surface()
    
    # Add TMA above surface center
    tma_symbols, tma_positions = create_tma_xyz()
    
    # Offset TMA above surface
    surf_z_max = max(pos[2] for pos in surf_positions)
    tma_offset = [7.0, 7.0, surf_z_max + 4.0]  # Center above surface
    
    # Combine structures
    combined_symbols = surf_symbols + tma_symbols
    combined_positions = surf_positions + [
        [pos[0] + tma_offset[0], pos[1] + tma_offset[1], pos[2] + tma_offset[2]]
        for pos in tma_positions
    ]
    
    write_xyz('initial_mld_system.xyz', combined_symbols, combined_positions,
              "Initial MLD system: TMA above hydroxylated Si surface")
    
    return combined_symbols, combined_positions

def main():
    """Create all initial structures"""
    print("Creating initial structures for MLD simulation...")
    
    print("  Creating TMA molecule...")
    create_tma_xyz()
    
    print("  Creating 2-butyne-1,4-diol...")
    create_diol_xyz()
    
    print("  Creating Si surface...")
    create_simple_surface()
    
    print("  Creating initial MLD system...")
    create_initial_mld_system()
    
    print("\nStructures created successfully!")
    print("Files generated:")
    print("  - tma_molecule.xyz")
    print("  - butyne_diol_molecule.xyz")
    print("  - simple_si_surface.xyz")
    print("  - initial_mld_system.xyz")

if __name__ == "__main__":
    main()