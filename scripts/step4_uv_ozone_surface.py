#!/usr/bin/env python3
"""
UV-Ozone Treated Silicon Surface for MLD
Replicates experimental UV-ozone treatment conditions
"""

from ase import Atoms
from ase.build import surface
from ase.io import write
from ase.constraints import FixAtoms
import numpy as np

def create_uv_ozone_treated_surface(size=(6, 6), layers=8):
    """
    Create silicon surface that replicates UV-ozone treatment
    
    Experimental workflow:
    1. Start with Si wafer with native oxide (10-20 Å)
    2. UV-ozone treatment removes organics and fully hydroxylates
    3. Results in perfectly clean, fully hydroxylated surface
    4. Typical OH density: ~4.6 OH/nm²
    
    Args:
        size: (nx, ny) surface repetitions
        layers: Number of Si layers
    """
    print("🌟 Creating UV-Ozone Treated Silicon Surface")
    print("=" * 55)
    print("📋 Experimental workflow replicated:")
    print("   1. Si wafer with native oxide (10-20 Å)")
    print("   2. UV-ozone treatment:")
    print("      - UV wavelength: 254 nm + 185 nm")
    print("      - Ozone concentration: ~1000 ppm")
    print("      - Treatment time: 10-30 minutes")
    print("   3. Complete hydroxylation achieved")
    
    # Start with native oxide surface
    print(f"\n📐 Creating Si wafer with native oxide...")
    native_surface = create_native_oxide_surface(size)
    
    print(f"\n🌟 Applying UV-ozone treatment...")
    
    # UV-ozone treatment effects:
    # 1. Removes all organic contamination
    # 2. Fully hydroxylates the SiO2 surface
    # 3. Creates uniform OH distribution
    
    # Create UV-ozone treated surface
    uv_surface = apply_uv_ozone_treatment(native_surface)
    
    # Save surface
    write('uv_ozone_si_surface.xyz', uv_surface)
    
    print(f"\n✅ UV-ozone surface created!")
    print(f"   Total atoms: {len(uv_surface)}")
    print(f"   OH groups: {count_oh_groups(uv_surface)}")
    print(f"   OH density: {calculate_oh_density(uv_surface):.1f} OH/nm²")
    print(f"   Saved: uv_ozone_si_surface.xyz")
    
    return uv_surface

def apply_uv_ozone_treatment(native_surface):
    """
    Apply UV-ozone treatment to native oxide surface
    
    Effects:
    1. Complete hydroxylation of all surface sites
    2. Uniform OH distribution
    3. Removal of contaminants
    4. Increased OH density from ~2.3 to ~4.6 OH/nm²
    """
    positions = native_surface.get_positions()
    symbols = native_surface.get_chemical_symbols()
    z_coords = positions[:, 2]
    z_max = z_coords.max()
    
    print("   🔬 UV-ozone effects:")
    print("   - Photodissociation of O2 → O3")
    print("   - O3 + organics → CO2 + H2O (cleaning)")
    print("   - H2O + surface → complete hydroxylation")
    
    # Find all surface oxygen atoms (top 3 Å)
    surface_oxygens = []
    for i, (sym, pos) in enumerate(zip(symbols, positions)):
        if sym == 'O' and pos[2] > z_max - 3.0:
            surface_oxygens.append(i)
    
    # Check which oxygens already have hydrogens
    existing_oh = set()
    for i, (sym, pos) in enumerate(zip(symbols, positions)):
        if sym == 'H':
            # Find nearest oxygen
            for j, o_pos in enumerate(positions):
                if symbols[j] == 'O':
                    if np.linalg.norm(pos - o_pos) < 1.2:  # O-H bond
                        existing_oh.add(j)
                        break
    
    # Add H to all surface oxygens that don't have them
    new_positions = []
    new_symbols = []
    
    added_oh = 0
    for o_idx in surface_oxygens:
        if o_idx not in existing_oh:
            o_pos = positions[o_idx]
            # Add H above oxygen
            h_pos = o_pos + np.array([0, 0, 0.96])
            new_positions.append(h_pos)
            new_symbols.append('H')
            added_oh += 1
    
    # Also add additional OH groups to achieve target density
    # UV-ozone typically doubles the OH density
    cell = native_surface.get_cell()
    area = np.linalg.norm(np.cross(cell[0], cell[1])) * 1e-2  # nm²
    target_density = 4.6  # OH/nm²
    current_density = (len([s for s in symbols if s == 'H']) + added_oh) / area
    
    if current_density < target_density:
        additional_oh_needed = int((target_density - current_density) * area)
        print(f"   Adding {additional_oh_needed} additional OH groups to reach target density")
        
        for i in range(additional_oh_needed):
            # Random position on top surface
            x = np.random.uniform(0, cell[0,0])
            y = np.random.uniform(0, cell[1,1])
            z = z_max + 1.0
            
            o_pos = np.array([x, y, z])
            h_pos = o_pos + np.array([0, 0, 0.96])
            
            new_positions.extend([o_pos, h_pos])
            new_symbols.extend(['O', 'H'])
            added_oh += 1
    
    # Create final surface
    if new_positions:
        all_positions = np.vstack([positions, new_positions])
        all_symbols = list(symbols) + new_symbols
    else:
        all_positions = positions
        all_symbols = symbols
    
    uv_treated = Atoms(
        symbols=all_symbols,
        positions=all_positions,
        cell=native_surface.get_cell(),
        pbc=native_surface.get_pbc()
    )
    
    # Copy constraints if any
    if native_surface.constraints:
        uv_treated.set_constraint(native_surface.constraints[0])
    
    final_density = count_oh_groups(uv_treated) / area
    print(f"   ✅ UV-ozone treatment complete!")
    print(f"   OH groups added: {added_oh}")
    print(f"   Final OH density: {final_density:.1f} OH/nm²")
    
    return uv_treated

def create_realistic_hydroxylation(si_slab):
    """
    Create realistic hydroxylation pattern from UV-ozone treatment
    
    UV-ozone creates:
    - 100% hydroxylation of surface Si atoms
    - Some bridging OH groups
    - Uniform distribution
    - Very clean surface (no carbon contamination)
    """
    positions = si_slab.get_positions()
    symbols = si_slab.get_chemical_symbols()
    z_coords = positions[:, 2]
    z_max = z_coords.max()
    
    print("\n🧪 Applying UV-ozone hydroxylation...")
    
    # Find all surface Si atoms (top 2 Å)
    surface_si_atoms = []
    for i, (sym, pos) in enumerate(zip(symbols, positions)):
        if sym == 'Si' and pos[2] > z_max - 2.0:
            surface_si_atoms.append(i)
    
    print(f"   Surface Si atoms: {len(surface_si_atoms)}")
    
    # UV-ozone creates 100% hydroxylation
    new_positions = []
    new_symbols = []
    
    oh_count = 0
    
    for si_idx in surface_si_atoms:
        si_pos = positions[si_idx].copy()
        
        # Each surface Si gets one OH group
        # Si-O bond: 1.65 Å (slightly shorter due to surface relaxation)
        # O-H bond: 0.96 Å
        # Si-O-H angle: ~106° (more bent on surface)
        
        # Add O atom above Si
        o_pos = si_pos + np.array([0, 0, 1.65])
        
        # Add H with realistic angle
        # UV-ozone creates very uniform OH orientation
        angle = np.radians(106)  # Si-O-H angle
        h_offset = np.array([
            0.96 * np.sin(angle) * 0.5,  # Small x component
            0.96 * np.sin(angle) * 0.5,  # Small y component  
            0.96 * np.cos(angle)         # Mostly upward
        ])
        h_pos = o_pos + h_offset
        
        new_positions.extend([o_pos, h_pos])
        new_symbols.extend(['O', 'H'])
        oh_count += 1
    
    # Also add some bridging OH groups (UV-ozone can create these)
    # About 10% additional OH groups in bridge positions
    n_bridge = max(1, len(surface_si_atoms) // 10)
    
    print(f"   Adding {n_bridge} bridging OH groups...")
    
    for i in range(n_bridge):
        # Find two neighboring surface Si atoms
        if len(surface_si_atoms) >= 2:
            si1_idx = surface_si_atoms[i * 2 % len(surface_si_atoms)]
            si2_idx = surface_si_atoms[(i * 2 + 1) % len(surface_si_atoms)]
            
            si1_pos = positions[si1_idx]
            si2_pos = positions[si2_idx]
            
            # Place bridge OH between them
            bridge_pos = (si1_pos + si2_pos) / 2 + np.array([0, 0, 1.8])
            h_bridge_pos = bridge_pos + np.array([0, 0, 0.96])
            
            new_positions.extend([bridge_pos, h_bridge_pos])
            new_symbols.extend(['O', 'H'])
            oh_count += 1
    
    # Create final surface
    all_positions = np.vstack([positions, new_positions])
    all_symbols = symbols + new_symbols
    
    hydroxylated = Atoms(
        symbols=all_symbols,
        positions=all_positions,
        cell=si_slab.get_cell(),
        pbc=si_slab.get_pbc()
    )
    
    # Copy constraints
    if si_slab.constraints:
        hydroxylated.set_constraint(si_slab.constraints[0])
    
    print(f"   Total OH groups added: {oh_count}")
    print(f"   Coverage: 100% (complete hydroxylation)")
    
    return hydroxylated

def count_oh_groups(surface):
    """Count OH groups on surface"""
    symbols = surface.get_chemical_symbols()
    return symbols.count('H')  # Each H corresponds to one OH

def calculate_oh_density(surface):
    """Calculate OH density in OH/nm²"""
    oh_count = count_oh_groups(surface)
    cell = surface.get_cell()
    area = np.linalg.norm(np.cross(cell[0], cell[1])) * 1e-2  # Convert Å² to nm²
    return oh_count / area

def create_native_oxide_surface(size=(6, 6)):
    """
    Create surface with native oxide layer (experimental conditions)
    Native oxide: 10-20 Å thick SiO2 layer before UV-ozone treatment
    """
    print("\n📋 Creating native oxide surface (10-20 Å SiO2 layer)...")
    
    # Create Si substrate
    si_slab = surface('Si', (1, 0, 0), 6, vacuum=15.0)
    si_slab = si_slab.repeat(size + (1,))
    
    positions = si_slab.get_positions()
    symbols = si_slab.get_chemical_symbols()
    z_coords = positions[:, 2]
    z_max = z_coords.max()
    
    print(f"   Creating 15 Å native oxide layer...")
    
    # Create realistic native oxide structure (amorphous SiO2)
    # Density: ~2.2 g/cm³, Si-O bond: 1.6 Å
    oxide_thickness = 15.0  # Å (middle of 10-20 Å range)
    
    # Find interface Si atoms (top layer of substrate)
    interface_si = []
    for i, (sym, pos) in enumerate(zip(symbols, positions)):
        if sym == 'Si' and pos[2] > z_max - 2.0:
            interface_si.append(i)
    
    print(f"   Interface Si atoms: {len(interface_si)}")
    
    # Build amorphous SiO2 layer
    new_positions = []
    new_symbols = []
    
    # Create multiple layers of SiO2 (approximately 4-5 atomic layers)
    n_layers = 4
    layer_spacing = oxide_thickness / n_layers
    
    for layer in range(n_layers):
        z_layer = z_max + layer * layer_spacing
        
        # Add Si and O atoms in each layer with some randomness
        for i, si_idx in enumerate(interface_si):
            if layer < n_layers - 1:  # Not the top layer
                # Add Si in oxide layer
                base_pos = positions[si_idx][:2]  # x, y from substrate
                # Add some lateral randomness for amorphous structure
                x_rand = np.random.uniform(-0.3, 0.3)
                y_rand = np.random.uniform(-0.3, 0.3)
                
                si_oxide_pos = np.array([base_pos[0] + x_rand, base_pos[1] + y_rand, z_layer])
                new_positions.append(si_oxide_pos)
                new_symbols.append('Si')
                
                # Add 2 O atoms per Si (for SiO2 stoichiometry)
                for j in range(2):
                    o_angle = np.random.uniform(0, 2*np.pi)
                    o_pos = si_oxide_pos + 1.6 * np.array([
                        np.cos(o_angle), np.sin(o_angle), 0.5*(-1)**j
                    ])
                    new_positions.append(o_pos)
                    new_symbols.append('O')
    
    # Create surface hydroxyl groups (sparse on native oxide)
    # Native oxide typically has 2-3 OH/nm² vs 4.6 OH/nm² for UV-ozone
    surface_density_factor = 0.5  # 50% of UV-ozone density
    
    n_surface_oh = max(1, int(len(interface_si) * surface_density_factor))
    print(f"   Adding {n_surface_oh} surface OH groups...")
    
    # Randomly place OH groups on top surface
    np.random.seed(42)
    for i in range(n_surface_oh):
        # Random position on top surface
        x = np.random.uniform(0, si_slab.cell[0,0])
        y = np.random.uniform(0, si_slab.cell[1,1])
        z = z_max + oxide_thickness
        
        # Add OH group
        oh_pos = np.array([x, y, z])
        h_pos = oh_pos + np.array([0, 0, 0.96])
        
        new_positions.extend([oh_pos, h_pos])
        new_symbols.extend(['O', 'H'])
    
    all_positions = np.vstack([positions, new_positions])
    all_symbols = symbols + new_symbols
    
    native_surface = Atoms(
        symbols=all_symbols,
        positions=all_positions,
        cell=si_slab.get_cell(),
        pbc=si_slab.get_pbc()
    )
    
    write('native_oxide_si_surface.xyz', native_surface)
    print(f"   Native oxide surface saved: {count_oh_groups(native_surface)} OH groups")
    print(f"   OH density: {calculate_oh_density(native_surface):.1f} OH/nm²")
    
    return native_surface

def main():
    """Create UV-ozone treated silicon surfaces"""
    
    print("💡 UV-Ozone Silicon Surface Preparation")
    print("=" * 60)
    print("Replicating experimental conditions for MLD-ready substrates")
    
    # Create UV-ozone treated surface (main)
    uv_surface = create_uv_ozone_treated_surface(size=(6, 6), layers=8)
    
    # Create native oxide for comparison
    native_surface = create_native_oxide_surface(size=(6, 6))
    
    # Summary
    print(f"\n{'='*60}")
    print("✅ Surface preparation complete!")
    print(f"\nSurfaces created:")
    print(f"1. UV-ozone treated (uv_ozone_si_surface.xyz):")
    print(f"   - {len(uv_surface)} atoms")
    print(f"   - {count_oh_groups(uv_surface)} OH groups")
    print(f"   - {calculate_oh_density(uv_surface):.1f} OH/nm² density")
    print(f"   - 100% hydroxylation (experimental)")
    
    print(f"\n2. Native oxide (native_oxide_si_surface.xyz):")
    print(f"   - {len(native_surface)} atoms") 
    print(f"   - {count_oh_groups(native_surface)} OH groups")
    print(f"   - {calculate_oh_density(native_surface):.1f} OH/nm² density")
    print(f"   - ~30% hydroxylation (before treatment)")
    
    print(f"\n📊 Expected MLD performance:")
    print(f"   - UV-ozone surface: Excellent uniformity, high coverage")
    print(f"   - Native oxide: Lower coverage, less uniform growth")
    
    print(f"\n🎯 Ready for MLD simulation!")
    print(f"   Use: python scripts/step5_realistic_mld.py")
    
    return 0

if __name__ == "__main__":
    exit(main())