#!/usr/bin/env python3
"""
Small Surface Generator for DFT MLD Studies
Creates literature-validated small surface models suitable for DFT calculations

Based on best practices from:
- Halls & Raghavachari (2004): Si9O12H12 clusters
- Mameli et al. (2017): 2×2 Si(100) slabs
- Shirazi & Elliott (2013): Small representative models
"""

from ase import Atoms
from ase.build import surface, bulk
from ase.io import write
from ase.constraints import FixAtoms
import numpy as np
from typing import Tuple, List

class SmallSurfaceGenerator:
    """Generate DFT-appropriate small surface models"""
    
    def __init__(self):
        """Initialize generator with literature-based parameters"""
        
        # Size constraints based on literature
        self.MAX_ATOMS = 200        # Hard computational limit
        self.WARN_ATOMS = 150       # Warning threshold
        self.TARGET_ATOMS = 80      # Typical DFT MLD size
        
        # Surface parameters
        self.SI_LATTICE = 5.431     # Å, experimental Si lattice
        self.VACUUM = 15.0          # Å, minimum for DFT
        
        print("🔬 Small Surface Generator for DFT MLD")
        print(f"   Target size: ~{self.TARGET_ATOMS} atoms")
        print(f"   Maximum size: {self.MAX_ATOMS} atoms")
        print(f"   Based on literature best practices")
    
    def create_si100_2x2_slab(self, layers: int = 4) -> Atoms:
        """
        Create Si(100) 2×2 slab - Literature standard
        
        References:
        - Mameli et al. (2017): Used 2×2 and 3×3 slabs
        - Standard in ALD/MLD DFT studies
        
        Args:
            layers: Number of Si layers (4-6 typical)
            
        Returns:
            ASE Atoms object of clean Si(100) slab
        """
        print(f"\n📐 Creating Si(100) 2×2 slab ({layers} layers)")
        
        # Create Si(100) surface using ASE
        si_slab = surface('Si', (1, 0, 0), layers, vacuum=self.VACUUM)
        
        # Make 2×2 supercell (literature standard size)
        si_slab = si_slab.repeat((2, 2, 1))
        
        # Fix bottom layers (literature standard: fix bottom half)
        positions = si_slab.get_positions()
        z_coords = positions[:, 2]
        z_min = z_coords.min()
        
        # Fix bottom 2 layers
        fix_height = z_min + 5.0  # Approximately 2 layers
        fix_mask = [pos[2] < fix_height for pos in positions]
        n_fixed = sum(fix_mask)
        
        si_slab.set_constraint(FixAtoms(mask=fix_mask))
        
        print(f"   Atoms: {len(si_slab)}")
        print(f"   Cell: {si_slab.cell[0,0]:.1f} × {si_slab.cell[1,1]:.1f} Å²")
        print(f"   Fixed atoms: {n_fixed} (bottom layers)")
        
        if len(si_slab) > self.WARN_ATOMS:
            print(f"   ⚠️  Warning: {len(si_slab)} atoms > {self.WARN_ATOMS}")
        
        return si_slab
    
    def create_si100_3x3_slab(self, layers: int = 4) -> Atoms:
        """
        Create Si(100) 3×3 slab - For coverage studies
        
        Larger than typical, use only if needed for coverage effects
        """
        print(f"\n📐 Creating Si(100) 3×3 slab ({layers} layers)")
        
        si_slab = surface('Si', (1, 0, 0), layers, vacuum=self.VACUUM)
        si_slab = si_slab.repeat((3, 3, 1))
        
        # Fix bottom layers
        positions = si_slab.get_positions()
        z_coords = positions[:, 2]
        z_min = z_coords.min()
        
        fix_height = z_min + 5.0
        fix_mask = [pos[2] < fix_height for pos in positions]
        n_fixed = sum(fix_mask)
        
        si_slab.set_constraint(FixAtoms(mask=fix_mask))
        
        print(f"   Atoms: {len(si_slab)}")
        print(f"   Cell: {si_slab.cell[0,0]:.1f} × {si_slab.cell[1,1]:.1f} Å²")
        print(f"   Fixed atoms: {n_fixed}")
        
        if len(si_slab) > self.MAX_ATOMS:
            print(f"   ❌ ERROR: {len(si_slab)} atoms > {self.MAX_ATOMS} limit!")
            print(f"   This system is too large for practical DFT")
            return None
        elif len(si_slab) > self.WARN_ATOMS:
            print(f"   ⚠️  Warning: Large system ({len(si_slab)} atoms)")
        
        return si_slab
    
    def hydroxylate_surface(self, slab: Atoms, coverage: float = 1.0) -> Atoms:
        """
        Add OH groups to surface Si atoms
        
        Args:
            slab: Clean Si surface
            coverage: OH coverage (0.0 to 1.0)
            
        Returns:
            Hydroxylated surface
        """
        print(f"\n🧪 Hydroxylating surface (coverage: {coverage:.0%})")
        
        positions = slab.get_positions()
        symbols = slab.get_chemical_symbols()
        z_coords = positions[:, 2]
        z_max = z_coords.max()
        
        # Find surface Si atoms (top 2 Å)
        surface_si = []
        for i, (sym, pos) in enumerate(zip(symbols, positions)):
            if sym == 'Si' and pos[2] > z_max - 2.0:
                surface_si.append(i)
        
        print(f"   Surface Si atoms: {len(surface_si)}")
        
        # Select sites for hydroxylation
        n_oh = int(len(surface_si) * coverage)
        if coverage < 1.0:
            np.random.seed(42)  # Reproducible
            oh_sites = np.random.choice(surface_si, size=n_oh, replace=False)
        else:
            oh_sites = surface_si
        
        # Add OH groups
        new_positions = []
        new_symbols = []
        
        for si_idx in oh_sites:
            si_pos = positions[si_idx]
            
            # Add O above Si (Si-O bond: 1.65 Å)
            o_pos = si_pos + np.array([0, 0, 1.65])
            
            # Add H above O (O-H bond: 0.96 Å)
            # Random orientation for realistic surface
            theta = np.random.uniform(0, 2*np.pi)
            phi = np.radians(20)  # Slight tilt from vertical
            
            h_offset = 0.96 * np.array([
                np.sin(phi) * np.cos(theta),
                np.sin(phi) * np.sin(theta),
                np.cos(phi)
            ])
            h_pos = o_pos + h_offset
            
            new_positions.extend([o_pos, h_pos])
            new_symbols.extend(['O', 'H'])
        
        # Create hydroxylated surface
        all_positions = np.vstack([positions, new_positions])
        all_symbols = symbols + new_symbols
        
        hydroxylated = Atoms(
            symbols=all_symbols,
            positions=all_positions,
            cell=slab.get_cell(),
            pbc=slab.get_pbc()
        )
        
        # Copy constraints
        if slab.constraints:
            hydroxylated.set_constraint(slab.constraints[0])
        
        print(f"   OH groups added: {n_oh}")
        print(f"   Total atoms: {len(hydroxylated)}")
        
        if len(hydroxylated) > self.MAX_ATOMS:
            print(f"   ❌ ERROR: {len(hydroxylated)} atoms > {self.MAX_ATOMS} limit!")
            return None
        
        return hydroxylated
    
    def create_native_oxide_layer(self, si_slab: Atoms, thickness: float = 1.0) -> Atoms:
        """
        Add thin native oxide layer (simplified)
        
        Args:
            si_slab: Si substrate
            thickness: Oxide thickness in Å
            
        Returns:
            Si with thin SiO2 layer
        """
        print(f"\n🔬 Adding native oxide layer ({thickness:.1f} Å)")
        
        positions = si_slab.get_positions()
        symbols = si_slab.get_chemical_symbols()
        z_coords = positions[:, 2]
        z_max = z_coords.max()
        
        # Find interface Si atoms
        interface_si = []
        for i, (sym, pos) in enumerate(zip(symbols, positions)):
            if sym == 'Si' and pos[2] > z_max - 1.0:
                interface_si.append(i)
        
        # Add thin oxide layer (simplified as SiO2 units)
        new_positions = []
        new_symbols = []
        
        for si_idx in interface_si[::2]:  # Every other Si atom
            si_pos = positions[si_idx]
            
            # Add O above Si
            o1_pos = si_pos + np.array([0.8, 0, thickness])
            o2_pos = si_pos + np.array([-0.8, 0, thickness])
            
            new_positions.extend([o1_pos, o2_pos])
            new_symbols.extend(['O', 'O'])
        
        # Combine
        all_positions = np.vstack([positions, new_positions])
        all_symbols = symbols + new_symbols
        
        oxide_surface = Atoms(
            symbols=all_symbols,
            positions=all_positions,
            cell=si_slab.get_cell(),
            pbc=si_slab.get_pbc()
        )
        
        # Copy constraints
        if si_slab.constraints:
            oxide_surface.set_constraint(si_slab.constraints[0])
        
        print(f"   Oxide oxygens added: {len(new_symbols)}")
        print(f"   Total atoms: {len(oxide_surface)}")
        
        return oxide_surface
    
    def validate_surface_size(self, surface: Atoms, name: str) -> bool:
        """Validate surface is appropriate for DFT"""
        n_atoms = len(surface)
        
        print(f"\n✅ Validating {name}")
        print(f"   Atoms: {n_atoms}")
        
        if n_atoms > self.MAX_ATOMS:
            print(f"   ❌ TOO LARGE: {n_atoms} > {self.MAX_ATOMS}")
            print(f"   This will fail in DFT calculations")
            return False
        elif n_atoms > self.WARN_ATOMS:
            print(f"   ⚠️  Large: {n_atoms} > {self.WARN_ATOMS}")
            print(f"   May be slow but should work")
            return True
        else:
            print(f"   ✅ Good size: {n_atoms} atoms")
            return True

def main():
    """Generate small surfaces for DFT MLD studies"""
    
    generator = SmallSurfaceGenerator()
    
    print("\n" + "="*60)
    print("CREATING SMALL SURFACES FOR DFT MLD")
    print("="*60)
    
    surfaces_created = []
    
    # 1. Si(100) 2×2 slab - 4 layers (RECOMMENDED)
    print(f"\n{'='*50}")
    print("1. Si(100) 2×2 Slab (4 layers) - RECOMMENDED")
    print("="*50)
    
    si_2x2_4L = generator.create_si100_2x2_slab(layers=4)
    if generator.validate_surface_size(si_2x2_4L, "Si(100) 2×2 4L"):
        write('small_si100_2x2_4L.xyz', si_2x2_4L)
        surfaces_created.append(('small_si100_2x2_4L.xyz', len(si_2x2_4L)))
        print(f"   💾 Saved: small_si100_2x2_4L.xyz")
    
    # 2. Hydroxylated version
    si_2x2_4L_OH = generator.hydroxylate_surface(si_2x2_4L, coverage=1.0)
    if si_2x2_4L_OH and generator.validate_surface_size(si_2x2_4L_OH, "Hydroxylated 2×2"):
        write('small_si100_2x2_4L_OH.xyz', si_2x2_4L_OH)
        surfaces_created.append(('small_si100_2x2_4L_OH.xyz', len(si_2x2_4L_OH)))
        print(f"   💾 Saved: small_si100_2x2_4L_OH.xyz")
    
    # 3. Si(100) 2×2 slab - 6 layers (Alternative)
    print(f"\n{'='*50}")
    print("2. Si(100) 2×2 Slab (6 layers) - Alternative")
    print("="*50)
    
    si_2x2_6L = generator.create_si100_2x2_slab(layers=6)
    if generator.validate_surface_size(si_2x2_6L, "Si(100) 2×2 6L"):
        write('small_si100_2x2_6L.xyz', si_2x2_6L)
        surfaces_created.append(('small_si100_2x2_6L.xyz', len(si_2x2_6L)))
        print(f"   💾 Saved: small_si100_2x2_6L.xyz")
    
    # 4. Hydroxylated 6-layer version
    si_2x2_6L_OH = generator.hydroxylate_surface(si_2x2_6L, coverage=1.0)
    if si_2x2_6L_OH and generator.validate_surface_size(si_2x2_6L_OH, "Hydroxylated 2×2 6L"):
        write('small_si100_2x2_6L_OH.xyz', si_2x2_6L_OH)
        surfaces_created.append(('small_si100_2x2_6L_OH.xyz', len(si_2x2_6L_OH)))
        print(f"   💾 Saved: small_si100_2x2_6L_OH.xyz")
    
    # 5. Si(100) 3×3 slab - Only if requested (may be too large)
    print(f"\n{'='*50}")
    print("3. Si(100) 3×3 Slab (4 layers) - For coverage studies")
    print("="*50)
    
    si_3x3_4L = generator.create_si100_3x3_slab(layers=4)
    if si_3x3_4L and generator.validate_surface_size(si_3x3_4L, "Si(100) 3×3 4L"):
        write('small_si100_3x3_4L.xyz', si_3x3_4L)
        surfaces_created.append(('small_si100_3x3_4L.xyz', len(si_3x3_4L)))
        print(f"   💾 Saved: small_si100_3x3_4L.xyz")
        
        # Hydroxylated 3×3 (check size carefully)
        si_3x3_4L_OH = generator.hydroxylate_surface(si_3x3_4L, coverage=1.0)
        if si_3x3_4L_OH and generator.validate_surface_size(si_3x3_4L_OH, "Hydroxylated 3×3"):
            write('small_si100_3x3_4L_OH.xyz', si_3x3_4L_OH)
            surfaces_created.append(('small_si100_3x3_4L_OH.xyz', len(si_3x3_4L_OH)))
            print(f"   💾 Saved: small_si100_3x3_4L_OH.xyz")
    
    # Summary
    print(f"\n{'='*60}")
    print("✅ SMALL SURFACE GENERATION COMPLETE")
    print("="*60)
    
    print(f"\nSurfaces created: {len(surfaces_created)}")
    for filename, n_atoms in surfaces_created:
        status = "✅ DFT-ready" if n_atoms <= generator.WARN_ATOMS else "⚠️  Large"
        print(f"  📄 {filename}: {n_atoms} atoms ({status})")
    
    print(f"\n🎯 Recommended for DFT MLD:")
    print(f"  • small_si100_2x2_4L_OH.xyz (best balance)")
    print(f"  • small_si100_2x2_6L_OH.xyz (more bulk)")
    
    print(f"\n💡 Next steps:")
    print(f"  1. Use geometry_validator.py to check structures")
    print(f"  2. Run dft_small_mld.py with these surfaces")
    print(f"  3. Start with smallest system first")
    
    return 0

if __name__ == "__main__":
    exit(main())