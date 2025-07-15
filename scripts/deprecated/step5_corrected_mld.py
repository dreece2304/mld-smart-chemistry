#!/usr/bin/env python3
"""
Corrected MLD Simulation with Proper Chemistry
Fixes the methyl group accumulation and implements proper Al-O-Al growth
"""

from ase import Atoms
from ase.io import read, write
from ase.constraints import FixAtoms
import numpy as np
from typing import List, Tuple, Optional
import time

class CorrectedMLDSimulator:
    """
    MLD simulator with proper chemical reactions:
    
    TMA Pulse: Si-OH + Al(CH3)3 → Si-O-Al(CH3)2 + CH4↑
    H2O Pulse: Si-O-Al(CH3)2 + H2O → Si-O-Al(OH)CH3 + CH4↑
    Next TMA: Si-O-Al(OH) + Al(CH3)3 → Si-O-Al-O-Al(CH3)2 + CH4↑
    
    This creates proper Al-O-Al bridge formation!
    """
    
    def __init__(self, surface: Atoms, temperature: float = 473.15):
        """Initialize corrected MLD simulator"""
        self.surface = surface.copy()
        self.temperature = temperature
        self.cycle_count = 0
        
        # Track reactive sites properly
        self.oh_sites = []  # Surface OH groups
        self.al_sites = []  # Al atoms with attached methyls
        
        self._identify_surface_sites()
        
        print(f"🔧 Corrected MLD Simulator initialized")
        print(f"   Surface: {len(self.surface)} atoms")
        print(f"   Initial OH sites: {len(self.oh_sites)}")
    
    def _identify_surface_sites(self):
        """Find reactive OH groups on surface"""
        positions = self.surface.get_positions()
        symbols = self.surface.get_chemical_symbols()
        
        self.oh_sites = []
        
        # Find O-H pairs on surface (top 3 Å)
        z_coords = positions[:, 2]
        z_max = z_coords.max()
        
        for i, sym_i in enumerate(symbols):
            if sym_i == 'O' and positions[i][2] > z_max - 3.0:
                # Look for H bonded to this O
                for j, sym_j in enumerate(symbols):
                    if sym_j == 'H':
                        distance = np.linalg.norm(positions[i] - positions[j])
                        if 0.85 < distance < 1.15:  # O-H bond
                            self.oh_sites.append({
                                'O_idx': i,
                                'H_idx': j,
                                'position': positions[i],
                                'type': 'surface_OH',  # vs 'Al_OH' later
                                'available': True
                            })
                            break
    
    def _find_al_sites(self):
        """Find Al atoms that can react with H2O"""
        positions = self.surface.get_positions()
        symbols = self.surface.get_chemical_symbols()
        
        al_sites = []
        for i, sym in enumerate(symbols):
            if sym == 'Al':
                # Check if this Al has methyl groups
                al_pos = positions[i]
                bonded_carbons = []
                
                for j, sym_j in enumerate(symbols):
                    if sym_j == 'C':
                        c_pos = positions[j]
                        distance = np.linalg.norm(al_pos - c_pos)
                        if 1.8 < distance < 2.1:  # Al-C bond
                            bonded_carbons.append(j)
                
                if len(bonded_carbons) > 0:
                    al_sites.append({
                        'Al_idx': i,
                        'position': al_pos,
                        'bonded_carbons': bonded_carbons,
                        'methyl_count': len(bonded_carbons)
                    })
        
        return al_sites
    
    def tma_pulse(self) -> int:
        """
        TMA pulse with proper chemistry
        
        Reactions:
        1. Si-OH + Al(CH3)3 → Si-O-Al(CH3)2 + CH4↑
        2. Si-O-Al(OH) + Al(CH3)3 → Si-O-Al-O-Al(CH3)2 + CH4↑ (bridge formation!)
        """
        print(f"\n💨 TMA Pulse (Cycle {self.cycle_count + 1})")
        
        # Find available reaction sites
        available_oh = [site for site in self.oh_sites if site['available']]
        
        print(f"   Available OH sites: {len(available_oh)}")
        
        if not available_oh:
            print("   ⚠️  No available OH sites for TMA reaction")
            return 0
        
        reacted_sites = 0
        
        # React TMA with available OH sites
        for site in available_oh[:3]:  # Limit reactions per pulse (steric hindrance)
            success = self._react_tma_at_oh_site(site)
            if success:
                reacted_sites += 1
                site['available'] = False
        
        print(f"   ✅ TMA molecules reacted: {reacted_sites}")
        return reacted_sites
    
    def _react_tma_at_oh_site(self, oh_site):
        """
        React TMA at specific OH site
        Si-OH + Al(CH3)3 → Si-O-Al(CH3)2 + CH4↑
        """
        positions = self.surface.get_positions()
        symbols = self.surface.get_chemical_symbols()
        
        o_idx = oh_site['O_idx']
        h_idx = oh_site['H_idx']
        o_pos = positions[o_idx]
        
        # Remove the H (represents CH4 elimination)
        new_positions = []
        new_symbols = []
        
        for i, (pos, sym) in enumerate(zip(positions, symbols)):
            if i != h_idx:  # Skip the H being removed
                new_positions.append(pos)
                new_symbols.append(sym)
        
        # Add Al(CH3)2 group
        # Al position: 1.9 Å from O
        al_pos = o_pos + np.array([0, 0, 1.9])
        
        # Add two methyl groups in tetrahedral geometry
        # Al-C bond: 1.96 Å
        c1_angle = np.pi/3  # 60°
        c2_angle = -np.pi/3  # -60°
        
        c1_pos = al_pos + 1.96 * np.array([np.cos(c1_angle), np.sin(c1_angle), 0.3])
        c2_pos = al_pos + 1.96 * np.array([np.cos(c2_angle), np.sin(c2_angle), 0.3])
        
        # Add H atoms to methyls (simplified - 1 H per C for visualization)
        h1_pos = c1_pos + np.array([0, 0, 1.09])
        h2_pos = c2_pos + np.array([0, 0, 1.09])
        h3_pos = c1_pos + np.array([1.09, 0, 0])
        h4_pos = c2_pos + np.array([-1.09, 0, 0])
        
        # Add new atoms
        new_positions.extend([al_pos, c1_pos, c2_pos, h1_pos, h2_pos, h3_pos, h4_pos])
        new_symbols.extend(['Al', 'C', 'C', 'H', 'H', 'H', 'H'])
        
        # Update surface
        self.surface = Atoms(
            symbols=new_symbols,
            positions=new_positions,
            cell=self.surface.get_cell(),
            pbc=self.surface.get_pbc()
        )
        
        # Update constraints if present
        if self.surface.constraints:
            # Simplified: just keep the constraint structure
            pass
        
        return True
    
    def h2o_pulse(self) -> int:
        """
        H2O pulse with proper chemistry
        
        Reactions:
        1. Si-O-Al(CH3)2 + H2O → Si-O-Al(OH)CH3 + CH4↑
        2. Si-O-Al(CH3)(OH) + H2O → Si-O-Al(OH)2 + CH4↑
        """
        print(f"\n💧 H2O Pulse (Cycle {self.cycle_count + 1})")
        
        # Find Al sites that can react with H2O
        al_sites = self._find_al_sites()
        
        print(f"   Al sites with methyls: {len(al_sites)}")
        
        if not al_sites:
            print("   ⚠️  No Al-CH3 sites for H2O reaction")
            return 0
        
        reacted_sites = 0
        
        # React H2O with Al-CH3 groups
        for al_site in al_sites:
            if al_site['methyl_count'] > 0:
                success = self._react_h2o_at_al_site(al_site)
                if success:
                    reacted_sites += 1
        
        # Create new OH sites for next cycle
        self._create_new_oh_sites(reacted_sites)
        
        print(f"   ✅ H2O reactions: {reacted_sites}")
        print(f"   New OH sites created: {reacted_sites}")
        
        return reacted_sites
    
    def _react_h2o_at_al_site(self, al_site):
        """
        React H2O with Al-CH3 group
        Al-CH3 + H2O → Al-OH + CH4↑
        """
        positions = self.surface.get_positions()
        symbols = self.surface.get_chemical_symbols()
        
        al_idx = al_site['Al_idx']
        al_pos = positions[al_idx]
        
        # Remove one methyl group (CH4 elimination)
        if al_site['bonded_carbons']:
            c_to_remove = al_site['bonded_carbons'][0]  # Remove first methyl
            
            # Find H atoms bonded to this C
            h_to_remove = []
            c_pos = positions[c_to_remove]
            
            for i, (pos, sym) in enumerate(zip(positions, symbols)):
                if sym == 'H':
                    distance = np.linalg.norm(pos - c_pos)
                    if distance < 1.2:  # C-H bond
                        h_to_remove.append(i)
            
            # Remove C and its H atoms
            indices_to_remove = [c_to_remove] + h_to_remove
            
            new_positions = []
            new_symbols = []
            
            for i, (pos, sym) in enumerate(zip(positions, symbols)):
                if i not in indices_to_remove:
                    new_positions.append(pos)
                    new_symbols.append(sym)
            
            # Add OH group to Al
            oh_pos = al_pos + np.array([1.8, 0, 0])  # Al-O bond
            h_pos = oh_pos + np.array([0, 0, 0.96])  # O-H bond
            
            new_positions.extend([oh_pos, h_pos])
            new_symbols.extend(['O', 'H'])
            
            # Update surface
            self.surface = Atoms(
                symbols=new_symbols,
                positions=new_positions,
                cell=self.surface.get_cell(),
                pbc=self.surface.get_pbc()
            )
            
            return True
        
        return False
    
    def _create_new_oh_sites(self, count):
        """Create new OH sites for next TMA reaction"""
        # Update OH sites list to include new Al-OH groups
        positions = self.surface.get_positions()
        symbols = self.surface.get_chemical_symbols()
        
        # Find new O-H pairs (Al-OH groups)
        new_oh_sites = []
        
        for i, sym_i in enumerate(symbols):
            if sym_i == 'O':
                # Check if this O is bonded to Al (not Si)
                o_pos = positions[i]
                bonded_to_al = False
                
                for j, sym_j in enumerate(symbols):
                    if sym_j == 'Al':
                        al_pos = positions[j]
                        if np.linalg.norm(o_pos - al_pos) < 2.0:  # Al-O bond
                            bonded_to_al = True
                            break
                
                if bonded_to_al:
                    # Look for H bonded to this O
                    for k, sym_k in enumerate(symbols):
                        if sym_k == 'H':
                            h_pos = positions[k]
                            if np.linalg.norm(o_pos - h_pos) < 1.2:  # O-H bond
                                new_oh_sites.append({
                                    'O_idx': i,
                                    'H_idx': k,
                                    'position': o_pos,
                                    'type': 'Al_OH',
                                    'available': True
                                })
                                break
        
        # Add new OH sites
        self.oh_sites.extend(new_oh_sites)
    
    def run_cycle(self):
        """Run one complete corrected MLD cycle"""
        print(f"\n{'='*60}")
        print(f"🔄 Corrected MLD Cycle {self.cycle_count + 1}")
        print(f"{'='*60}")
        
        # TMA pulse
        tma_reactions = self.tma_pulse()
        
        # Save after TMA
        write(f'corrected_cycle{self.cycle_count+1}_after_tma.xyz', self.surface)
        
        # H2O pulse
        h2o_reactions = self.h2o_pulse()
        
        # Save after H2O
        write(f'corrected_cycle{self.cycle_count+1}_after_h2o.xyz', self.surface)
        
        self.cycle_count += 1
        
        # Analysis
        symbols = self.surface.get_chemical_symbols()
        composition = {
            'Al': symbols.count('Al'),
            'C': symbols.count('C'),
            'O': symbols.count('O'),
            'H': symbols.count('H'),
            'Si': symbols.count('Si')
        }
        
        print(f"\n📊 Cycle Summary:")
        print(f"   TMA reactions: {tma_reactions}")
        print(f"   H2O reactions: {h2o_reactions}")
        print(f"   Total atoms: {len(self.surface)}")
        print(f"   Composition: Al:{composition['Al']}, C:{composition['C']}, O:{composition['O']}")
        
        if composition['Al'] > 0:
            c_per_al = composition['C'] / composition['Al']
            print(f"   C/Al ratio: {c_per_al:.1f} (should decrease over cycles)")
        
        return tma_reactions, h2o_reactions
    
    def run_corrected_mld(self, n_cycles: int = 5):
        """Run corrected MLD process"""
        print(f"\n🚀 Starting Corrected MLD Process")
        print(f"   Expected behavior:")
        print(f"   - Cycle 1: C/Al ≈ 2 (Al(CH3)2 formation)")
        print(f"   - Later cycles: C/Al decreases (CH4 elimination)")
        print(f"   - Al-O-Al bridge formation")
        
        for i in range(n_cycles):
            tma_count, h2o_count = self.run_cycle()
            
            if tma_count == 0 and h2o_count == 0:
                print(f"\n⚠️  No reactions in cycle {i+1}, stopping")
                break
        
        # Final analysis
        write('corrected_mld_final.xyz', self.surface)
        
        symbols = self.surface.get_chemical_symbols()
        final_comp = {
            'Al': symbols.count('Al'),
            'C': symbols.count('C'),
            'O': symbols.count('O'),
            'Si': symbols.count('Si')
        }
        
        print(f"\n{'='*60}")
        print(f"✅ Corrected MLD Complete!")
        print(f"   Cycles: {self.cycle_count}")
        print(f"   Final composition: {final_comp}")
        
        if final_comp['Al'] > 0:
            final_c_per_al = final_comp['C'] / final_comp['Al']
            print(f"   Final C/Al: {final_c_per_al:.1f}")
            
            if final_c_per_al < 1.5:
                print(f"   ✅ Good methyl elimination (low C/Al)")
            else:
                print(f"   ⚠️  High C/Al - may need more H2O reactions")
        
        return self.surface

def main():
    """Run corrected MLD simulation"""
    
    print("🔧 Corrected MLD Simulation")
    print("Fixes methyl accumulation and implements proper Al-O-Al growth")
    print("=" * 70)
    
    # Load surface
    surface_files = ['uv_ozone_si_surface.xyz', 'surface_si100_hydroxylated.xyz', 'simple_surface.xyz']
    surface = None
    
    for surf_file in surface_files:
        if os.path.exists(surf_file):
            surface = read(surf_file)
            print(f"✅ Loaded surface: {surf_file} ({len(surface)} atoms)")
            break
    
    if surface is None:
        print("❌ No surface found! Run surface creation script first.")
        return 1
    
    # Create corrected simulator
    simulator = CorrectedMLDSimulator(surface)
    
    # Run corrected MLD
    final_surface = simulator.run_corrected_mld(n_cycles=5)
    
    print(f"\n🎉 Corrected MLD simulation complete!")
    print(f"   Output files: corrected_*.xyz")
    print(f"   Use surface_visualizer.py to analyze results")
    
    return 0

if __name__ == "__main__":
    import os
    exit(main())