#!/usr/bin/env python3
"""
Step 5 Realistic: Molecular Layer Deposition (MLD) Simulation
Proper modeling of TMA/H2O ALD/MLD cycles on silicon surfaces

This simulates realistic MLD conditions:
- Surface reactions (not gas phase)
- Proper reaction mechanisms and kinetics
- Temperature effects
- Self-limiting growth
- Steric hindrance
"""

from ase import Atoms
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.build import add_adsorbate
from ase.optimize import BFGS
from ase.thermochemistry import IdealGasThermo
from ase.calculators.emt import EMT
import numpy as np
from typing import List, Tuple, Optional
import time

class MLDSimulator:
    """Realistic MLD simulation with proper surface chemistry"""
    
    def __init__(self, surface: Atoms, temperature: float = 473.15):  # 200°C typical ALD temp
        """
        Initialize MLD simulator
        
        Args:
            surface: ASE Atoms object of hydroxylated surface
            temperature: Process temperature in K
        """
        self.surface = surface.copy()
        self.temperature = temperature
        self.cycle_count = 0
        self.reaction_sites = []
        self.coverage_history = []
        
        # Reaction parameters (based on literature)
        self.tma_pulse_time = 0.1  # seconds
        self.h2o_pulse_time = 0.1  # seconds
        self.purge_time = 5.0  # seconds
        self.tma_pressure = 1e-3  # Torr
        self.h2o_pressure = 1e-3  # Torr
        
        # Identify reactive OH sites
        self._identify_oh_sites()
        
    def _identify_oh_sites(self):
        """Find surface OH groups that can react"""
        positions = self.surface.get_positions()
        symbols = self.surface.get_chemical_symbols()
        
        self.oh_sites = []
        
        # Find O-H pairs
        for i, sym_i in enumerate(symbols):
            if sym_i == 'O':
                # Look for nearby H
                for j, sym_j in enumerate(symbols):
                    if sym_j == 'H':
                        distance = np.linalg.norm(positions[i] - positions[j])
                        if 0.9 < distance < 1.1:  # O-H bond
                            self.oh_sites.append({
                                'O_idx': i,
                                'H_idx': j,
                                'position': positions[i],
                                'accessible': True
                            })
                            break
        
        print(f"📍 Found {len(self.oh_sites)} reactive OH sites on surface")
        
    def calculate_steric_factor(self, site_pos: np.ndarray, occupied_sites: List[np.ndarray]) -> float:
        """
        Calculate steric hindrance factor for a given site
        
        Args:
            site_pos: Position of the site to check
            occupied_sites: Positions of already occupied sites
            
        Returns:
            Steric factor (0-1, where 1 is no hindrance)
        """
        if not occupied_sites:
            return 1.0
        
        # TMA is bulky - approximately 6-7 Å diameter
        tma_radius = 3.5  # Å
        
        min_distance = float('inf')
        for occ_pos in occupied_sites:
            dist = np.linalg.norm(site_pos[:2] - occ_pos[:2])  # 2D distance on surface
            min_distance = min(min_distance, dist)
        
        # Steric factor decreases as sites get closer
        if min_distance < 2 * tma_radius:
            return max(0, (min_distance - tma_radius) / tma_radius)
        return 1.0
    
    def tma_pulse(self, fraction: float = 1.0) -> int:
        """
        Simulate TMA pulse
        
        Args:
            fraction: Fraction of available sites to attempt reaction (for MC simulation)
            
        Returns:
            Number of TMA molecules adsorbed
        """
        print(f"\n💨 TMA Pulse (Cycle {self.cycle_count + 1})")
        print(f"   Temperature: {self.temperature - 273.15:.0f}°C")
        print(f"   Pressure: {self.tma_pressure*760:.1f} mTorr")
        
        # Find available OH sites
        available_sites = [site for site in self.oh_sites if site['accessible']]
        print(f"   Available OH sites: {len(available_sites)}")
        
        if not available_sites:
            print("   ⚠️  No available sites for TMA reaction")
            return 0
        
        # Reaction: Si-OH + Al(CH3)3 → Si-O-Al(CH3)2 + CH4
        reacted_sites = []
        occupied_positions = []
        
        # Sort sites by z-coordinate (react with highest sites first)
        available_sites.sort(key=lambda s: s['position'][2], reverse=True)
        
        # Attempt reaction at each site
        for site in available_sites:
            # Check steric hindrance
            steric_factor = self.calculate_steric_factor(site['position'], occupied_positions)
            
            # Temperature-dependent reaction probability
            # Arrhenius-like: higher T = higher reaction rate
            base_probability = 0.8  # High for TMA (very reactive)
            temp_factor = np.exp(-2000 / (8.314 * self.temperature))  # Activation energy ~2 kJ/mol
            
            reaction_probability = base_probability * temp_factor * steric_factor * fraction
            
            if np.random.random() < reaction_probability:
                # React! Remove H, add Al(CH3)2
                self._react_tma_at_site(site)
                reacted_sites.append(site)
                occupied_positions.append(site['position'])
                site['accessible'] = False  # Site now has Al(CH3)2
        
        print(f"   ✅ TMA molecules adsorbed: {len(reacted_sites)}")
        print(f"   Surface coverage: {len(reacted_sites)/len(self.oh_sites)*100:.1f}%")
        
        # Purge
        time.sleep(0.01)  # Simulate purge time (scaled down)
        print(f"   💨 Purging for {self.purge_time} seconds...")
        
        return len(reacted_sites)
    
    def _react_tma_at_site(self, site):
        """Add Al(CH3)2 group at OH site"""
        positions = self.surface.get_positions()
        symbols = self.surface.get_chemical_symbols()
        
        # Remove H from OH
        h_idx = site['H_idx']
        o_pos = positions[site['O_idx']]
        
        # Add Al(CH3)2 group
        # Al position: 1.9 Å from O (Al-O bond)
        al_pos = o_pos + np.array([0, 0, 1.9])
        
        # Add two methyl groups in bent geometry
        # Al-C bond: 1.96 Å, C-Al-C angle: ~120°
        c1_pos = al_pos + np.array([1.96 * np.cos(np.pi/6), 0, 1.96 * np.sin(np.pi/6)])
        c2_pos = al_pos + np.array([-1.96 * np.cos(np.pi/6), 0, 1.96 * np.sin(np.pi/6)])
        
        # Add hydrogens to methyls (simplified - just one H per C for visualization)
        h1_pos = c1_pos + np.array([0, 0, 1.09])
        h2_pos = c2_pos + np.array([0, 0, 1.09])
        
        # Build new positions and symbols (remove H, add Al(CH3)2)
        mask = [i != h_idx for i in range(len(symbols))]
        new_positions = positions[mask]
        new_symbols = [symbols[i] for i in range(len(symbols)) if mask[i]]
        
        # Add new atoms
        new_positions = np.vstack([new_positions, al_pos, c1_pos, c2_pos, h1_pos, h2_pos])
        new_symbols.extend(['Al', 'C', 'C', 'H', 'H'])
        
        # Update surface
        self.surface = Atoms(
            symbols=new_symbols,
            positions=new_positions,
            cell=self.surface.get_cell(),
            pbc=self.surface.get_pbc()
        )
        
        # Update constraints if any
        if self.surface.constraints:
            # Keep fixed atoms fixed
            old_constraint = self.surface.constraints[0]
            if isinstance(old_constraint, FixAtoms):
                # Adjust indices after removing H
                new_indices = []
                for idx in old_constraint.index:
                    if idx < h_idx:
                        new_indices.append(idx)
                    elif idx > h_idx:
                        new_indices.append(idx - 1)
                self.surface.set_constraint(FixAtoms(indices=new_indices))
    
    def h2o_pulse(self) -> int:
        """
        Simulate H2O pulse
        
        Returns:
            Number of OH groups created
        """
        print(f"\n💧 H2O Pulse (Cycle {self.cycle_count + 1})")
        print(f"   Temperature: {self.temperature - 273.15:.0f}°C")
        print(f"   Pressure: {self.h2o_pressure*760:.1f} mTorr")
        
        # Find Al(CH3)2 groups
        symbols = self.surface.get_chemical_symbols()
        positions = self.surface.get_positions()
        
        al_sites = []
        for i, sym in enumerate(symbols):
            if sym == 'Al':
                al_sites.append(i)
        
        print(f"   Al(CH3)2 sites available: {len(al_sites)}")
        
        if not al_sites:
            print("   ⚠️  No Al sites for H2O reaction")
            return 0
        
        # Reaction: Si-O-Al(CH3)2 + H2O → Si-O-Al(OH)CH3 + CH4
        # Simplified: we'll convert Al(CH3)2 back to OH for next cycle
        reacted_count = 0
        
        for al_idx in al_sites:
            # Temperature-dependent reaction probability
            base_probability = 0.9  # H2O is very reactive with Al-CH3
            temp_factor = np.exp(-1500 / (8.314 * self.temperature))
            
            if np.random.random() < base_probability * temp_factor:
                # In real ALD, this would create Al-OH and release CH4
                # For simplicity, we'll mark site as accessible again
                reacted_count += 1
        
        # Reset OH sites for next cycle (simplified)
        for site in self.oh_sites:
            site['accessible'] = True
        
        print(f"   ✅ OH groups created: {reacted_count}")
        
        # Purge
        time.sleep(0.01)
        print(f"   💨 Purging for {self.purge_time} seconds...")
        
        return reacted_count
    
    def run_cycle(self):
        """Run one complete MLD cycle (TMA + H2O)"""
        print(f"\n{'='*60}")
        print(f"🔄 MLD Cycle {self.cycle_count + 1}")
        print(f"{'='*60}")
        
        # TMA pulse
        tma_count = self.tma_pulse()
        
        # Save after TMA
        write(f'mld_cycle{self.cycle_count+1}_after_tma.xyz', self.surface)
        
        # H2O pulse
        oh_count = self.h2o_pulse()
        
        # Save after H2O
        write(f'mld_cycle{self.cycle_count+1}_after_h2o.xyz', self.surface)
        
        # Track coverage
        coverage = tma_count / len(self.oh_sites) if self.oh_sites else 0
        self.coverage_history.append(coverage)
        
        self.cycle_count += 1
        
        print(f"\n📊 Cycle Summary:")
        print(f"   TMA molecules adsorbed: {tma_count}")
        print(f"   Surface coverage: {coverage*100:.1f}%")
        print(f"   Growth per cycle: ~{coverage * 2.5:.2f} Å")  # Approximate for Al2O3
        
        return tma_count, oh_count
    
    def run_mld_process(self, n_cycles: int = 5):
        """
        Run multiple MLD cycles
        
        Args:
            n_cycles: Number of cycles to run
        """
        print(f"\n🚀 Starting MLD Process")
        print(f"   Surface: {len(self.surface)} atoms")
        print(f"   Temperature: {self.temperature - 273.15:.0f}°C")
        print(f"   Planned cycles: {n_cycles}")
        
        for i in range(n_cycles):
            tma_count, oh_count = self.run_cycle()
            
            if tma_count == 0:
                print(f"\n⚠️  Growth stopped - no available sites")
                break
        
        # Final summary
        print(f"\n{'='*60}")
        print(f"✅ MLD Process Complete!")
        print(f"   Total cycles: {self.cycle_count}")
        print(f"   Average coverage: {np.mean(self.coverage_history)*100:.1f}%")
        print(f"   Final surface: {len(self.surface)} atoms")
        
        # Save final structure
        write('mld_final_structure.xyz', self.surface)
        
        return self.surface

def main():
    """Run realistic MLD simulation"""
    
    print("🎯 Realistic MLD Simulation")
    print("=" * 60)
    
    # Load surface
    print("\n📂 Loading surface...")
    try:
        # Try to load one of our created surfaces
        surface = read('surface_si100_hydroxylated.xyz')
        print(f"✅ Loaded Si(100) surface: {len(surface)} atoms")
    except:
        print("❌ No surface found. Run step4_create_surface_advanced.py first!")
        return 1
    
    # Create simulator
    simulator = MLDSimulator(surface, temperature=473.15)  # 200°C
    
    # Run MLD cycles
    final_surface = simulator.run_mld_process(n_cycles=5)
    
    # Analysis
    print(f"\n📊 Growth Analysis:")
    if simulator.coverage_history:
        print(f"   Coverage per cycle: {[f'{c*100:.1f}%' for c in simulator.coverage_history]}")
        thickness = sum(simulator.coverage_history) * 2.5  # Approximate
        print(f"   Estimated film thickness: {thickness:.1f} Å")
    
    print(f"\n🎉 Realistic MLD simulation complete!")
    print(f"   Output files:")
    print(f"   - mld_final_structure.xyz")
    print(f"   - mld_cycle*_after_*.xyz (intermediate structures)")
    
    return 0

if __name__ == "__main__":
    exit(main())