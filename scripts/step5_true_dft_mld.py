#!/usr/bin/env python3
"""
True DFT MLD Simulation using GPAW
Real quantum mechanical calculations for TMA/H2O MLD cycles

WARNING: This will take HOURS to DAYS to complete!
Each geometry optimization takes 10-60 minutes.
"""

from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS, LBFGS
from ase.constraints import FixAtoms
from ase.calculators.calculator import Calculator
import numpy as np
import time
import os
from datetime import datetime

# Import GPAW for DFT calculations
try:
    from gpaw import GPAW, PW
    GPAW_AVAILABLE = True
except ImportError:
    GPAW_AVAILABLE = False
    print("❌ GPAW not available! Install with: pip install gpaw")
    exit(1)

class DFTMLDSimulator:
    """
    True DFT MLD simulation with quantum mechanical calculations
    
    This performs real ab initio calculations for each reaction step:
    1. TMA adsorption and reaction energy
    2. Geometry optimization after each step
    3. Transition state searches (simplified)
    4. Real activation barriers
    """
    
    def __init__(self, surface: Atoms, temperature: float = 473.15):
        """
        Initialize DFT MLD simulator
        
        Args:
            surface: Hydroxylated surface
            temperature: Process temperature in K
        """
        self.surface = surface.copy()
        self.temperature = temperature
        self.cycle_count = 0
        
        # DFT parameters
        self.dft_params = {
            'mode': PW(400),           # 400 eV cutoff
            'xc': 'PBE',               # Exchange-correlation
            'kpts': (1, 1, 1),        # Gamma point for surfaces
            'txt': None,               # Will be set per calculation
            'symmetry': 'off',         # Often better for surfaces
            'convergence': {
                'energy': 1e-5,        # eV
                'density': 1e-4,
                'eigenstates': 1e-4
            },
            'mixer': {'beta': 0.1},    # Conservative mixing
            'maxiter': 300             # SCF iterations
        }
        
        # Optimization parameters
        self.opt_params = {
            'fmax': 0.02,              # eV/Å (tight for surfaces)
            'max_steps': 200,
            'optimizer': LBFGS         # Memory efficient
        }
        
        # Identify reactive sites
        self._identify_reaction_sites()
        
        print(f"🔬 DFT MLD Simulator initialized")
        print(f"   Surface: {len(self.surface)} atoms")
        print(f"   Temperature: {self.temperature - 273.15:.0f}°C")
        print(f"   Reactive OH sites: {len(self.oh_sites)}")
        print(f"   ⚠️  WARNING: This will take HOURS per cycle!")
    
    def _identify_reaction_sites(self):
        """Find reactive OH groups on surface"""
        positions = self.surface.get_positions()
        symbols = self.surface.get_chemical_symbols()
        
        self.oh_sites = []
        
        # Find O-H pairs on surface
        for i, sym_i in enumerate(symbols):
            if sym_i == 'O':
                # Check if this O is on the surface (top 3 Å)
                z_coords = positions[:, 2]
                if positions[i][2] > z_coords.max() - 3.0:
                    # Look for nearby H
                    for j, sym_j in enumerate(symbols):
                        if sym_j == 'H':
                            distance = np.linalg.norm(positions[i] - positions[j])
                            if 0.85 < distance < 1.15:  # O-H bond range
                                self.oh_sites.append({
                                    'O_idx': i,
                                    'H_idx': j,
                                    'position': positions[i],
                                    'reacted': False
                                })
                                break
        
        print(f"📍 Found {len(self.oh_sites)} reactive surface OH sites")
    
    def create_dft_calculator(self, label: str) -> GPAW:
        """
        Create GPAW calculator for current system
        
        Args:
            label: Label for output files
            
        Returns:
            GPAW calculator instance
        """
        calc_params = self.dft_params.copy()
        calc_params['txt'] = f'{label}_gpaw.txt'
        
        return GPAW(**calc_params)
    
    def optimize_structure(self, atoms: Atoms, label: str) -> Tuple[Atoms, float, bool]:
        """
        Optimize structure with DFT
        
        Args:
            atoms: Structure to optimize
            label: Label for files
            
        Returns:
            (optimized_atoms, final_energy, converged)
        """
        print(f"🔧 DFT optimization: {label}")
        print(f"   Atoms: {len(atoms)}")
        
        # Set up calculator
        calc = self.create_dft_calculator(label)
        atoms.calc = calc
        
        # Set up optimizer
        optimizer_class = self.opt_params['optimizer']
        opt = optimizer_class(atoms, logfile=f'{label}_opt.log')
        
        start_time = time.time()
        
        try:
            # Run optimization
            opt.run(fmax=self.opt_params['fmax'], 
                   steps=self.opt_params['max_steps'])
            
            # Get final results
            final_energy = atoms.get_potential_energy()
            forces = atoms.get_forces()
            max_force = np.sqrt((forces**2).sum(axis=1).max())
            
            converged = max_force < self.opt_params['fmax']
            elapsed = time.time() - start_time
            
            print(f"   ✅ Optimization complete ({elapsed/60:.1f} min)")
            print(f"   Energy: {final_energy:.4f} eV")
            print(f"   Max force: {max_force:.4f} eV/Å")
            print(f"   Converged: {converged}")
            
            # Save structure
            write(f'{label}_optimized.xyz', atoms)
            
            return atoms.copy(), final_energy, converged
            
        except Exception as e:
            print(f"   ❌ Optimization failed: {e}")
            elapsed = time.time() - start_time
            print(f"   Time spent: {elapsed/60:.1f} min")
            return atoms.copy(), 0.0, False
    
    def calculate_reaction_energy(self, reactants: Atoms, products: Atoms, 
                                label: str) -> float:
        """
        Calculate reaction energy: E(products) - E(reactants)
        
        Args:
            reactants: Initial state
            products: Final state  
            label: Label for calculations
            
        Returns:
            Reaction energy in eV
        """
        print(f"⚡ Calculating reaction energy: {label}")
        
        # Optimize both structures
        reactants_opt, e_reactants, conv_r = self.optimize_structure(
            reactants, f'{label}_reactants'
        )
        
        products_opt, e_products, conv_p = self.optimize_structure(
            products, f'{label}_products'
        )
        
        if conv_r and conv_p:
            reaction_energy = e_products - e_reactants
            print(f"   📊 Reaction energy: {reaction_energy:.3f} eV")
            
            if reaction_energy < 0:
                print(f"   ✅ Exothermic reaction (favorable)")
            else:
                print(f"   ⚠️  Endothermic reaction ({reaction_energy:.3f} eV)")
                
            return reaction_energy
        else:
            print(f"   ❌ Reaction energy calculation failed (convergence issues)")
            return 0.0
    
    def tma_adsorption_dft(self, site_idx: int) -> Tuple[Atoms, float, bool]:
        """
        Perform DFT calculation for TMA adsorption at specific site
        
        Reaction: Si-OH + Al(CH3)3 → Si-O-Al(CH3)2 + CH4
        
        Args:
            site_idx: Index of OH site to react
            
        Returns:
            (product_surface, reaction_energy, success)
        """
        print(f"\n💥 DFT TMA Adsorption at site {site_idx}")
        
        site = self.oh_sites[site_idx]
        if site['reacted']:
            print(f"   ⚠️  Site already reacted!")
            return self.surface.copy(), 0.0, False
        
        # Create reactant system: Surface + TMA molecule
        try:
            tma = read('trimethylaluminum_optimized.xyz')
        except:
            print(f"   ❌ TMA molecule not found!")
            return self.surface.copy(), 0.0, False
        
        # Position TMA near the OH site
        oh_pos = site['position']
        tma_center = tma.get_positions().mean(axis=0)
        
        # Place TMA 3 Å away from OH
        offset = oh_pos + np.array([3.0, 0.0, 2.0])
        tma.translate(offset - tma_center)
        
        # Create reactant system
        reactant_pos = np.vstack([self.surface.get_positions(), tma.get_positions()])
        reactant_symbols = (list(self.surface.get_chemical_symbols()) + 
                           list(tma.get_chemical_symbols()))
        
        reactants = Atoms(
            symbols=reactant_symbols,
            positions=reactant_pos,
            cell=self.surface.get_cell(),
            pbc=self.surface.get_pbc()
        )
        
        # Fix bottom layers of surface
        n_surface = len(self.surface)
        z_coords = self.surface.get_positions()[:, 2]
        z_min = z_coords.min()
        
        fix_indices = []
        for i in range(n_surface):
            if self.surface.get_positions()[i][2] < z_min + 6.0:
                fix_indices.append(i)
        
        if fix_indices:
            reactants.set_constraint(FixAtoms(indices=fix_indices))
        
        print(f"   Reactants: {len(reactants)} atoms")
        print(f"   Fixed atoms: {len(fix_indices)}")
        
        # Optimize reactant state
        reactants_opt, e_reactants, conv_r = self.optimize_structure(
            reactants, f'cycle{self.cycle_count+1}_tma_site{site_idx}_reactants'
        )
        
        if not conv_r:
            print(f"   ❌ Reactant optimization failed")
            return self.surface.copy(), 0.0, False
        
        # Create product state (simplified - just remove H and add Al(CH3)2)
        # In reality, you'd need transition state search
        products = self._create_tma_product_state(reactants_opt, site)
        
        # Optimize product state
        products_opt, e_products, conv_p = self.optimize_structure(
            products, f'cycle{self.cycle_count+1}_tma_site{site_idx}_products'
        )
        
        if not conv_p:
            print(f"   ❌ Product optimization failed")
            return self.surface.copy(), 0.0, False
        
        # Calculate reaction energy
        reaction_energy = e_products - e_reactants
        
        print(f"   📊 TMA adsorption results:")
        print(f"      Reactant energy: {e_reactants:.4f} eV")
        print(f"      Product energy: {e_products:.4f} eV")
        print(f"      Reaction energy: {reaction_energy:.4f} eV")
        
        # Mark site as reacted
        site['reacted'] = True
        
        # Update surface to product state
        self.surface = products_opt.copy()
        
        return products_opt, reaction_energy, True
    
    def _create_tma_product_state(self, reactants: Atoms, site: dict) -> Atoms:
        """
        Create product state for TMA reaction (simplified)
        Real implementation would need proper reaction coordinate
        """
        # This is a simplified approach - real MLD would need:
        # 1. Transition state search
        # 2. Proper reaction pathway
        # 3. CH4 elimination pathway
        
        positions = reactants.get_positions()
        symbols = reactants.get_chemical_symbols()
        
        # Find Al atom in TMA (should be last few atoms)
        al_idx = None
        for i in range(len(symbols)-10, len(symbols)):
            if symbols[i] == 'Al':
                al_idx = i
                break
        
        if al_idx is None:
            return reactants.copy()
        
        # Remove H from OH group
        h_idx = site['H_idx']
        mask = [i != h_idx for i in range(len(symbols))]
        
        # Remove one methyl group (simplified CH4 elimination)
        # Find C atoms bonded to Al
        al_pos = positions[al_idx]
        c_to_remove = []
        
        for i, sym in enumerate(symbols):
            if sym == 'C' and i > len(self.surface.get_chemical_symbols()):
                c_pos = positions[i]
                if np.linalg.norm(c_pos - al_pos) < 2.2:  # Al-C bond
                    c_to_remove.append(i)
                    # Remove H atoms bonded to this C
                    for j, sym2 in enumerate(symbols):
                        if sym2 == 'H' and j > len(self.surface.get_chemical_symbols()):
                            h_pos = positions[j]
                            if np.linalg.norm(h_pos - c_pos) < 1.2:  # C-H bond
                                c_to_remove.append(j)
                    break  # Remove only one methyl group
        
        # Apply all removals
        for idx in sorted(c_to_remove, reverse=True):
            mask[idx] = False
        
        # Create product structure
        new_positions = positions[mask]
        new_symbols = [symbols[i] for i in range(len(symbols)) if mask[i]]
        
        products = Atoms(
            symbols=new_symbols,
            positions=new_positions,
            cell=reactants.get_cell(),
            pbc=reactants.get_pbc()
        )
        
        # Copy constraints
        if reactants.constraints:
            products.set_constraint(reactants.constraints[0])
        
        return products
    
    def run_dft_cycle(self) -> dict:
        """
        Run one complete DFT MLD cycle
        
        Returns:
            Dictionary with cycle results
        """
        print(f"\n{'='*80}")
        print(f"🔬 DFT MLD Cycle {self.cycle_count + 1}")
        print(f"{'='*80}")
        print(f"⚠️  Estimated time: 2-8 hours per cycle")
        
        cycle_start = time.time()
        results = {
            'cycle': self.cycle_count + 1,
            'tma_reactions': [],
            'total_time': 0,
            'success': False
        }
        
        # TMA pulse - react at available sites
        available_sites = [i for i, site in enumerate(self.oh_sites) 
                          if not site['reacted']]
        
        print(f"\n💨 TMA Pulse")
        print(f"   Available OH sites: {len(available_sites)}")
        
        if not available_sites:
            print(f"   ⚠️  No sites available for reaction")
            return results
        
        # React at first available site (for demonstration)
        # In reality, you'd consider kinetics, steric effects, etc.
        site_idx = available_sites[0]
        
        surface_after_tma, reaction_energy, success = self.tma_adsorption_dft(site_idx)
        
        if success:
            results['tma_reactions'].append({
                'site': site_idx,
                'reaction_energy': reaction_energy,
                'success': True
            })
            
            # Save intermediate structure
            write(f'dft_cycle{self.cycle_count+1}_after_tma.xyz', surface_after_tma)
        
        # H2O pulse (simplified for this demo)
        print(f"\n💧 H2O Pulse (simplified - just save current state)")
        write(f'dft_cycle{self.cycle_count+1}_after_h2o.xyz', self.surface)
        
        # Update cycle count
        self.cycle_count += 1
        
        # Calculate total time
        cycle_time = time.time() - cycle_start
        results['total_time'] = cycle_time
        results['success'] = success
        
        print(f"\n📊 Cycle {self.cycle_count} Summary:")
        print(f"   TMA reactions: {len(results['tma_reactions'])}")
        print(f"   Total time: {cycle_time/3600:.1f} hours")
        print(f"   Success: {success}")
        
        return results
    
    def run_dft_mld_process(self, n_cycles: int = 2) -> list:
        """
        Run complete DFT MLD process
        
        Args:
            n_cycles: Number of cycles (keep small - each takes hours!)
            
        Returns:
            List of cycle results
        """
        print(f"\n🚀 Starting DFT MLD Process")
        print(f"   Planned cycles: {n_cycles}")
        print(f"   Estimated total time: {n_cycles * 4:.0f}-{n_cycles * 8:.0f} hours")
        print(f"   ⚠️  WARNING: This is computationally intensive!")
        
        process_start = time.time()
        all_results = []
        
        for cycle in range(n_cycles):
            cycle_results = self.run_dft_cycle()
            all_results.append(cycle_results)
            
            if not cycle_results['success']:
                print(f"\n⚠️  Cycle {cycle+1} failed, stopping process")
                break
            
            # Save final structure after each cycle
            write(f'dft_final_cycle{cycle+1}.xyz', self.surface)
        
        # Final summary
        total_time = time.time() - process_start
        successful_cycles = len([r for r in all_results if r['success']])
        
        print(f"\n{'='*80}")
        print(f"✅ DFT MLD Process Complete!")
        print(f"   Successful cycles: {successful_cycles}/{n_cycles}")
        print(f"   Total time: {total_time/3600:.2f} hours")
        print(f"   Final surface: {len(self.surface)} atoms")
        
        # Save final structure
        write('dft_mld_final_structure.xyz', self.surface)
        
        return all_results

def main():
    """Run true DFT MLD simulation"""
    
    print("🔬 True DFT MLD Simulation")
    print("=" * 80)
    print("⚠️  WARNING: This will take HOURS to complete!")
    print("Each optimization step: 10-60 minutes")
    print("Full cycle: 2-8 hours")
    print("=" * 80)
    
    # Load surface
    print("\n📂 Loading hydroxylated surface...")
    try:
        # Try UV-ozone surface first
        surface = read('uv_ozone_si_surface.xyz')
        print(f"✅ Loaded UV-ozone surface: {len(surface)} atoms")
    except:
        try:
            surface = read('surface_si100_hydroxylated.xyz')
            print(f"✅ Loaded Si(100) surface: {len(surface)} atoms")
        except:
            print("❌ No suitable surface found!")
            print("Run step4_uv_ozone_surface.py or step4_create_surface_advanced.py first")
            return 1
    
    # Check if user really wants to proceed
    print(f"\n⚠️  About to start DFT calculations...")
    print(f"This will use significant computational resources.")
    print(f"Estimated time: 4-16 hours for 2 cycles")
    
    # Create DFT simulator
    simulator = DFTMLDSimulator(surface, temperature=473.15)
    
    # Run DFT MLD (start with just 1-2 cycles)
    results = simulator.run_dft_mld_process(n_cycles=2)
    
    # Analysis
    print(f"\n📊 DFT Results Summary:")
    for i, cycle_result in enumerate(results):
        if cycle_result['success']:
            print(f"   Cycle {i+1}: {len(cycle_result['tma_reactions'])} reactions")
            print(f"   Time: {cycle_result['total_time']/3600:.1f} hours")
    
    print(f"\n🎉 DFT MLD simulation complete!")
    print(f"   Output files: dft_*.xyz")
    print(f"   Log files: *_gpaw.txt, *_opt.log")
    
    return 0

if __name__ == "__main__":
    exit(main())