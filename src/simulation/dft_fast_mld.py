#!/usr/bin/env python3
"""
Fast DFT MLD Simulator with Pre-Optimized Structures
Skips redundant optimizations by using pre-calculated structures and energies

Key features:
- Uses pre-optimized molecules (TMA, H2O)
- Optional surface optimization skip
- Loads pre-calculated energies from files
- Focuses on reaction steps only
"""

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS, LBFGS
from ase.constraints import FixAtoms
import time
import os
import json
from datetime import datetime
from typing import Tuple, Optional, Dict, List

# Import required modules
try:
    from gpaw import GPAW, PW
    GPAW_AVAILABLE = True
except ImportError:
    GPAW_AVAILABLE = False
    print("❌ GPAW not available! Install with: pip install gpaw")

from geometry_validator import GeometryValidator

class FastDFTMLDSimulator:
    """Fast DFT MLD simulator using pre-optimized structures"""
    
    def __init__(self, max_atoms: int = 200):
        """Initialize with pre-optimized structure support"""
        self.MAX_ATOMS = max_atoms
        self.WARN_ATOMS = max_atoms * 0.75
        
        # Conservative DFT parameters
        self.dft_params = {
            'mode': PW(300),
            'xc': 'PBE',
            'kpts': (1, 1, 1),
            'txt': None,
            'symmetry': 'off',
            'convergence': {
                'energy': 1e-5,
                'density': 1e-4,
                'eigenstates': 1e-4
            },
            'mixer': {'beta': 0.1},
            'maxiter': 300,
            'nbands': 'nao'
        }
        
        # Optimization parameters
        self.opt_params = {
            'fmax': 0.03,
            'max_steps': 100,
            'optimizer': LBFGS
        }
        
        # Pre-calculated energies storage
        self.energy_cache_file = "pre_calculated_energies.json"
        self.energy_cache = self.load_energy_cache()
        
        # Initialize geometry validator
        self.validator = GeometryValidator()
        
        print("🚀 Fast DFT MLD Simulator")
        print(f"   Maximum atoms: {self.MAX_ATOMS}")
        print(f"   Using pre-optimized structures when available")
        
    def load_energy_cache(self) -> Dict[str, float]:
        """Load pre-calculated energies from file"""
        if os.path.exists(self.energy_cache_file):
            try:
                with open(self.energy_cache_file, 'r') as f:
                    cache = json.load(f)
                print(f"   📂 Loaded {len(cache)} pre-calculated energies")
                return cache
            except:
                return {}
        return {}
    
    def save_energy_cache(self):
        """Save energy cache to file"""
        with open(self.energy_cache_file, 'w') as f:
            json.dump(self.energy_cache, f, indent=2)
    
    def get_cached_energy(self, label: str) -> Optional[float]:
        """Get pre-calculated energy if available"""
        if label in self.energy_cache:
            print(f"   ⚡ Using cached energy for {label}: {self.energy_cache[label]:.4f} eV")
            return self.energy_cache[label]
        return None
    
    def validate_system_size(self, atoms: Atoms, name: str) -> bool:
        """Validate system is appropriate for DFT"""
        n_atoms = len(atoms)
        
        print(f"\n🔍 Size validation: {name}")
        print(f"   Atoms: {n_atoms}")
        
        if n_atoms > self.MAX_ATOMS:
            print(f"   ❌ TOO LARGE: {n_atoms} > {self.MAX_ATOMS}")
            return False
        elif n_atoms > self.WARN_ATOMS:
            print(f"   ⚠️  Large system: {n_atoms} > {self.WARN_ATOMS:.0f}")
        else:
            print(f"   ✅ Good size: {n_atoms} atoms")
        
        return True
    
    def create_dft_calculator(self, atoms: Atoms, label: str) -> Optional[GPAW]:
        """Create GPAW calculator with appropriate settings"""
        if not GPAW_AVAILABLE:
            return None
        
        # Adjust vacuum
        positions = atoms.get_positions()
        if len(positions) > 0:
            size = positions.max(axis=0) - positions.min(axis=0)
            vacuum_needed = max(10.0, 0.3 * size.max())
            
            atoms_centered = atoms.copy()
            atoms_centered.center(vacuum=vacuum_needed)
            
            print(f"   Vacuum: {vacuum_needed:.1f} Å")
        
        # Create calculator
        params = self.dft_params.copy()
        params['txt'] = f'{label}_dft.txt'
        
        try:
            calc = GPAW(**params)
            print(f"   ✅ DFT calculator created")
            return calc
        except Exception as e:
            print(f"   ❌ Calculator creation failed: {e}")
            return None
    
    def optimize_structure(self, atoms: Atoms, label: str, 
                         skip_if_exists: bool = True) -> Tuple[Atoms, float, bool]:
        """
        Optimize structure with option to skip if already done
        """
        print(f"\n🔧 DFT optimization: {label}")
        
        # Check if we should skip
        if skip_if_exists:
            # Check for existing optimized file
            opt_file = f"{label}_optimized.xyz"
            if os.path.exists(opt_file):
                print(f"   📂 Found existing optimized structure: {opt_file}")
                opt_atoms = read(opt_file)
                
                # Check for cached energy
                cached_energy = self.get_cached_energy(label)
                if cached_energy is not None:
                    return opt_atoms, cached_energy, True
                
                # Calculate energy for existing structure
                print(f"   ⚡ Calculating energy for existing structure...")
                calc = self.create_dft_calculator(opt_atoms, label)
                if calc:
                    opt_atoms.calc = calc
                    try:
                        energy = opt_atoms.get_potential_energy()
                        self.energy_cache[label] = energy
                        self.save_energy_cache()
                        print(f"   ✅ Energy: {energy:.4f} eV")
                        return opt_atoms, energy, True
                    except Exception as e:
                        print(f"   ❌ Energy calculation failed: {e}")
        
        # Size check
        if not self.validate_system_size(atoms, label):
            return atoms.copy(), 0.0, False
        
        # Normal optimization
        calc = self.create_dft_calculator(atoms, label)
        if not calc:
            return atoms.copy(), 0.0, False
        
        atoms.calc = calc
        
        # Set up optimization
        opt_log = f'{label}_opt.log'
        if self.opt_params['optimizer'] == LBFGS:
            opt = LBFGS(atoms, logfile=opt_log)
        else:
            opt = BFGS(atoms, logfile=opt_log)
        
        print(f"   📝 Log: {opt_log}")
        print(f"   🎯 Target fmax: {self.opt_params['fmax']} eV/Å")
        
        # Run optimization
        start_time = time.time()
        try:
            converged = opt.run(fmax=self.opt_params['fmax'], 
                              steps=self.opt_params['max_steps'])
            
            fmax = opt.get_residual()
            energy = atoms.get_potential_energy()
            elapsed = time.time() - start_time
            
            print(f"   ✅ Optimization complete")
            print(f"   Energy: {energy:.4f} eV")
            print(f"   Max force: {fmax:.4f} eV/Å")
            print(f"   Converged: {converged}")
            print(f"   Time: {elapsed/60:.1f} minutes")
            
            # Save optimized structure and energy
            write(f'{label}_optimized.xyz', atoms)
            self.energy_cache[label] = energy
            self.save_energy_cache()
            
            return atoms.copy(), energy, converged
            
        except Exception as e:
            elapsed = time.time() - start_time
            print(f"   ❌ Optimization failed: {e}")
            print(f"   Time: {elapsed/60:.1f} minutes")
            return atoms.copy(), 0.0, False
    
    def load_optimized_molecule(self, molecule_type: str) -> Optional[Tuple[Atoms, float]]:
        """Load pre-optimized molecule and its energy"""
        molecule_files = {
            'tma': ['trimethylaluminum_optimized.xyz', 'tma_optimized.xyz'],
            'h2o': ['water_optimized.xyz', 'h2o_optimized.xyz'],
            'water': ['water_optimized.xyz', 'h2o_optimized.xyz']
        }
        
        if molecule_type.lower() not in molecule_files:
            return None
        
        for mol_file in molecule_files[molecule_type.lower()]:
            if os.path.exists(mol_file):
                print(f"\n📂 Loading pre-optimized {molecule_type}: {mol_file}")
                mol = read(mol_file)
                
                # Get energy
                cached_energy = self.get_cached_energy(molecule_type)
                if cached_energy:
                    return mol, cached_energy
                
                # Calculate energy if not cached
                print(f"   Calculating energy...")
                calc = self.create_dft_calculator(mol, molecule_type)
                if calc:
                    mol.calc = calc
                    try:
                        energy = mol.get_potential_energy()
                        self.energy_cache[molecule_type] = energy
                        self.save_energy_cache()
                        print(f"   ✅ Energy: {energy:.4f} eV")
                        return mol, energy
                    except:
                        pass
                
                return mol, 0.0
        
        return None
    
    def run_mld_cycle(self, surface_file: str, skip_surface_opt: bool = True):
        """
        Run MLD cycle with pre-optimized molecules
        
        Args:
            surface_file: Path to surface structure
            skip_surface_opt: Skip surface optimization if True
        """
        print("\n" + "="*60)
        print("FAST MLD CYCLE SIMULATION")
        print("="*60)
        
        # Load surface
        if not os.path.exists(surface_file):
            print(f"❌ Surface file not found: {surface_file}")
            return
        
        surface = read(surface_file)
        print(f"\n📂 Loaded surface: {surface_file}")
        print(f"   Atoms: {len(surface)}")
        
        # Surface optimization (optional)
        if skip_surface_opt:
            print("\n⚡ Skipping surface optimization (using as-is)")
            surface_opt = surface
            surface_energy = self.get_cached_energy('surface') or 0.0
            if surface_energy == 0.0:
                print("   ⚠️  No cached surface energy available")
        else:
            surface_opt, surface_energy, _ = self.optimize_structure(
                surface, "surface", skip_if_exists=True
            )
        
        # Load pre-optimized TMA
        tma_data = self.load_optimized_molecule('tma')
        if not tma_data:
            print("❌ No pre-optimized TMA found!")
            return
        tma_opt, tma_energy = tma_data
        
        # Load pre-optimized H2O
        h2o_data = self.load_optimized_molecule('h2o')
        if not h2o_data:
            print("❌ No pre-optimized H2O found!")
            return
        h2o_opt, h2o_energy = h2o_data
        
        # Now focus on the actual MLD reactions
        print("\n" + "="*60)
        print("MLD REACTION STEPS")
        print("="*60)
        
        # Step 1: TMA adsorption
        print("\n🧪 Step 1: TMA Adsorption")
        # Position TMA above surface
        tma_positioned = tma_opt.copy()
        tma_positioned.translate([0, 0, 4.0])  # 4Å above surface
        
        # Combine and optimize
        combined_tma = surface_opt + tma_positioned
        if self.validate_system_size(combined_tma, "Surface+TMA"):
            print("   Optimizing TMA adsorption...")
            cycle1_tma, cycle1_tma_energy, _ = self.optimize_structure(
                combined_tma, "cycle1_tma_adsorbed", skip_if_exists=False
            )
            
            # Calculate binding energy
            if surface_energy != 0 and tma_energy != 0 and cycle1_tma_energy != 0:
                binding_tma = cycle1_tma_energy - surface_energy - tma_energy
                print(f"   💡 TMA binding energy: {binding_tma:.4f} eV")
        
        print("\n✅ Fast MLD simulation complete!")
        print("   Check *_optimized.xyz files for structures")
        print("   Check pre_calculated_energies.json for energies")

def main():
    """Main execution"""
    print("🚀 Fast DFT MLD Simulator")
    print("   Uses pre-optimized structures to save time")
    
    # Initialize simulator
    simulator = FastDFTMLDSimulator(max_atoms=200)
    
    # Check for available surfaces (prioritize smaller ones)
    surface_options = [
        'small_si100_2x2_4L.xyz',        # 128 atoms (preferred)
        'small_si100_2x2_3L_OH.xyz',     # ~80 atoms with OH (backup)
        'small_si100_2x2_4L_OH.xyz',     # 160 atoms (last resort)
        'small_si100_2x2_6L.xyz'
    ]
    
    surface_file = None
    for surf in surface_options:
        if os.path.exists(surf):
            surface_file = surf
            break
    
    if not surface_file:
        print("❌ No suitable surface file found!")
        print("   Run create_small_surfaces.py first")
        return
    
    # Run simulation
    simulator.run_mld_cycle(
        surface_file=surface_file,
        skip_surface_opt=True  # Skip surface optimization
    )

if __name__ == "__main__":
    main()