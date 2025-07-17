#!/usr/bin/env python3
"""
DFT Small System MLD Simulator
Proper DFT MLD simulation with literature-validated small systems

Prevents the "3171 atoms" failure by using appropriately sized systems:
- Maximum 200 atoms (hard limit)
- Literature-validated geometries
- Progressive complexity approach
- Comprehensive error handling

Based on literature best practices from:
- Halls & Raghavachari (2004)
- Mameli et al. (2017) 
- Shirazi & Elliott (2013)
"""

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS, LBFGS
from ase.constraints import FixAtoms
import time
import os
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

class SmallDFTMLDSimulator:
    """DFT MLD simulator for appropriately sized systems"""
    
    def __init__(self, max_atoms: int = 200):
        """
        Initialize with strict size limits
        
        Args:
            max_atoms: Maximum atoms allowed (literature: 50-200)
        """
        self.MAX_ATOMS = max_atoms
        self.WARN_ATOMS = max_atoms * 0.75  # Warning threshold
        
        # Conservative DFT parameters for reliability
        self.dft_params = {
            'mode': PW(300),           # Conservative cutoff
            'xc': 'PBE',               # Standard functional
            'kpts': (1, 1, 1),        # Gamma point
            'txt': None,               # Set per calculation
            'symmetry': 'off',         # Safer for surfaces
            'convergence': {
                'energy': 1e-5,        # Tight convergence
                'density': 1e-4,
                'eigenstates': 1e-4
            },
            'mixer': {'beta': 0.1},    # Conservative mixing
            'maxiter': 300,            # More SCF iterations
            'nbands': 'nao'           # Automatic bands
        }
        
        # Optimization parameters
        self.opt_params = {
            'fmax': 0.03,              # Tight force convergence
            'max_steps': 100,          # Conservative step limit
            'optimizer': LBFGS         # Memory efficient
        }
        
        # Initialize geometry validator
        self.validator = GeometryValidator()
        
        print("🔬 Small DFT MLD Simulator")
        print(f"   Maximum atoms: {self.MAX_ATOMS}")
        print(f"   Warning threshold: {self.WARN_ATOMS:.0f}")
        print(f"   Conservative DFT parameters")
        
        if not GPAW_AVAILABLE:
            print("   ❌ GPAW not available - simulation will fail")
        else:
            print("   ✅ GPAW available")
    
    def validate_system_size(self, atoms: Atoms, name: str) -> bool:
        """Validate system is appropriate for DFT"""
        
        n_atoms = len(atoms)
        
        print(f"\n🔍 Size validation: {name}")
        print(f"   Atoms: {n_atoms}")
        
        if n_atoms > self.MAX_ATOMS:
            print(f"   ❌ TOO LARGE: {n_atoms} > {self.MAX_ATOMS}")
            print(f"   This will cause DFT failure")
            return False
        elif n_atoms > self.WARN_ATOMS:
            print(f"   ⚠️  Large system: {n_atoms} > {self.WARN_ATOMS:.0f}")
            print(f"   Will be slow but should work")
        else:
            print(f"   ✅ Good size: {n_atoms} atoms")
        
        return True
    
    def create_dft_calculator(self, atoms: Atoms, label: str) -> Optional[GPAW]:
        """Create GPAW calculator with appropriate settings"""
        
        if not GPAW_AVAILABLE:
            return None
        
        # Adjust vacuum based on system size
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
    
    def optimize_structure(self, atoms: Atoms, label: str) -> Tuple[Atoms, float, bool]:
        """
        Optimize structure with DFT and comprehensive error handling
        
        Returns:
            (optimized_atoms, energy, converged)
        """
        print(f"\n🔧 DFT optimization: {label}")
        
        # Size check first
        if not self.validate_system_size(atoms, label):
            return atoms.copy(), 0.0, False
        
        # Geometry validation
        validation = self.validator.validate_structure(atoms, label)
        if not validation['valid']:
            print(f"   ⚠️  Attempting to fix geometry issues...")
            atoms, fixed = self.validator.fix_structure(atoms, label)
            if not fixed:
                print(f"   ❌ Could not fix geometry - optimization may fail")
        
        # Create calculator
        calc = self.create_dft_calculator(atoms, label)
        if calc is None:
            return atoms.copy(), 0.0, False
        
        atoms.calc = calc
        
        # Set up optimizer
        opt_class = self.opt_params['optimizer']
        opt = opt_class(atoms, logfile=f'{label}_opt.log')
        
        # Run optimization with monitoring
        start_time = time.time()
        converged = False
        
        try:
            print(f"   Starting optimization...")
            print(f"   Target: fmax < {self.opt_params['fmax']} eV/Å")
            
            # Custom convergence monitoring
            step_count = [0]  # Use list for closure
            
            def monitor_convergence():
                step_count[0] += 1
                if step_count[0] % 10 == 0:
                    try:
                        forces = atoms.get_forces()
                        fmax = np.sqrt((forces**2).sum(axis=1).max())
                        energy = atoms.get_potential_energy()
                        print(f"     Step {step_count[0]}: E={energy:.4f} eV, Fmax={fmax:.4f} eV/Å")
                    except:
                        pass
            
            # Run optimization
            opt.run(fmax=self.opt_params['fmax'], 
                   steps=self.opt_params['max_steps'])
            
            # Check convergence
            forces = atoms.get_forces()
            fmax = np.sqrt((forces**2).sum(axis=1).max())
            converged = fmax < self.opt_params['fmax']
            
            energy = atoms.get_potential_energy()
            elapsed = time.time() - start_time
            
            print(f"   ✅ Optimization complete")
            print(f"   Energy: {energy:.4f} eV")
            print(f"   Max force: {fmax:.4f} eV/Å")
            print(f"   Converged: {converged}")
            print(f"   Time: {elapsed/60:.1f} minutes")
            
            # Save optimized structure
            write(f'{label}_optimized.xyz', atoms)
            
            return atoms.copy(), energy, converged
            
        except Exception as e:
            elapsed = time.time() - start_time
            print(f"   ❌ Optimization failed: {e}")
            print(f"   Time: {elapsed/60:.1f} minutes")
            
            # Try to save partial result
            try:
                energy = atoms.get_potential_energy()
                write(f'{label}_partial.xyz', atoms)
                print(f"   💾 Saved partial result")
                return atoms.copy(), energy, False
            except:
                return atoms.copy(), 0.0, False
    
    def load_small_surface(self) -> Optional[Atoms]:
        """Load appropriately sized surface for DFT"""
        
        print(f"\n📂 Loading small surface for DFT...")
        
        # Priority order: smallest to largest
        surface_files = [
            'cluster_si4o6h6.xyz',               # ~16 atoms (best)
            'cluster_surface_patch_small.xyz',   # ~30 atoms  
            'small_si100_2x2_4L_OH.xyz',        # ~80 atoms
            'small_si100_2x2_6L_OH.xyz',        # ~120 atoms
        ]
        
        for surf_file in surface_files:
            if os.path.exists(surf_file):
                try:
                    surface = read(surf_file)
                    
                    if self.validate_system_size(surface, surf_file):
                        print(f"   ✅ Loaded: {surf_file}")
                        print(f"   Surface atoms: {len(surface)}")
                        return surface
                    else:
                        print(f"   ⚠️  {surf_file} too large, trying next...")
                        
                except Exception as e:
                    print(f"   ❌ Error loading {surf_file}: {e}")
        
        print(f"   ❌ No suitable small surface found!")
        print(f"   Run create_small_surfaces.py or cluster_models.py first")
        return None
    
    def load_tma_molecule(self) -> Optional[Atoms]:
        """Load TMA molecule"""
        
        print(f"\n📂 Loading TMA molecule...")
        
        tma_files = [
            'trimethylaluminum_optimized.xyz',
            'test_tma_dft.xyz'
        ]
        
        for tma_file in tma_files:
            if os.path.exists(tma_file):
                try:
                    tma = read(tma_file)
                    print(f"   ✅ Loaded TMA: {tma_file} ({len(tma)} atoms)")
                    return tma
                except Exception as e:
                    print(f"   ❌ Error loading {tma_file}: {e}")
        
        print(f"   ⚠️  No TMA file found, creating simple structure...")
        return self.create_simple_tma()
    
    def create_simple_tma(self) -> Atoms:
        """Create simple TMA structure"""
        
        positions = [
            [0.0, 0.0, 0.0],  # Al center
        ]
        symbols = ['Al']
        
        # Three methyl groups
        for i, angle in enumerate([0, 2*np.pi/3, 4*np.pi/3]):
            # Carbon
            c_pos = np.array([1.96 * np.cos(angle), 1.96 * np.sin(angle), 0.2])
            positions.append(c_pos)
            symbols.append('C')
            
            # Three hydrogens per carbon
            h_positions = [
                c_pos + np.array([0, 0, 1.09]),
                c_pos + np.array([0.9, 0, -0.5]),
                c_pos + np.array([-0.45, 0.8, -0.5])
            ]
            positions.extend(h_positions)
            symbols.extend(['H', 'H', 'H'])
        
        tma = Atoms(symbols=symbols, positions=positions)
        tma.center(vacuum=8.0)
        
        print(f"   ✅ Created simple TMA: {len(tma)} atoms")
        return tma
    
    def setup_tma_surface_system(self, surface: Atoms, tma: Atoms) -> Tuple[Atoms, bool]:
        """
        Create TMA + surface system with proper separation
        
        Returns:
            (combined_system, is_valid)
        """
        print(f"\n🔗 Setting up TMA + surface system...")
        
        # Check total size first
        total_atoms = len(surface) + len(tma)
        print(f"   Surface: {len(surface)} atoms")
        print(f"   TMA: {len(tma)} atoms") 
        print(f"   Total: {total_atoms} atoms")
        
        if total_atoms > self.MAX_ATOMS:
            print(f"   ❌ Combined system too large: {total_atoms} > {self.MAX_ATOMS}")
            return surface.copy(), False
        
        # Position TMA above surface
        surface_positions = surface.get_positions()
        surface_z_max = surface_positions[:, 2].max()
        
        # Place TMA 4 Å above surface (literature standard)
        tma_copy = tma.copy()
        tma_positions = tma_copy.get_positions()
        tma_center = tma_positions.mean(axis=0)
        
        # Position TMA above surface center
        surface_center = surface_positions.mean(axis=0)
        target_position = surface_center.copy()
        target_position[2] = surface_z_max + 4.0  # 4 Å separation
        
        offset = target_position - tma_center
        tma_copy.translate(offset)
        
        # Combine systems
        combined_positions = np.vstack([
            surface.get_positions(),
            tma_copy.get_positions()
        ])
        
        combined_symbols = (
            list(surface.get_chemical_symbols()) +
            list(tma_copy.get_chemical_symbols())
        )
        
        combined = Atoms(
            symbols=combined_symbols,
            positions=combined_positions,
            cell=surface.get_cell(),
            pbc=surface.get_pbc()
        )
        
        # Copy constraints from surface
        if surface.constraints:
            combined.set_constraint(surface.constraints[0])
        
        # Add more vacuum
        combined.center(vacuum=12.0)
        
        print(f"   ✅ Combined system created")
        print(f"   TMA-surface separation: 4.0 Å")
        
        return combined, True
    
    def run_single_reaction_study(self) -> dict:
        """
        Run single TMA adsorption study with DFT
        
        Returns:
            Dictionary with results
        """
        print(f"\n{'='*70}")
        print("🔬 DFT MLD SINGLE REACTION STUDY")
        print("="*70)
        print("Literature-based approach with small systems")
        
        results = {
            'success': False,
            'surface_atoms': 0,
            'tma_atoms': 0,
            'total_atoms': 0,
            'surface_energy': None,
            'tma_energy': None,
            'combined_energy': None,
            'adsorption_energy': None,
            'computation_time': 0
        }
        
        start_time = time.time()
        
        # 1. Load surface
        surface = self.load_small_surface()
        if surface is None:
            return results
        
        results['surface_atoms'] = len(surface)
        
        # 2. Load TMA
        tma = self.load_tma_molecule()
        if tma is None:
            return results
        
        results['tma_atoms'] = len(tma)
        
        # 3. Setup combined system
        combined, is_valid = self.setup_tma_surface_system(surface, tma)
        if not is_valid:
            return results
        
        results['total_atoms'] = len(combined)
        
        # 4. Optimize surface alone
        print(f"\n📊 Step 1: Surface optimization")
        surface_opt, surface_energy, surf_converged = self.optimize_structure(
            surface, "surface_alone"
        )
        
        if surf_converged:
            results['surface_energy'] = surface_energy
            print(f"   ✅ Surface energy: {surface_energy:.4f} eV")
        else:
            print(f"   ⚠️  Surface optimization had issues")
        
        # 5. Optimize TMA alone
        print(f"\n📊 Step 2: TMA optimization")
        tma_opt, tma_energy, tma_converged = self.optimize_structure(
            tma, "tma_alone"
        )
        
        if tma_converged:
            results['tma_energy'] = tma_energy
            print(f"   ✅ TMA energy: {tma_energy:.4f} eV")
        else:
            print(f"   ⚠️  TMA optimization had issues")
        
        # 6. Optimize combined system
        print(f"\n📊 Step 3: Combined system optimization")
        combined_opt, combined_energy, comb_converged = self.optimize_structure(
            combined, "tma_surface_combined"
        )
        
        if comb_converged:
            results['combined_energy'] = combined_energy
            print(f"   ✅ Combined energy: {combined_energy:.4f} eV")
        else:
            print(f"   ⚠️  Combined optimization had issues")
        
        # 7. Calculate adsorption energy
        if all([results['surface_energy'], results['tma_energy'], results['combined_energy']]):
            adsorption_energy = (results['combined_energy'] - 
                               results['surface_energy'] - 
                               results['tma_energy'])
            
            results['adsorption_energy'] = adsorption_energy
            
            print(f"\n📊 Energetic Analysis:")
            print(f"   Surface energy: {results['surface_energy']:.4f} eV")
            print(f"   TMA energy: {results['tma_energy']:.4f} eV")
            print(f"   Combined energy: {results['combined_energy']:.4f} eV")
            print(f"   Adsorption energy: {adsorption_energy:.4f} eV")
            
            if adsorption_energy < 0:
                print(f"   ✅ Exothermic adsorption (favorable)")
            else:
                print(f"   ⚠️  Endothermic adsorption")
            
            # Compare to literature values (-1.5 to -2.5 eV typical)
            if -3.0 < adsorption_energy < -1.0:
                print(f"   ✅ Reasonable compared to literature")
            else:
                print(f"   ⚠️  Unusual compared to literature (-1.5 to -2.5 eV)")
        
        results['computation_time'] = time.time() - start_time
        results['success'] = True
        
        return results

def main():
    """Run small DFT MLD simulation"""
    
    if not GPAW_AVAILABLE:
        print("❌ GPAW not available! Install with: pip install gpaw")
        return 1
    
    print("🔬 Small DFT MLD Simulation")
    print("Literature-validated approach with appropriate system sizes")
    print("="*70)
    
    # Create simulator with size limits
    simulator = SmallDFTMLDSimulator(max_atoms=200)
    
    # Run single reaction study
    results = simulator.run_single_reaction_study()
    
    # Final summary
    print(f"\n{'='*70}")
    print("🏁 DFT MLD SIMULATION COMPLETE")
    print("="*70)
    
    if results['success']:
        print(f"✅ Simulation successful!")
        print(f"   Total computation time: {results['computation_time']/3600:.2f} hours")
        print(f"   System size: {results['total_atoms']} atoms")
        
        if results['adsorption_energy']:
            print(f"   Adsorption energy: {results['adsorption_energy']:.3f} eV")
        
        print(f"\n📁 Output files:")
        print(f"   - *_optimized.xyz (final structures)")
        print(f"   - *_dft.txt (DFT calculation logs)")
        print(f"   - *_opt.log (optimization logs)")
        
    else:
        print(f"❌ Simulation failed")
        print(f"   Check system sizes and DFT setup")
    
    print(f"\n💡 Next steps:")
    print(f"   • Analyze results with surface_visualizer.py")
    print(f"   • Compare energetics to literature values")
    print(f"   • Scale up to larger systems if needed")
    
    return 0 if results['success'] else 1

if __name__ == "__main__":
    exit(main())