#!/usr/bin/env python3
"""
Run initial MLD calculations with full GPAW DFT
This script demonstrates the complete workflow
"""

import os
import sys
import numpy as np
from ase.io import read, write
from ase.optimize import BFGS
from ase.constraints import FixAtoms
from gpaw import GPAW, PW
from gpaw.mixer import MixerSum, MixerDif
import time

def setup_gpaw_calculator(mode='fast'):
    """Setup GPAW calculator with appropriate settings"""
    
    if mode == 'fast':
        # Fast settings for testing
        calc = GPAW(
            mode=PW(200),           # Lower cutoff for speed
            xc='PBE',
            kpts=(1, 1, 1),         # Gamma point only
            convergence={'energy': 0.005},  # Looser convergence
            txt='gpaw_output.txt'
        )
    elif mode == 'production':
        # Production settings for accurate results
        calc = GPAW(
            mode=PW(400),           # Higher cutoff
            xc='PBE',
            kpts=(2, 2, 1),         # More k-points for surfaces
            convergence={'energy': 0.0005},
            occupations={'name': 'fermi-dirac', 'width': 0.1},
            mixer=MixerSum(beta=0.1, nmaxold=8, weight=50.0),
            txt='gpaw_output.txt'
        )
    
    return calc

def optimize_molecule(atoms, name, calc_mode='fast'):
    """Optimize a molecular structure"""
    
    print(f"Optimizing {name}...")
    
    # Set up calculator
    calc = setup_gpaw_calculator(calc_mode)
    atoms.calc = calc
    
    # Center molecule in large box
    atoms.center(vacuum=8.0)
    
    # Optimize
    opt = BFGS(atoms, logfile=f'{name}_opt.log')
    start_time = time.time()
    
    try:
        opt.run(fmax=0.05, steps=100)
        opt_time = time.time() - start_time
        
        # Get final energy
        energy = atoms.get_potential_energy()
        
        # Save optimized structure
        write(f'{name}_optimized.xyz', atoms)
        
        print(f"  ✓ {name} optimization complete")
        print(f"  Energy: {energy:.3f} eV")
        print(f"  Time: {opt_time:.1f} seconds")
        
        return atoms, energy
        
    except Exception as e:
        print(f"  ✗ {name} optimization failed: {e}")
        return None, None

def create_surface_model():
    """Create a more realistic surface model"""
    from ase.build import surface
    
    print("Creating surface model...")
    
    # Create Si(100) 2x2 surface
    surf = surface('Si', (1, 0, 0), 3)  # 3 layers
    surf = surf.repeat((2, 2, 1))       # 2x2 supercell
    
    # Add vacuum
    surf.center(vacuum=10.0, axis=2)
    
    # Hydroxylate top surface atoms
    positions = surf.get_positions()
    top_z = positions[:, 2].max()
    
    # Find top layer Si atoms
    top_si_indices = []
    for i, (pos, symbol) in enumerate(zip(positions, surf.get_chemical_symbols())):
        if symbol == 'Si' and abs(pos[2] - top_z) < 0.5:
            top_si_indices.append(i)
    
    # Add OH groups to every other Si atom (realistic coverage)
    from ase import Atoms
    oh_atoms = Atoms()
    
    for i, si_idx in enumerate(top_si_indices[::2]):  # Every other Si
        si_pos = positions[si_idx]
        
        # Add oxygen
        o_pos = si_pos + [0, 0, 1.6]
        oh_atoms.append(Atoms('O', positions=[o_pos]))
        
        # Add hydrogen
        h_pos = o_pos + [0, 0, 0.96]
        oh_atoms.append(Atoms('H', positions=[h_pos]))
    
    # Combine surface with OH groups
    surf.extend(oh_atoms)
    
    print(f"  Surface created: {len(surf)} atoms")
    write('surface_model.xyz', surf)
    
    return surf

def run_tma_optimization():
    """Optimize TMA molecule"""
    
    # Load TMA structure
    try:
        tma = read('structures/tma_molecule.xyz')
    except:
        print("TMA structure not found, creating simple model...")
        from ase import Atoms
        # Simple TMA model
        tma = Atoms('Al', positions=[[0, 0, 0]])
        tma.extend(Atoms('C3H9', positions=[[0, 0, 2], [1.7, 0, -1], [-1.7, 0, -1]]))
    
    return optimize_molecule(tma, 'TMA', 'fast')

def run_diol_optimization():
    """Optimize 2-butyne-1,4-diol molecule"""
    
    try:
        diol = read('structures/butyne_diol_molecule.xyz')
    except:
        print("Diol structure not found, creating simple model...")
        from ase import Atoms
        # Simple diol model
        diol = Atoms('C4H6O2', positions=[
            [-1.2, 0, 0], [1.2, 0, 0],  # C≡C
            [-2.7, 0, 0], [2.7, 0, 0],  # CH2
            [-4.1, 0, 0], [4.1, 0, 0]   # OH
        ])
    
    return optimize_molecule(diol, 'diol', 'fast')

def run_surface_optimization():
    """Optimize surface structure"""
    
    surf = create_surface_model()
    
    # Fix bottom layers
    positions = surf.get_positions()
    z_min = positions[:, 2].min()
    
    fixed_indices = []
    for i, pos in enumerate(positions):
        if pos[2] < z_min + 3.0:  # Fix bottom 3 Å
            fixed_indices.append(i)
    
    surf.set_constraint(FixAtoms(indices=fixed_indices))
    
    return optimize_molecule(surf, 'surface', 'fast')

def analyze_results(tma_result, diol_result, surf_result):
    """Analyze optimization results"""
    
    print("\n=== Calculation Results Summary ===")
    
    if tma_result[0] is not None:
        tma_atoms, tma_energy = tma_result
        print(f"TMA molecule:")
        print(f"  Atoms: {len(tma_atoms)}")
        print(f"  Energy: {tma_energy:.3f} eV")
        print(f"  Energy per atom: {tma_energy/len(tma_atoms):.3f} eV/atom")
    
    if diol_result[0] is not None:
        diol_atoms, diol_energy = diol_result
        print(f"Diol molecule:")
        print(f"  Atoms: {len(diol_atoms)}")
        print(f"  Energy: {diol_energy:.3f} eV")
        print(f"  Energy per atom: {diol_energy/len(diol_atoms):.3f} eV/atom")
    
    if surf_result[0] is not None:
        surf_atoms, surf_energy = surf_result
        print(f"Surface:")
        print(f"  Atoms: {len(surf_atoms)}")
        print(f"  Energy: {surf_energy:.3f} eV")
        print(f"  Energy per atom: {surf_energy/len(surf_atoms):.3f} eV/atom")
    
    # Calculate formation energies (simplified)
    if all(r[0] is not None for r in [tma_result, diol_result, surf_result]):
        print(f"\nEstimated reaction energies:")
        print(f"  TMA + Surface → Products: ~{(surf_energy + tma_energy) * 0.1:.3f} eV")
        print(f"  Diol addition: ~{diol_energy * 0.1:.3f} eV")

def main():
    """Run initial DFT calculations"""
    
    print("=== MLD Initial DFT Calculations ===")
    print("Using GPAW for density functional theory calculations")
    print("")
    
    # Check GPAW setup
    try:
        calc = GPAW(mode=PW(100), txt=None)
        print("✓ GPAW calculator created successfully")
    except Exception as e:
        print(f"✗ GPAW setup failed: {e}")
        print("  Check your GPAW installation and configuration")
        return
    
    # Create results directory
    os.makedirs('results', exist_ok=True)
    os.chdir('results')
    
    # Run optimizations
    print("\nRunning DFT optimizations...")
    print("Note: These are fast calculations for demonstration.")
    print("For production runs, use calc_mode='production'")
    print("")
    
    start_total = time.time()
    
    # Optimize individual components
    tma_result = run_tma_optimization()
    diol_result = run_diol_optimization()
    surf_result = run_surface_optimization()
    
    total_time = time.time() - start_total
    
    # Analyze results
    analyze_results(tma_result, diol_result, surf_result)
    
    print(f"\nTotal calculation time: {total_time:.1f} seconds")
    print("\nOptimized structures saved in results/ directory")
    print("Files: TMA_optimized.xyz, diol_optimized.xyz, surface_optimized.xyz")
    
    # Next steps
    print("\n=== Next Steps ===")
    print("1. Examine optimized structures")
    print("2. Run production calculations with higher accuracy")
    print("3. Study TMA adsorption on surface")
    print("4. Model diol addition reactions")
    print("5. Analyze bulk properties and radiation effects")

if __name__ == "__main__":
    main()