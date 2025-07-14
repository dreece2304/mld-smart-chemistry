#!/usr/bin/env python3
"""
Simple DFT test with GPAW
Tests basic functionality before running full calculations
"""

import os
import time
from ase.io import read, write
from ase.optimize import BFGS
from gpaw import GPAW, PW

def test_gpaw_basic():
    """Test basic GPAW functionality"""
    
    print("=== GPAW Basic Test ===")
    
    # Create simple molecule for testing
    from ase import Atoms
    h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
    h2.center(vacuum=5.0)
    
    # Set up simple calculator
    calc = GPAW(
        mode=PW(200),
        xc='PBE',
        txt='h2_test.txt'
    )
    
    h2.calc = calc
    
    print("Testing H2 molecule...")
    start_time = time.time()
    
    try:
        energy = h2.get_potential_energy()
        forces = h2.get_forces()
        
        calc_time = time.time() - start_time
        
        print(f"✓ GPAW calculation successful!")
        print(f"  Energy: {energy:.3f} eV")
        print(f"  Max force: {abs(forces).max():.3f} eV/Å")
        print(f"  Time: {calc_time:.1f} seconds")
        
        return True
        
    except Exception as e:
        print(f"✗ GPAW calculation failed: {e}")
        return False

def test_molecule_optimization():
    """Test optimization of our TMA molecule"""
    
    print("\n=== TMA Molecule Optimization Test ===")
    
    # Load TMA
    try:
        tma = read('structures/tma_molecule.xyz')
        print(f"✓ TMA loaded: {len(tma)} atoms")
    except Exception as e:
        print(f"✗ Could not load TMA: {e}")
        return False
    
    # Center in box
    tma.center(vacuum=6.0)
    
    # Set up calculator (very basic settings for speed)
    calc = GPAW(
        mode=PW(150),  # Low cutoff for speed
        xc='PBE',
        convergence={'energy': 0.01},  # Loose convergence
        maxiter=50,    # Limit iterations
        txt='tma_test.txt'
    )
    
    tma.calc = calc
    
    print("Running optimization...")
    start_time = time.time()
    
    try:
        # Quick optimization
        opt = BFGS(tma, logfile='tma_opt_test.log')
        opt.run(fmax=0.1, steps=20)  # Loose convergence, few steps
        
        calc_time = time.time() - start_time
        energy = tma.get_potential_energy()
        
        print(f"✓ TMA optimization successful!")
        print(f"  Final energy: {energy:.3f} eV")
        print(f"  Energy per atom: {energy/len(tma):.3f} eV/atom")
        print(f"  Time: {calc_time:.1f} seconds")
        
        # Save optimized structure
        write('tma_optimized_test.xyz', tma)
        print("  Optimized structure saved as tma_optimized_test.xyz")
        
        return True
        
    except Exception as e:
        print(f"✗ TMA optimization failed: {e}")
        return False

def test_surface_calculation():
    """Test surface calculation"""
    
    print("\n=== Surface Calculation Test ===")
    
    # Create minimal surface
    from ase.build import surface
    surf = surface('Si', (1, 0, 0), 2)  # Just 2 layers
    surf.center(vacuum=8.0, axis=2)
    
    print(f"Created surface: {len(surf)} atoms")
    
    # Fix bottom layer
    from ase.constraints import FixAtoms
    fixed_indices = [i for i in range(len(surf)//2)]
    surf.set_constraint(FixAtoms(indices=fixed_indices))
    
    # Set up calculator
    calc = GPAW(
        mode=PW(150),
        xc='PBE',
        kpts=(2, 2, 1),  # Some k-points for surface
        convergence={'energy': 0.01},
        maxiter=30,
        txt='surface_test.txt'
    )
    
    surf.calc = calc
    
    print("Running surface calculation...")
    start_time = time.time()
    
    try:
        energy = surf.get_potential_energy()
        calc_time = time.time() - start_time
        
        print(f"✓ Surface calculation successful!")
        print(f"  Energy: {energy:.3f} eV")
        print(f"  Energy per atom: {energy/len(surf):.3f} eV/atom")
        print(f"  Time: {calc_time:.1f} seconds")
        
        write('surface_test.xyz', surf)
        print("  Surface structure saved as surface_test.xyz")
        
        return True
        
    except Exception as e:
        print(f"✗ Surface calculation failed: {e}")
        return False

def main():
    """Run all tests"""
    
    print("DFT Functionality Test Suite")
    print("=" * 40)
    
    # Create test directory
    os.makedirs('test_results', exist_ok=True)
    os.chdir('test_results')
    
    results = []
    
    # Test 1: Basic GPAW
    results.append(test_gpaw_basic())
    
    # Test 2: Molecule optimization
    if results[-1]:  # Only if basic test passed
        results.append(test_molecule_optimization())
    else:
        print("Skipping molecule test due to basic test failure")
        results.append(False)
    
    # Test 3: Surface calculation
    if any(results):  # If any test passed
        results.append(test_surface_calculation())
    else:
        print("Skipping surface test due to previous failures")
        results.append(False)
    
    # Summary
    print("\n" + "=" * 40)
    print("TEST SUMMARY")
    print("=" * 40)
    print(f"Basic GPAW test:     {'PASS' if results[0] else 'FAIL'}")
    print(f"Molecule optimization: {'PASS' if results[1] else 'FAIL'}")
    print(f"Surface calculation:   {'PASS' if results[2] else 'FAIL'}")
    
    if all(results):
        print("\n✓ All tests passed! Ready for full MLD calculations.")
        print("\nNext steps:")
        print("1. Run full optimization: python calculations/mld_deposition_workflow.py")
        print("2. Analyze structures: python analysis/bulk_structure_analysis.py")
        print("3. Study radiation effects: python calculations/radiation/uv_ebeam_exposure.py")
    elif any(results):
        print("\n⚠ Some tests passed. Basic functionality works.")
        print("Check failed tests and adjust settings if needed.")
    else:
        print("\n✗ All tests failed. Check GPAW installation and configuration.")
    
    print(f"\nTest files saved in: {os.getcwd()}")

if __name__ == "__main__":
    main()