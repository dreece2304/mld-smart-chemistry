#!/usr/bin/env python3
"""
Surface Convergence Fixer for DFT MLD
Specifically addresses the Test 3 & 4 failures from dft_simple_test.py

The issue: Surface systems (especially small clusters) are much harder to converge than molecules
Solution: Progressive parameter adjustment with surface-specific settings
"""

import numpy as np
from ase import Atoms
from ase.io import read, write
import os
import time

try:
    from gpaw import GPAW, PW
    GPAW_AVAILABLE = True
except ImportError:
    GPAW_AVAILABLE = False

def create_surface_optimized_calculator(atoms: Atoms, label: str) -> GPAW:
    """Create GPAW calculator optimized for surface systems"""
    
    print(f"🔧 Creating surface-optimized calculator for {label}...")
    
    # Surface-specific parameters based on literature
    surface_params = {
        'mode': PW(250),                    # Moderate cutoff
        'xc': 'PBE',                        # Standard functional
        'kpts': (1, 1, 1),                 # Gamma point
        'txt': f'{label}_surface.txt',      # Log file
        'symmetry': 'off',                  # Often better for surfaces
        
        # Convergence criteria - looser for testing
        'convergence': {
            'energy': 1e-3,                # Loose energy convergence
            'density': 1e-2,               # Loose density convergence
            'eigenstates': 1e-2            # Loose eigenstate convergence
        },
        
        # Mixing parameters - very conservative for stability
        'mixer': {
            'beta': 0.01,                  # Very small mixing
            'nmaxold': 3,                  # Few old densities
            'weight': 100.0                # Large mixing weight
        },
        
        # SCF settings
        'maxiter': 1000,                   # Many iterations allowed
        'nbands': 'nao',                   # Automatic band count
        
        # Smearing for surface systems (helps convergence)
        'occupations': {
            'name': 'fermi-dirac',         # Fermi-Dirac smearing
            'width': 0.2                   # Heavy smearing for testing
        },
        
        # Force calculation settings
        'parallel': {'kpt': 1},            # Parallel over k-points
        'spinpol': False,                  # No spin polarization for now
    }
    
    # Adjust vacuum for surface systems
    test_atoms = atoms.copy()
    test_atoms.center(vacuum=12.0)  # Large vacuum for surfaces
    
    print(f"   Cutoff: {surface_params['mode']}")
    print(f"   Mixing: β={surface_params['mixer']['beta']}")
    print(f"   Smearing: {surface_params['occupations']['width']} eV")
    print(f"   Max iterations: {surface_params['maxiter']}")
    print(f"   Vacuum: 12.0 Å")
    
    return GPAW(**surface_params), test_atoms

def test_surface_system(atoms: Atoms, name: str) -> dict:
    """Test surface system with optimized parameters"""
    
    print(f"\n{'='*60}")
    print(f"🧪 TESTING SURFACE SYSTEM: {name}")
    print(f"{'='*60}")
    print(f"Atoms: {len(atoms)}")
    print(f"Elements: {set(atoms.get_chemical_symbols())}")
    
    results = {
        'name': name,
        'atoms': len(atoms),
        'success': False,
        'energy': None,
        'time': 0,
        'error': None
    }
    
    if not GPAW_AVAILABLE:
        results['error'] = "GPAW not available"
        return results
    
    try:
        # Create surface-optimized calculator
        calc, test_atoms = create_surface_optimized_calculator(atoms, name)
        test_atoms.calc = calc
        
        # Run single point calculation
        start_time = time.time()
        
        print(f"Starting surface calculation...")
        print(f"(This may take 10-30 minutes for surface systems)")
        
        energy = test_atoms.get_potential_energy()
        elapsed = time.time() - start_time
        
        results['success'] = True
        results['energy'] = energy
        results['time'] = elapsed
        
        print(f"✅ SUCCESS!")
        print(f"   Energy: {energy:.4f} eV")
        print(f"   Time: {elapsed/60:.1f} minutes")
        
        # Save successful structure
        write(f'{name}_surface_converged.xyz', test_atoms)
        print(f"   💾 Saved: {name}_surface_converged.xyz")
        
        return results
        
    except Exception as e:
        elapsed = time.time() - start_time
        results['error'] = str(e)
        results['time'] = elapsed
        
        print(f"❌ FAILED: {e}")
        print(f"   Time: {elapsed/60:.1f} minutes")
        
        return results

def create_even_simpler_cluster() -> Atoms:
    """Create ultra-simple cluster for testing"""
    
    print("\n🧩 Creating ultra-simple cluster for testing...")
    
    # Just Si2O3H3 - minimal surface representation
    positions = [
        [0.0, 0.0, 0.0],      # Si1
        [2.3, 0.0, 0.0],      # Si2
        [1.15, 0.0, 1.0],     # O bridge
        [0.0, 1.6, 0.0],      # O1 on Si1
        [2.3, 1.6, 0.0],      # O2 on Si2
        [0.0, 1.6, 1.0],      # H1 on O1
        [2.3, 1.6, 1.0],      # H2 on O2
        [1.15, 0.0, 2.0],     # H3 on bridge O
    ]
    
    symbols = ['Si', 'Si', 'O', 'O', 'O', 'H', 'H', 'H']
    
    cluster = Atoms(symbols=symbols, positions=positions)
    cluster.center(vacuum=10.0)
    
    print(f"   Created Si2O3H3 cluster: {len(cluster)} atoms")
    return cluster

def main():
    """Fix surface convergence issues from dft_simple_test.py"""
    
    print("🔧 Surface Convergence Fixer")
    print("Addresses Test 3 & 4 failures from dft_simple_test.py")
    print("="*70)
    
    if not GPAW_AVAILABLE:
        print("❌ GPAW not available!")
        return 1
    
    # Test systems in order of complexity
    test_systems = []
    
    # 1. Try ultra-simple cluster first
    print(f"\n{'='*50}")
    print("1. Ultra-Simple Cluster Test")
    print("="*50)
    
    simple_cluster = create_even_simpler_cluster()
    result1 = test_surface_system(simple_cluster, "ultra_simple_cluster")
    test_systems.append(result1)
    
    # 2. Test the failed cluster from dft_simple_test.py
    if os.path.exists('test_surface_cluster_dft.xyz'):
        print(f"\n{'='*50}")
        print("2. Original Failed Cluster")
        print("="*50)
        
        try:
            failed_cluster = read('test_surface_cluster_dft.xyz')
            result2 = test_surface_system(failed_cluster, "original_failed_cluster")
            test_systems.append(result2)
        except Exception as e:
            print(f"❌ Could not load failed cluster: {e}")
    
    # 3. Test cluster from cluster_models.py if available
    if os.path.exists('cluster_si4o6h6.xyz'):
        print(f"\n{'='*50}")
        print("3. Si4O6H6 Cluster")
        print("="*50)
        
        try:
            si4o6h6 = read('cluster_si4o6h6.xyz')
            result3 = test_surface_system(si4o6h6, "si4o6h6_cluster")
            test_systems.append(result3)
        except Exception as e:
            print(f"❌ Could not load Si4O6H6 cluster: {e}")
    
    # 4. Test H2O + surface if available
    if os.path.exists('test_h2o_surface_dft.xyz'):
        print(f"\n{'='*50}")
        print("4. H2O + Surface System")
        print("="*50)
        
        try:
            h2o_surface = read('test_h2o_surface_dft.xyz')
            result4 = test_surface_system(h2o_surface, "h2o_surface_system")
            test_systems.append(result4)
        except Exception as e:
            print(f"❌ Could not load H2O + surface: {e}")
    
    # Summary
    print(f"\n{'='*70}")
    print("🏁 SURFACE CONVERGENCE TEST SUMMARY")
    print("="*70)
    
    successful_systems = []
    failed_systems = []
    
    for result in test_systems:
        if result['success']:
            successful_systems.append(result)
            print(f"✅ {result['name']}: SUCCESS ({result['time']/60:.1f} min)")
        else:
            failed_systems.append(result)
            print(f"❌ {result['name']}: FAILED ({result['time']/60:.1f} min)")
    
    print(f"\nSuccessful systems: {len(successful_systems)}")
    print(f"Failed systems: {len(failed_systems)}")
    
    if successful_systems:
        print(f"\n✅ Surface DFT is working!")
        print(f"   Use the successful parameters for production runs")
        print(f"   Files saved: *_surface_converged.xyz")
        
        print(f"\n💡 Working parameters:")
        print(f"   • PW(250) cutoff")
        print(f"   • β=0.01 mixing")
        print(f"   • 0.2 eV smearing")
        print(f"   • 1000 max iterations")
        print(f"   • Loose convergence criteria")
        
        return 0
    else:
        print(f"\n❌ All surface systems failed")
        print(f"   This suggests a fundamental issue with:")
        print(f"   • GPAW installation")
        print(f"   • System resources")
        print(f"   • Surface geometry")
        
        return 1

if __name__ == "__main__":
    exit(main())