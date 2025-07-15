#!/usr/bin/env python3
"""
Ultra-Simple DFT Test for MLD Setup
Tests basic DFT functionality with minimal systems before full MLD

Test hierarchy:
1. Single H2O molecule (3 atoms)
2. Single TMA molecule (13 atoms) 
3. Small surface cluster (20-30 atoms)
4. H2O approaching surface (50 atoms max)

This ensures DFT setup works before attempting large MLD systems.
"""

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS
from ase.constraints import FixAtoms
import time
import os
from typing import Tuple, Optional

# Import GPAW for DFT
try:
    from gpaw import GPAW, PW
    GPAW_AVAILABLE = True
except ImportError:
    GPAW_AVAILABLE = False
    print("❌ GPAW not available! Install with: pip install gpaw")

# Import progress tracking and smart initial guess
try:
    from add_progress_to_dft import DFTProgress
    from smart_initial_guess import InitialGuessOptimizer, smart_single_point
    SMART_FEATURES = True
except ImportError:
    SMART_FEATURES = False
    print("⚠️  Smart features not available - using basic mode")

class SimpleDFTTest:
    """Ultra-simple DFT tests to verify setup"""
    
    def __init__(self):
        """Initialize with very conservative DFT parameters"""
        
        # Ultra-conservative DFT parameters for testing
        self.dft_params = {
            'mode': PW(200),           # Very low cutoff for speed
            'xc': 'PBE',               # Standard functional
            'kpts': (1, 1, 1),        # Gamma point only
            'txt': None,               # Will be set per calculation
            'symmetry': 'off',         # Often safer
            'convergence': {
                'energy': 1e-3,        # Very loose convergence for testing
                'density': 1e-2,
                'eigenstates': 1e-2
            },
            'mixer': {'beta': 0.01, 'nmaxold': 3},   # Very conservative mixing
            'maxiter': 500,            # More iterations for difficult systems
            'occupations': {'name': 'fermi-dirac', 'width': 0.1},  # Heavy smearing
            'nbands': 'nao'            # Automatic band count
        }
        
        # Size limits for testing
        self.MAX_TEST_ATOMS = 50
        
        # Initialize progress tracking
        if SMART_FEATURES:
            self.progress = DFTProgress()
            self.optimizer = InitialGuessOptimizer()
            print("🧪 Smart DFT Test Suite")
            print("   Testing DFT setup with minimal systems")
            print("   Using smart initial guess and progress tracking")
        else:
            self.progress = None
            self.optimizer = None
            print("🧪 Simple DFT Test Suite")
            print("   Testing DFT setup with minimal systems")
            print("   Conservative parameters for reliability")
        
        if not GPAW_AVAILABLE:
            print("   ❌ GPAW not available - tests will fail")
            return
        else:
            print("   ✅ GPAW available")
    
    def create_dft_calculator(self, label: str) -> Optional[GPAW]:
        """Create ultra-conservative GPAW calculator"""
        
        if not GPAW_AVAILABLE:
            return None
        
        params = self.dft_params.copy()
        params['txt'] = f'{label}_test.txt'
        
        try:
            calc = GPAW(**params)
            return calc
        except Exception as e:
            print(f"❌ Failed to create GPAW calculator: {e}")
            return None
    
    def test_h2o_molecule(self) -> bool:
        """Test 1: Single H2O molecule (3 atoms)"""
        
        print(f"\n{'='*50}")
        print("🧪 TEST 1: H2O Molecule (3 atoms)")
        print("="*50)
        
        # Create H2O molecule with good geometry
        h2o = Atoms('H2O',
                   positions=[
                       [0.0, 0.0, 0.0],      # O
                       [0.757, 0.586, 0.0],  # H1
                       [-0.757, 0.586, 0.0]  # H2
                   ])
        
        # Large vacuum for molecule
        h2o.center(vacuum=8.0)
        
        print(f"Atoms: {len(h2o)}")
        print(f"Cell size: {h2o.cell.diagonal()}")
        
        # Use smart calculation if available
        if SMART_FEATURES:
            if self.progress:
                self.progress.print_status("Starting H2O calculation with smart features...", 'info')
            
            success, energy, elapsed = smart_single_point(h2o, 'h2o_test')
            
            if success:
                if self.progress:
                    self.progress.print_status(f"H2O single point successful! Energy: {energy:.4f} eV", 'success')
                print(f"   Time: {elapsed:.1f} seconds")
                write('test_h2o_dft.xyz', h2o)
                return True
            else:
                if self.progress:
                    self.progress.print_status("H2O single point failed", 'error')
                print(f"   Time: {elapsed:.1f} seconds")
                return False
        
        # Fallback to basic calculation
        calc = self.create_dft_calculator('h2o_test')
        if calc is None:
            return False
        
        h2o.calc = calc
        
        # Test single point calculation
        start_time = time.time()
        
        try:
            energy = h2o.get_potential_energy()
            elapsed = time.time() - start_time
            
            print(f"✅ H2O single point successful!")
            print(f"   Energy: {energy:.4f} eV")
            print(f"   Time: {elapsed:.1f} seconds")
            
            # Save result
            write('test_h2o_dft.xyz', h2o)
            
            return True
            
        except Exception as e:
            elapsed = time.time() - start_time
            print(f"❌ H2O single point failed: {e}")
            print(f"   Time: {elapsed:.1f} seconds")
            return False
    
    def test_tma_molecule(self) -> bool:
        """Test 2: TMA molecule (13 atoms)"""
        
        print(f"\n{'='*50}")
        print("🧪 TEST 2: TMA Molecule (13 atoms)")
        print("="*50)
        
        # Try to load optimized TMA, or create simple one
        try:
            tma = read('trimethylaluminum_optimized.xyz')
            print(f"Loaded TMA from file: {len(tma)} atoms")
        except:
            # Create simple TMA structure
            print("Creating simple TMA structure...")
            tma = self.create_simple_tma()
        
        # Ensure good vacuum
        tma.center(vacuum=8.0)
        
        print(f"Atoms: {len(tma)}")
        print(f"Elements: {set(tma.get_chemical_symbols())}")
        
        # Set up DFT
        calc = self.create_dft_calculator('tma_test')
        if calc is None:
            return False
        
        tma.calc = calc
        
        # Test single point
        start_time = time.time()
        
        try:
            energy = tma.get_potential_energy()
            elapsed = time.time() - start_time
            
            print(f"✅ TMA single point successful!")
            print(f"   Energy: {energy:.4f} eV")
            print(f"   Time: {elapsed:.1f} seconds")
            
            write('test_tma_dft.xyz', tma)
            return True
            
        except Exception as e:
            elapsed = time.time() - start_time
            print(f"❌ TMA single point failed: {e}")
            print(f"   Time: {elapsed:.1f} seconds")
            return False
    
    def create_simple_tma(self) -> Atoms:
        """Create simple TMA structure if file not available"""
        
        # Al at center
        al_pos = np.array([0.0, 0.0, 0.0])
        
        # Three methyl groups in trigonal arrangement
        positions = [al_pos]
        symbols = ['Al']
        
        # Methyl positions (Al-C = 1.96 Å)
        angles = np.array([0, 2*np.pi/3, 4*np.pi/3])
        
        for angle in angles:
            # Carbon position
            c_pos = al_pos + 1.96 * np.array([np.cos(angle), np.sin(angle), 0.2])
            positions.append(c_pos)
            symbols.append('C')
            
            # Three hydrogens per carbon (simplified - tetrahedral)
            h_positions = [
                c_pos + np.array([0, 0, 1.09]),
                c_pos + np.array([0.9, 0, -0.5]),
                c_pos + np.array([-0.9, 0, -0.5])
            ]
            
            positions.extend(h_positions)
            symbols.extend(['H', 'H', 'H'])
        
        return Atoms(symbols=symbols, positions=positions)
    
    def test_small_surface(self) -> bool:
        """Test 3: Small surface cluster (~30 atoms)"""
        
        print(f"\n{'='*50}")
        print("🧪 TEST 3: Small Surface Cluster (~30 atoms)")
        print("="*50)
        
        # Create minimal Si4O6H6 cluster (typical in literature)
        surface_cluster = self.create_si4o6h6_cluster()
        
        print(f"Atoms: {len(surface_cluster)}")
        print(f"Elements: {set(surface_cluster.get_chemical_symbols())}")
        
        if len(surface_cluster) > self.MAX_TEST_ATOMS:
            print(f"❌ Cluster too large: {len(surface_cluster)} > {self.MAX_TEST_ATOMS}")
            return False
        
        # Set up DFT
        calc = self.create_dft_calculator('surface_cluster_test')
        if calc is None:
            return False
        
        surface_cluster.calc = calc
        
        # Test single point
        start_time = time.time()
        
        try:
            energy = surface_cluster.get_potential_energy()
            elapsed = time.time() - start_time
            
            print(f"✅ Surface cluster single point successful!")
            print(f"   Energy: {energy:.4f} eV")
            print(f"   Time: {elapsed:.1f} seconds")
            
            write('test_surface_cluster_dft.xyz', surface_cluster)
            return True
            
        except Exception as e:
            elapsed = time.time() - start_time
            print(f"❌ Surface cluster single point failed: {e}")
            print(f"   Time: {elapsed:.1f} seconds")
            return False
    
    def create_si4o6h6_cluster(self) -> Atoms:
        """Create Si4O6H6 cluster - standard in literature"""
        
        # Simple tetrahedral Si4 cluster with bridging oxygens
        # Based on Halls & Raghavachari (2004) approach
        
        positions = []
        symbols = []
        
        # Central Si4 tetrahedron (simplified)
        si_positions = [
            [0.0, 0.0, 0.0],
            [2.3, 0.0, 0.0],
            [1.15, 2.0, 0.0],
            [1.15, 0.67, 1.6]
        ]
        
        positions.extend(si_positions)
        symbols.extend(['Si'] * 4)
        
        # Bridging oxygens between Si atoms
        o_positions = [
            [1.15, 0.0, 0.0],    # Between Si1-Si2
            [0.58, 1.0, 0.0],    # Between Si1-Si3
            [1.73, 1.0, 0.0],    # Between Si2-Si3
            [0.58, 0.33, 0.8],   # Between Si1-Si4
            [1.73, 0.33, 0.8],   # Between Si2-Si4
            [1.15, 1.33, 0.8]    # Between Si3-Si4
        ]
        
        positions.extend(o_positions)
        symbols.extend(['O'] * 6)
        
        # Terminal H atoms on oxygens
        h_positions = [
            [1.15, 0.0, 1.0],    # On O1
            [0.0, 1.0, 0.0],     # On O2
            [2.3, 1.0, 0.0],     # On O3
            [0.0, 0.33, 0.8],    # On O4
            [2.3, 0.33, 0.8],    # On O5
            [1.15, 2.0, 0.8]     # On O6
        ]
        
        positions.extend(h_positions)
        symbols.extend(['H'] * 6)
        
        cluster = Atoms(symbols=symbols, positions=positions)
        cluster.center(vacuum=10.0)
        
        return cluster
    
    def test_h2o_surface_approach(self) -> bool:
        """Test 4: H2O approaching surface (~40 atoms)"""
        
        print(f"\n{'='*50}")
        print("🧪 TEST 4: H2O + Surface Approach (~40 atoms)")
        print("="*50)
        
        # Create small surface
        surface = self.create_si4o6h6_cluster()
        
        # Create H2O
        h2o = Atoms('H2O',
                   positions=[
                       [0.0, 0.0, 5.0],      # O (5 Å above surface)
                       [0.757, 0.586, 5.0], # H1
                       [-0.757, 0.586, 5.0] # H2
                   ])
        
        # Combine systems
        combined_positions = np.vstack([
            surface.get_positions(),
            h2o.get_positions()
        ])
        
        combined_symbols = (
            list(surface.get_chemical_symbols()) +
            list(h2o.get_chemical_symbols())
        )
        
        combined = Atoms(
            symbols=combined_symbols,
            positions=combined_positions
        )
        
        combined.center(vacuum=8.0)
        
        print(f"Atoms: {len(combined)}")
        print(f"Elements: {set(combined.get_chemical_symbols())}")
        
        if len(combined) > self.MAX_TEST_ATOMS:
            print(f"❌ System too large: {len(combined)} > {self.MAX_TEST_ATOMS}")
            return False
        
        # Set up DFT
        calc = self.create_dft_calculator('h2o_surface_test')
        if calc is None:
            return False
        
        combined.calc = calc
        
        # Test single point
        start_time = time.time()
        
        try:
            energy = combined.get_potential_energy()
            elapsed = time.time() - start_time
            
            print(f"✅ H2O + surface single point successful!")
            print(f"   Energy: {energy:.4f} eV")
            print(f"   Time: {elapsed:.1f} seconds")
            
            write('test_h2o_surface_dft.xyz', combined)
            return True
            
        except Exception as e:
            elapsed = time.time() - start_time
            print(f"❌ H2O + surface single point failed: {e}")
            print(f"   Time: {elapsed:.1f} seconds")
            return False
    
    def run_all_tests(self) -> dict:
        """Run all DFT tests in order"""
        
        print(f"\n{'='*60}")
        print("🧪 RUNNING ALL DFT TESTS")
        print("="*60)
        print("Testing DFT setup with progressively complex systems")
        
        if not GPAW_AVAILABLE:
            print("❌ GPAW not available - cannot run tests")
            return {'all_passed': False, 'error': 'GPAW not available'}
        
        tests = [
            ('H2O molecule', self.test_h2o_molecule),
            ('TMA molecule', self.test_tma_molecule),
            ('Surface cluster', self.test_small_surface),
            ('H2O + surface', self.test_h2o_surface_approach)
        ]
        
        results = {}
        all_passed = True
        
        for test_name, test_func in tests:
            print(f"\n{'─'*40}")
            print(f"Running: {test_name}")
            print("─"*40)
            
            try:
                success = test_func()
                results[test_name] = success
                
                if success:
                    print(f"✅ {test_name}: PASSED")
                else:
                    print(f"❌ {test_name}: FAILED")
                    all_passed = False
                    
            except Exception as e:
                print(f"❌ {test_name}: ERROR - {e}")
                results[test_name] = False
                all_passed = False
        
        # Summary
        print(f"\n{'='*60}")
        print("🏁 TEST SUMMARY")
        print("="*60)
        
        for test_name, passed in results.items():
            status = "✅ PASS" if passed else "❌ FAIL"
            print(f"  {test_name}: {status}")
        
        if all_passed:
            print(f"\n🎉 ALL TESTS PASSED!")
            print(f"   DFT setup is working correctly")
            print(f"   Ready for larger MLD calculations")
        else:
            print(f"\n⚠️  SOME TESTS FAILED")
            print(f"   Check GPAW installation and parameters")
            print(f"   Do not proceed to large MLD calculations yet")
        
        results['all_passed'] = all_passed
        return results

def main():
    """Run simple DFT tests"""
    
    print("🧪 Simple DFT Test Suite for MLD")
    print("Testing basic DFT functionality before full simulations")
    
    tester = SimpleDFTTest()
    results = tester.run_all_tests()
    
    if results['all_passed']:
        print(f"\n💡 Next steps:")
        print(f"  1. DFT is working - proceed to create_small_surfaces.py")
        print(f"  2. Run geometry_validator.py on small surfaces")
        print(f"  3. Use dft_small_mld.py for actual MLD calculations")
        return 0
    else:
        print(f"\n💡 Troubleshooting:")
        print(f"  1. Check GPAW installation: pip install gpaw")
        print(f"  2. Verify ASE compatibility")
        print(f"  3. Check computational resources")
        return 1

if __name__ == "__main__":
    exit(main())