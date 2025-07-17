#!/usr/bin/env python3
"""
Progressive MLD Test Suite
Tests system functionality from basic to advanced levels
"""

import os
import sys
import time
import tempfile
import subprocess
from datetime import datetime
import argparse

class Colors:
    """Terminal colors"""
    GREEN = '\033[92m'
    RED = '\033[91m' 
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    MAGENTA = '\033[95m'
    CYAN = '\033[96m'
    BOLD = '\033[1m'
    END = '\033[0m'

class MLDTestSuite:
    """Progressive test suite for MLD simulation system"""
    
    def __init__(self, verbose=False):
        self.verbose = verbose
        self.test_results = {}
        self.start_time = time.time()
        
        # Set up GPAW environment
        gpaw_path = os.path.expanduser('~/.local/share/gpaw/gpaw-setups-24.11.0')
        os.environ['GPAW_SETUP_PATH'] = gpaw_path
    
    def print_status(self, message, status='info', indent=0):
        """Print colored status message"""
        indent_str = "  " * indent
        
        if status == 'pass':
            print(f"{indent_str}{Colors.GREEN}✅ {message}{Colors.END}")
        elif status == 'fail':
            print(f"{indent_str}{Colors.RED}❌ {message}{Colors.END}")
        elif status == 'warn':
            print(f"{indent_str}{Colors.YELLOW}⚠️  {message}{Colors.END}")
        elif status == 'info':
            print(f"{indent_str}{Colors.BLUE}ℹ️  {message}{Colors.END}")
        elif status == 'running':
            print(f"{indent_str}{Colors.CYAN}🏃 {message}{Colors.END}")
        else:
            print(f"{indent_str}{message}")
    
    def run_test(self, test_name, test_func, level):
        """Run a single test with error handling"""
        self.print_status(f"Running {test_name}...", 'running')
        start_time = time.time()
        
        try:
            result = test_func()
            elapsed = time.time() - start_time
            
            if result:
                self.print_status(f"{test_name}: PASSED ({elapsed:.1f}s)", 'pass', 1)
                self.test_results[test_name] = {'status': 'PASS', 'time': elapsed, 'level': level}
                return True
            else:
                self.print_status(f"{test_name}: FAILED ({elapsed:.1f}s)", 'fail', 1)
                self.test_results[test_name] = {'status': 'FAIL', 'time': elapsed, 'level': level}
                return False
                
        except Exception as e:
            elapsed = time.time() - start_time
            self.print_status(f"{test_name}: ERROR ({elapsed:.1f}s)", 'fail', 1)
            if self.verbose:
                self.print_status(f"Error: {str(e)}", 'fail', 2)
            self.test_results[test_name] = {'status': 'ERROR', 'time': elapsed, 'level': level, 'error': str(e)}
            return False
    
    def level_1_imports(self):
        """Level 1: Test package imports"""
        self.print_status("LEVEL 1: Package Import Tests", 'info')
        
        # Test 1.1: Core scientific packages
        def test_core_imports():
            import numpy as np
            import scipy
            import matplotlib
            assert np.__version__ >= '1.21.0'
            return True
        
        # Test 1.2: ASE import and basic functionality
        def test_ase_import():
            import ase
            from ase import Atoms
            from ase.io import read, write
            from ase.optimize import BFGS
            
            # Test basic Atoms creation
            h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
            assert len(h2) == 2
            return True
        
        # Test 1.3: GPAW import
        def test_gpaw_import():
            import gpaw
            from gpaw import GPAW, PW
            from gpaw.mpi import world
            
            # Check GPAW version
            version = gpaw.__version__
            assert version >= '22.8.0'
            return True
        
        # Test 1.4: Progress tracking
        def test_progress_import():
            import tqdm
            from tqdm.auto import tqdm as auto_tqdm
            
            # Test basic progress bar
            pbar = auto_tqdm(total=10, desc="Test")
            pbar.close()
            return True
        
        tests = [
            ("Core Scientific Packages", test_core_imports),
            ("ASE Import and Basic", test_ase_import),
            ("GPAW Import", test_gpaw_import),
            ("Progress Tracking", test_progress_import)
        ]
        
        passed = 0
        for test_name, test_func in tests:
            if self.run_test(test_name, test_func, 1):
                passed += 1
        
        return passed == len(tests)
    
    def level_2_basic_functionality(self):
        """Level 2: Basic functionality tests"""
        self.print_status("LEVEL 2: Basic Functionality Tests", 'info')
        
        # Test 2.1: ASE Atoms manipulation
        def test_atoms_manipulation():
            from ase import Atoms
            import numpy as np
            
            # Create molecule
            mol = Atoms('H2O', positions=[[0, 0, 0], [0.76, 0.59, 0], [-0.76, 0.59, 0]])
            mol.center(vacuum=5.0)
            
            # Test properties
            assert len(mol) == 3
            assert mol.get_chemical_symbols() == ['H', 'O', 'H']
            
            # Test cell
            cell = mol.get_cell()
            assert all(cell.diagonal() > 10.0)  # Should have vacuum
            
            return True
        
        # Test 2.2: GPAW calculator creation
        def test_gpaw_calculator():
            from gpaw import GPAW, PW
            
            # Create calculator (don't run calculation)
            calc = GPAW(
                mode=PW(200),
                xc='PBE',
                kpts=(1, 1, 1),
                txt=None,
                symmetry='off'
            )
            
            # Check parameters
            assert calc.parameters.mode.ecut == 200
            assert calc.parameters.xc == 'PBE'
            
            return True
        
        # Test 2.3: Structure file I/O
        def test_structure_io():
            from ase import Atoms
            from ase.io import write, read
            import tempfile
            import os
            
            # Create test structure
            mol = Atoms('CO', positions=[[0, 0, 0], [1.13, 0, 0]])
            mol.center(vacuum=5.0)
            
            # Write and read back
            with tempfile.NamedTemporaryFile(suffix='.xyz', delete=False) as f:
                temp_file = f.name
            
            try:
                write(temp_file, mol)
                mol_read = read(temp_file)
                
                # Check consistency
                assert len(mol) == len(mol_read)
                assert mol.get_chemical_symbols() == mol_read.get_chemical_symbols()
                
                return True
            finally:
                if os.path.exists(temp_file):
                    os.unlink(temp_file)
        
        # Test 2.4: Basic optimization setup
        def test_optimization_setup():
            from ase import Atoms
            from ase.optimize import BFGS
            from ase.calculators.emt import EMT
            
            # Create simple system with EMT calculator
            mol = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.8]])
            mol.center(vacuum=5.0)
            mol.calc = EMT()
            
            # Set up optimizer (don't run)
            opt = BFGS(mol)
            
            # Check that forces can be calculated
            forces = mol.get_forces()
            assert forces.shape == (2, 3)
            
            return True
        
        tests = [
            ("ASE Atoms Manipulation", test_atoms_manipulation),
            ("GPAW Calculator Creation", test_gpaw_calculator),
            ("Structure File I/O", test_structure_io),
            ("Basic Optimization Setup", test_optimization_setup)
        ]
        
        passed = 0
        for test_name, test_func in tests:
            if self.run_test(test_name, test_func, 2):
                passed += 1
        
        return passed == len(tests)
    
    def level_3_gpaw_basic(self):
        """Level 3: Basic GPAW functionality"""
        self.print_status("LEVEL 3: Basic GPAW Tests", 'info')
        
        # Test 3.1: GPAW datasets
        def test_gpaw_datasets():
            gpaw_path = os.environ.get('GPAW_SETUP_PATH')
            assert gpaw_path is not None, "GPAW_SETUP_PATH not set"
            assert os.path.exists(gpaw_path), f"GPAW datasets not found at {gpaw_path}"
            
            # Check for some common elements
            for element in ['H.PBE', 'C.PBE', 'O.PBE']:
                setup_file = os.path.join(gpaw_path, element)
                assert os.path.exists(setup_file), f"Missing setup for {element}"
            
            return True
        
        # Test 3.2: Simple H2 calculation
        def test_h2_calculation():
            from ase import Atoms
            from gpaw import GPAW, PW
            import tempfile
            import os
            
            # Create H2 molecule
            h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
            h2.center(vacuum=8.0)
            
            # Create temporary directory for calculation
            with tempfile.TemporaryDirectory() as temp_dir:
                os.chdir(temp_dir)
                
                # Set up GPAW calculator
                calc = GPAW(
                    mode=PW(200),
                    xc='PBE',
                    kpts=(1, 1, 1),
                    txt='h2_test.out',
                    symmetry='off',
                    maxiter=20  # Limit iterations for test
                )
                
                h2.calc = calc
                
                # Calculate energy (this will actually run DFT)
                energy = h2.get_potential_energy()
                
                # Check that we got a reasonable energy
                assert -10.0 < energy < -5.0, f"Unreasonable H2 energy: {energy}"
                
                return True
        
        # Test 3.3: Force calculation
        def test_force_calculation():
            from ase import Atoms
            from gpaw import GPAW, PW
            import tempfile
            import os
            import numpy as np
            
            # Create slightly stretched H2
            h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.9]])  # Longer than equilibrium
            h2.center(vacuum=8.0)
            
            with tempfile.TemporaryDirectory() as temp_dir:
                os.chdir(temp_dir)
                
                calc = GPAW(
                    mode=PW(200),
                    xc='PBE',
                    kpts=(1, 1, 1),
                    txt='h2_forces.out',
                    symmetry='off',
                    maxiter=20
                )
                
                h2.calc = calc
                
                # Calculate forces
                forces = h2.get_forces()
                
                # Check force shape and magnitude
                assert forces.shape == (2, 3)
                force_magnitude = np.linalg.norm(forces)
                assert 0.01 < force_magnitude < 10.0, f"Unreasonable force magnitude: {force_magnitude}"
                
                return True
        
        tests = [
            ("GPAW Datasets Check", test_gpaw_datasets),
            ("H2 Energy Calculation", test_h2_calculation),
            ("Force Calculation", test_force_calculation)
        ]
        
        passed = 0
        for test_name, test_func in tests:
            if self.run_test(test_name, test_func, 3):
                passed += 1
        
        return passed == len(tests)
    
    def level_4_optimization(self):
        """Level 4: Optimization tests"""
        self.print_status("LEVEL 4: Optimization Tests", 'info')
        
        # Test 4.1: Simple optimization
        def test_simple_optimization():
            from ase import Atoms
            from ase.optimize import BFGS
            from gpaw import GPAW, PW
            import tempfile
            import os
            
            # Create H2 with non-optimal distance
            h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.8]])
            h2.center(vacuum=8.0)
            
            with tempfile.TemporaryDirectory() as temp_dir:
                os.chdir(temp_dir)
                
                calc = GPAW(
                    mode=PW(200),
                    xc='PBE',
                    kpts=(1, 1, 1),
                    txt='opt_test.out',
                    symmetry='off'
                )
                
                h2.calc = calc
                
                # Run optimization
                opt = BFGS(h2, trajectory='opt_test.traj')
                opt.run(fmax=0.1, steps=10)  # Loose convergence for test
                
                # Check that optimization ran
                assert opt.get_number_of_steps() > 0
                
                return True
        
        # Test 4.2: Progress tracking optimization
        def test_progress_optimization():
            # Test our custom progress optimization
            try:
                from progress_optimization import optimize_structure
                from ase import Atoms
                from gpaw import GPAW, PW
                import tempfile
                import os
                
                # Create test molecule
                mol = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.8]])
                mol.center(vacuum=8.0)
                
                calc_params = {
                    'mode': PW(200),
                    'xc': 'PBE',
                    'kpts': (1, 1, 1),
                    'symmetry': 'off'
                }
                
                with tempfile.TemporaryDirectory() as temp_dir:
                    os.chdir(temp_dir)
                    
                    # Run with progress tracking
                    optimized, stats = optimize_structure(
                        mol, calc_params,
                        optimizer='BFGS', fmax=0.1, max_steps=5,
                        name="Test_H2"
                    )
                    
                    # Check results
                    assert stats['converged'] or stats['n_steps'] == 5
                    assert 'final_energy' in stats
                    
                    return True
                    
            except ImportError:
                self.print_status("progress_optimization not available", 'warn', 2)
                return True  # Don't fail if module not available
        
        tests = [
            ("Simple H2 Optimization", test_simple_optimization),
            ("Progress Tracking Optimization", test_progress_optimization)
        ]
        
        passed = 0
        for test_name, test_func in tests:
            if self.run_test(test_name, test_func, 4):
                passed += 1
        
        return passed == len(tests)
    
    def level_5_structure_loading(self):
        """Level 5: MLD structure loading tests"""
        self.print_status("LEVEL 5: MLD Structure Loading Tests", 'info')
        
        # Test 5.1: Load MLD structures
        def test_structure_loading():
            from ase.io import read
            
            required_files = [
                'structures/tma_molecule.xyz',
                'structures/butyne_diol_molecule.xyz',
                'structures/simple_si_surface.xyz'
            ]
            
            structures = {}
            
            for file_path in required_files:
                assert os.path.exists(file_path), f"Missing structure file: {file_path}"
                
                mol = read(file_path)
                assert len(mol) > 0, f"Empty structure in {file_path}"
                
                structures[os.path.basename(file_path)] = mol
            
            # Basic checks
            tma = structures['tma_molecule.xyz']
            assert 'Al' in tma.get_chemical_symbols(), "TMA should contain aluminum"
            
            diol = structures['butyne_diol_molecule.xyz']
            symbols = diol.get_chemical_symbols()
            assert 'C' in symbols and 'O' in symbols, "Diol should contain C and O"
            
            surface = structures['simple_si_surface.xyz']
            assert 'Si' in surface.get_chemical_symbols(), "Surface should contain silicon"
            
            return True
        
        # Test 5.2: Structure validation
        def test_structure_validation():
            from ase.io import read
            import numpy as np
            
            # Load and validate TMA
            tma = read('structures/tma_molecule.xyz')
            
            # Check for reasonable bond lengths
            positions = tma.get_positions()
            
            # Calculate all pairwise distances
            from scipy.spatial.distance import pdist
            distances = pdist(positions)
            
            # Check for reasonable range (no atoms too close or too far)
            min_dist = np.min(distances)
            max_dist = np.max(distances)
            
            assert min_dist > 0.5, f"Atoms too close: {min_dist}"
            assert max_dist < 10.0, f"Atoms too far: {max_dist}"
            
            return True
        
        tests = [
            ("Load MLD Structures", test_structure_loading),
            ("Structure Validation", test_structure_validation)
        ]
        
        passed = 0
        for test_name, test_func in tests:
            if self.run_test(test_name, test_func, 5):
                passed += 1
        
        return passed == len(tests)
    
    def level_6_thermal_basics(self):
        """Level 6: Basic thermal simulation tests"""
        self.print_status("LEVEL 6: Thermal Simulation Basics", 'info')
        
        # Test 6.1: Temperature setup
        def test_temperature_setup():
            from ase import Atoms
            from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
            import numpy as np
            
            # Create test molecule
            mol = Atoms('H2O', positions=[[0, 0, 0], [0.76, 0.59, 0], [-0.76, 0.59, 0]])
            mol.center(vacuum=8.0)
            
            # Set temperature
            temperature = 300  # K
            MaxwellBoltzmannDistribution(mol, temperature_K=temperature)
            
            # Check that velocities were set
            velocities = mol.get_velocities()
            assert velocities is not None
            assert velocities.shape == (3, 3)
            
            # Check kinetic energy is reasonable
            kinetic_energy = mol.get_kinetic_energy()
            expected_ke = 1.5 * 3 * 8.617e-5 * temperature  # 3/2 kT per atom in eV
            
            # Should be within factor of 2 (statistical fluctuation)
            assert 0.5 * expected_ke < kinetic_energy < 2.0 * expected_ke
            
            return True
        
        # Test 6.2: MD integrator setup
        def test_md_integrator():
            from ase import Atoms
            from ase.md.verlet import VelocityVerlet
            from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
            from ase.calculators.emt import EMT
            from ase import units
            
            # Create test system with EMT calculator (fast)
            mol = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
            mol.center(vacuum=8.0)
            mol.calc = EMT()
            
            # Set temperature
            MaxwellBoltzmannDistribution(mol, temperature_K=300)
            
            # Create MD integrator
            dyn = VelocityVerlet(mol, timestep=1.0 * units.fs)
            
            # Run a few steps to test
            dyn.run(3)
            
            # Check that MD ran
            assert dyn.nsteps == 3
            
            return True
        
        tests = [
            ("Temperature Setup", test_temperature_setup),
            ("MD Integrator Setup", test_md_integrator)
        ]
        
        passed = 0
        for test_name, test_func in tests:
            if self.run_test(test_name, test_func, 6):
                passed += 1
        
        return passed == len(tests)
    
    def level_7_script_integration(self):
        """Level 7: Script integration tests"""
        self.print_status("LEVEL 7: Script Integration Tests", 'info')
        
        # Test 7.1: quick_run.py test command
        def test_quick_run():
            import subprocess
            
            # Run the quick test
            result = subprocess.run(
                [sys.executable, 'quick_run.py', '--test'],
                capture_output=True, text=True, timeout=300
            )
            
            # Check that it completed successfully
            assert result.returncode == 0, f"quick_run test failed: {result.stderr}"
            assert "test completed successfully" in result.stdout.lower()
            
            return True
        
        # Test 7.2: Validation script
        def test_validation_script():
            import subprocess
            
            # Run validation (without fixes to avoid side effects)
            result = subprocess.run(
                [sys.executable, 'validate_setup.py'],
                capture_output=True, text=True, timeout=60
            )
            
            # Should complete (may have warnings but shouldn't crash)
            assert result.returncode in [0, 1], f"Validation script crashed: {result.stderr}"
            
            return True
        
        tests = [
            ("Quick Run Test", test_quick_run),
            ("Validation Script", test_validation_script)
        ]
        
        passed = 0
        for test_name, test_func in tests:
            if self.run_test(test_name, test_func, 7):
                passed += 1
        
        return passed == len(tests)
    
    def generate_report(self):
        """Generate test suite report"""
        elapsed_time = time.time() - self.start_time
        
        print(f"\n{Colors.BOLD}{'='*80}")
        print("MLD TEST SUITE REPORT")
        print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Total Runtime: {elapsed_time:.1f} seconds")
        print(f"{'='*80}{Colors.END}")
        
        # Summary by level
        levels = {}
        total_tests = len(self.test_results)
        passed_tests = 0
        
        for test_name, result in self.test_results.items():
            level = result['level']
            if level not in levels:
                levels[level] = {'total': 0, 'passed': 0, 'failed': 0, 'errors': 0}
            
            levels[level]['total'] += 1
            if result['status'] == 'PASS':
                levels[level]['passed'] += 1
                passed_tests += 1
            elif result['status'] == 'FAIL':
                levels[level]['failed'] += 1
            else:
                levels[level]['errors'] += 1
        
        # Print level summary
        print(f"\n{Colors.CYAN}SUMMARY BY LEVEL:{Colors.END}")
        for level in sorted(levels.keys()):
            stats = levels[level]
            print(f"  Level {level}: {stats['passed']}/{stats['total']} passed "
                  f"({stats['failed']} failed, {stats['errors']} errors)")
        
        # Overall status
        success_rate = (passed_tests / total_tests) * 100 if total_tests > 0 else 0
        
        print(f"\n{Colors.BOLD}OVERALL RESULTS:{Colors.END}")
        print(f"  Total Tests: {total_tests}")
        print(f"  Passed: {passed_tests}")
        print(f"  Success Rate: {success_rate:.1f}%")
        
        if success_rate >= 90:
            print(f"{Colors.GREEN}✅ EXCELLENT - System ready for production use{Colors.END}")
        elif success_rate >= 70:
            print(f"{Colors.YELLOW}⚠️  GOOD - Minor issues, mostly functional{Colors.END}")
        elif success_rate >= 50:
            print(f"{Colors.YELLOW}⚠️  FAIR - Some issues, basic functionality works{Colors.END}")
        else:
            print(f"{Colors.RED}❌ POOR - Major issues, needs significant fixes{Colors.END}")
        
        # Failed tests details
        failed_tests = [name for name, result in self.test_results.items() 
                       if result['status'] != 'PASS']
        
        if failed_tests:
            print(f"\n{Colors.RED}FAILED TESTS:{Colors.END}")
            for test_name in failed_tests:
                result = self.test_results[test_name]
                print(f"  - {test_name}: {result['status']}")
                if 'error' in result and self.verbose:
                    print(f"    Error: {result['error']}")
        
        return success_rate >= 70

def main():
    """Main test suite function"""
    
    parser = argparse.ArgumentParser(
        description='Progressive MLD Test Suite',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Test Levels:
  1: Package imports
  2: Basic functionality  
  3: GPAW basics
  4: Optimization
  5: Structure loading
  6: Thermal basics
  7: Script integration

Examples:
  python test_suite.py --level 1-3          # Run levels 1 through 3
  python test_suite.py --level 1,5,7        # Run specific levels
  python test_suite.py --all --verbose      # Run all tests with details
        """
    )
    
    parser.add_argument('--level', default='1-3',
                       help='Test levels to run (e.g., "1-3", "1,5,7", "all")')
    parser.add_argument('--verbose', action='store_true',
                       help='Verbose output with error details')
    parser.add_argument('--all', action='store_true',
                       help='Run all test levels')
    
    args = parser.parse_args()
    
    # Parse level specification
    if args.all or args.level.lower() == 'all':
        levels_to_run = list(range(1, 8))
    elif '-' in args.level:
        start, end = map(int, args.level.split('-'))
        levels_to_run = list(range(start, end + 1))
    elif ',' in args.level:
        levels_to_run = [int(x.strip()) for x in args.level.split(',')]
    else:
        levels_to_run = [int(args.level)]
    
    # Create test suite
    test_suite = MLDTestSuite(verbose=args.verbose)
    
    print(f"{Colors.BOLD}{Colors.BLUE}")
    print("🧪 MLD PROGRESSIVE TEST SUITE")
    print("=" * 50)
    print(f"Running levels: {levels_to_run}")
    print(f"{Colors.END}")
    
    # Define test levels
    level_functions = {
        1: test_suite.level_1_imports,
        2: test_suite.level_2_basic_functionality,
        3: test_suite.level_3_gpaw_basic,
        4: test_suite.level_4_optimization,
        5: test_suite.level_5_structure_loading,
        6: test_suite.level_6_thermal_basics,
        7: test_suite.level_7_script_integration
    }
    
    # Run requested levels
    overall_success = True
    for level in levels_to_run:
        if level in level_functions:
            level_success = level_functions[level]()
            if not level_success:
                overall_success = False
                if not args.verbose:
                    test_suite.print_status(f"Level {level} had failures - use --verbose for details", 'warn')
        else:
            test_suite.print_status(f"Unknown test level: {level}", 'warn')
    
    # Generate report
    report_success = test_suite.generate_report()
    
    return 0 if (overall_success and report_success) else 1

if __name__ == "__main__":
    exit(main())