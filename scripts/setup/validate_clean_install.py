#!/usr/bin/env python3
"""
Comprehensive validation for clean MLD installation
Tests all components systematically with detailed reporting
"""

import os
import sys
import time
import tempfile
import subprocess
import platform
from datetime import datetime
from pathlib import Path
import argparse

class Colors:
    """Terminal colors for output"""
    GREEN = '\033[92m'
    RED = '\033[91m' 
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    MAGENTA = '\033[95m'
    CYAN = '\033[96m'
    WHITE = '\033[97m'
    BOLD = '\033[1m'
    END = '\033[0m'

class CleanInstallValidator:
    """Comprehensive validator for clean MLD installation"""
    
    def __init__(self, verbose=False):
        self.verbose = verbose
        self.results = {}
        self.start_time = time.time()
        
        # Environment info
        self.env_info = {
            'platform': platform.system(),
            'python_version': platform.python_version(),
            'conda_env': os.environ.get('CONDA_DEFAULT_ENV', 'none'),
            'is_wsl': 'Microsoft' in platform.uname().release,
            'gpaw_setup_path': os.environ.get('GPAW_SETUP_PATH', ''),
            'mpl_backend': os.environ.get('MPLBACKEND', 'default')
        }
    
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
        elif status == 'step':
            print(f"{indent_str}{Colors.CYAN}🔧 {message}{Colors.END}")
        else:
            print(f"{indent_str}{message}")
    
    def test_environment_setup(self):
        """Test basic environment setup"""
        self.print_status("Testing Environment Setup", 'info')
        
        tests = []
        
        # Python version
        py_version = tuple(map(int, platform.python_version().split('.')))
        if py_version >= (3, 8) and py_version <= (3, 11):
            self.print_status(f"Python {platform.python_version()}", 'pass', 1)
            tests.append(True)
        else:
            self.print_status(f"Python {platform.python_version()} (not optimal)", 'warn', 1)
            tests.append(False)
        
        # Conda environment
        conda_env = self.env_info['conda_env']
        if conda_env == 'mld_modeling':
            self.print_status(f"Conda environment: {conda_env}", 'pass', 1)
            tests.append(True)
        else:
            self.print_status(f"Conda environment: {conda_env} (should be mld_modeling)", 'warn', 1)
            tests.append(False)
        
        # WSL detection
        if self.env_info['is_wsl']:
            self.print_status("Running in WSL2", 'info', 1)
            # Check WSL-specific settings
            mpl_backend = self.env_info['mpl_backend']
            if mpl_backend == 'Agg':
                self.print_status("Matplotlib backend: headless (good for WSL)", 'pass', 1)
                tests.append(True)
            else:
                self.print_status(f"Matplotlib backend: {mpl_backend} (may cause issues)", 'warn', 1)
                tests.append(False)
        else:
            self.print_status("Running on native Linux", 'info', 1)
            tests.append(True)
        
        self.results['environment'] = all(tests)
        return all(tests)
    
    def test_core_packages(self):
        """Test core scientific packages"""
        self.print_status("Testing Core Packages", 'info')
        
        # Core packages with minimum versions
        core_packages = {
            'numpy': '1.21.0',
            'scipy': '1.7.0', 
            'matplotlib': '3.5.0',
            'ase': '3.22.0',
            'gpaw': '22.8.0',
            'mpi4py': '3.1.0',
            'tqdm': '4.0.0'
        }
        
        passed_tests = []
        
        for package, min_version in core_packages.items():
            try:
                module = __import__(package)
                version = getattr(module, '__version__', 'unknown')
                self.print_status(f"{package}: {version}", 'pass', 1)
                passed_tests.append(True)
            except ImportError as e:
                self.print_status(f"{package}: not installed ({e})", 'fail', 1)
                passed_tests.append(False)
        
        self.results['core_packages'] = all(passed_tests)
        return all(passed_tests)
    
    def test_gpaw_setup(self):
        """Test GPAW installation and datasets"""
        self.print_status("Testing GPAW Setup", 'info')
        
        tests = []
        
        # GPAW import
        try:
            import gpaw
            self.print_status(f"GPAW import: {gpaw.__version__}", 'pass', 1)
            tests.append(True)
        except ImportError as e:
            self.print_status(f"GPAW import failed: {e}", 'fail', 1)
            tests.append(False)
            self.results['gpaw_setup'] = False
            return False
        
        # GPAW_SETUP_PATH
        gpaw_path = self.env_info['gpaw_setup_path']
        if gpaw_path and os.path.exists(gpaw_path):
            self.print_status(f"GPAW_SETUP_PATH: {gpaw_path}", 'pass', 1)
            tests.append(True)
            
            # Count dataset files (GPAW 2024 uses .gz compressed files)
            setup_files = list(Path(gpaw_path).glob('*.gz'))
            self.print_status(f"GPAW datasets: {len(setup_files)} files", 'pass', 1)
            
            # Check for key elements (look for .gz files)
            key_elements = ['H.PBE.gz', 'C.PBE.gz', 'O.PBE.gz', 'Al.PBE.gz', 'Si.PBE.gz']
            missing_elements = []
            
            for element in key_elements:
                element_file = Path(gpaw_path) / element
                if not element_file.exists():
                    missing_elements.append(element)
            
            if missing_elements:
                self.print_status(f"Missing datasets: {missing_elements}", 'warn', 1)
                tests.append(False)
            else:
                self.print_status("All key element datasets found", 'pass', 1)
                tests.append(True)
                
        else:
            self.print_status(f"GPAW_SETUP_PATH not set or invalid: {gpaw_path}", 'fail', 1)
            tests.append(False)
        
        self.results['gpaw_setup'] = all(tests)
        return all(tests)
    
    def test_mpi_functionality(self):
        """Test MPI parallel functionality"""
        self.print_status("Testing MPI Functionality", 'info')
        
        tests = []
        
        # Basic MPI import
        try:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            size = comm.Get_size()
            rank = comm.Get_rank()
            
            self.print_status(f"MPI size: {size}, rank: {rank}", 'pass', 1)
            tests.append(True)
        except ImportError as e:
            self.print_status(f"MPI import failed: {e}", 'fail', 1)
            tests.append(False)
            self.results['mpi'] = False
            return False
        
        # Test MPI with subprocess (simple test)
        try:
            result = subprocess.run([
                'mpirun', '-np', '2', 'python', '-c',
                'from mpi4py import MPI; print("Rank {}/{}".format(MPI.COMM_WORLD.Get_rank(), MPI.COMM_WORLD.Get_size()))'
            ], capture_output=True, text=True, timeout=10)
            
            if result.returncode == 0:
                lines = result.stdout.strip().split('\n')
                if len(lines) == 2:  # Should have 2 processes
                    self.print_status("MPI multiprocess test: passed", 'pass', 1)
                    tests.append(True)
                else:
                    self.print_status(f"MPI multiprocess test: unexpected output", 'warn', 1)
                    tests.append(False)
            else:
                self.print_status(f"MPI multiprocess test failed: {result.stderr}", 'warn', 1)
                tests.append(False)
                
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            self.print_status(f"MPI multiprocess test error: {e}", 'warn', 1)
            tests.append(False)
        
        self.results['mpi'] = all(tests)
        return all(tests)
    
    def test_visualization(self):
        """Test visualization capabilities (headless)"""
        self.print_status("Testing Visualization (Headless)", 'info')
        
        tests = []
        
        # Matplotlib headless test
        try:
            import matplotlib
            matplotlib.use('Agg')  # Force headless
            import matplotlib.pyplot as plt
            import numpy as np
            
            # Create simple test plot
            x = np.linspace(0, 10, 100)
            y = np.sin(x)
            
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.plot(x, y, 'b-', linewidth=2)
            ax.set_title('Validation Test Plot')
            ax.set_xlabel('X')
            ax.set_ylabel('sin(X)')
            ax.grid(True, alpha=0.3)
            
            # Save to temporary file
            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                fig.savefig(tmp.name, dpi=100, bbox_inches='tight')
                plt.close(fig)
                
                # Check file was created
                if os.path.exists(tmp.name) and os.path.getsize(tmp.name) > 0:
                    self.print_status("Matplotlib plotting: working", 'pass', 1)
                    tests.append(True)
                    os.unlink(tmp.name)  # Cleanup
                else:
                    self.print_status("Matplotlib plotting: file not created", 'fail', 1)
                    tests.append(False)
                    
        except Exception as e:
            self.print_status(f"Matplotlib test failed: {e}", 'fail', 1)
            tests.append(False)
        
        # Optional visualization packages (don't affect pass/fail)
        optional_viz = ['plotly', 'py3dmol', 'nglview']
        
        for package in optional_viz:
            try:
                __import__(package)
                self.print_status(f"{package}: available", 'pass', 1)
            except ImportError:
                self.print_status(f"{package}: not available (optional)", 'warn', 1)
        
        self.results['visualization'] = all(tests)
        return all(tests)
    
    def test_basic_dft(self):
        """Test basic DFT calculation"""
        self.print_status("Testing Basic DFT Functionality", 'info')
        
        try:
            from ase import Atoms
            from gpaw import GPAW, PW
            import numpy as np
            
            # Create simple H2 molecule
            h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
            h2.center(vacuum=6.0)
            
            # Set up very basic GPAW calculator
            with tempfile.TemporaryDirectory() as tmpdir:
                calc = GPAW(
                    mode=PW(300),  # Higher cutoff for better convergence
                    xc='PBE',
                    kpts=(1, 1, 1),
                    txt=os.path.join(tmpdir, 'h2_test.out'),
                    maxiter=50,  # More iterations for convergence
                    symmetry='off',  # Avoid symmetry issues
                    convergence={'energy': 0.1}  # Looser convergence for test
                )
                
                h2.calc = calc
                
                # Try to calculate energy
                start_time = time.time()
                energy = h2.get_potential_energy()
                calc_time = time.time() - start_time
                
                # Check if energy is reasonable for H2
                if -10.0 < energy < -5.0:
                    self.print_status(f"H2 DFT test: E={energy:.3f} eV ({calc_time:.1f}s)", 'pass', 1)
                    self.results['basic_dft'] = True
                    return True
                else:
                    self.print_status(f"H2 DFT test: unreasonable energy {energy:.3f} eV", 'warn', 1)
                    self.results['basic_dft'] = False
                    return False
                    
        except Exception as e:
            self.print_status(f"Basic DFT test failed: {e}", 'fail', 1)
            self.results['basic_dft'] = False
            return False
    
    def test_structure_files(self):
        """Test MLD structure files"""
        self.print_status("Testing MLD Structure Files", 'info')
        
        required_files = [
            'structures/tma_molecule.xyz',
            'structures/butyne_diol_molecule.xyz',
            'structures/simple_si_surface.xyz'
        ]
        
        tests = []
        
        for file_path in required_files:
            if os.path.exists(file_path):
                try:
                    from ase.io import read
                    atoms = read(file_path)
                    n_atoms = len(atoms)
                    symbols = set(atoms.get_chemical_symbols())
                    
                    self.print_status(f"{file_path}: {n_atoms} atoms ({', '.join(symbols)})", 'pass', 1)
                    tests.append(True)
                    
                except Exception as e:
                    self.print_status(f"{file_path}: read error - {e}", 'fail', 1)
                    tests.append(False)
            else:
                self.print_status(f"{file_path}: missing", 'fail', 1)
                tests.append(False)
        
        self.results['structure_files'] = all(tests)
        return all(tests)
    
    def test_progress_tracking(self):
        """Test progress tracking functionality"""
        self.print_status("Testing Progress Tracking", 'info')
        
        try:
            from tqdm.auto import tqdm
            import time
            
            # Test basic progress bar
            with tqdm(total=10, desc="Test Progress", leave=False) as pbar:
                for i in range(10):
                    time.sleep(0.01)  # Small delay
                    pbar.update(1)
            
            self.print_status("Progress bars: working", 'pass', 1)
            self.results['progress_tracking'] = True
            return True
            
        except Exception as e:
            self.print_status(f"Progress tracking test failed: {e}", 'fail', 1)
            self.results['progress_tracking'] = False
            return False
    
    def generate_report(self):
        """Generate comprehensive validation report"""
        elapsed_time = time.time() - self.start_time
        
        print(f"\n{Colors.BOLD}{'='*80}")
        print("CLEAN MLD INSTALLATION VALIDATION REPORT")
        print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Validation time: {elapsed_time:.1f} seconds")
        print(f"{'='*80}{Colors.END}")
        
        # Environment summary
        print(f"\n{Colors.CYAN}ENVIRONMENT INFORMATION:{Colors.END}")
        print(f"  Platform: {self.env_info['platform']}")
        print(f"  Python: {self.env_info['python_version']}")
        print(f"  Conda Environment: {self.env_info['conda_env']}")
        print(f"  WSL2: {'Yes' if self.env_info['is_wsl'] else 'No'}")
        print(f"  Matplotlib Backend: {self.env_info['mpl_backend']}")
        print(f"  GPAW Setup Path: {self.env_info['gpaw_setup_path']}")
        
        # Test results summary
        print(f"\n{Colors.CYAN}TEST RESULTS:{Colors.END}")
        
        total_tests = len(self.results)
        passed_tests = sum(self.results.values())
        
        for test_name, passed in self.results.items():
            status_icon = "✅" if passed else "❌"
            test_display = test_name.replace('_', ' ').title()
            print(f"  {status_icon} {test_display}")
        
        # Overall status
        success_rate = (passed_tests / total_tests * 100) if total_tests > 0 else 0
        
        print(f"\n{Colors.BOLD}OVERALL STATUS:{Colors.END}")
        print(f"  Tests passed: {passed_tests}/{total_tests}")
        print(f"  Success rate: {success_rate:.1f}%")
        
        if success_rate >= 90:
            print(f"{Colors.GREEN}🎉 EXCELLENT - Installation fully functional{Colors.END}")
            status = "excellent"
        elif success_rate >= 75:
            print(f"{Colors.YELLOW}✅ GOOD - Minor issues, mostly functional{Colors.END}")
            status = "good"
        elif success_rate >= 50:
            print(f"{Colors.YELLOW}⚠️  FAIR - Some issues, basic functionality works{Colors.END}")
            status = "fair"
        else:
            print(f"{Colors.RED}❌ POOR - Major issues, needs attention{Colors.END}")
            status = "poor"
        
        # Recommendations
        print(f"\n{Colors.BOLD}NEXT STEPS:{Colors.END}")
        
        if status == "excellent":
            print(f"  1. Run: python quick_run.py --test")
            print(f"  2. Start thermal MLD simulations")
            print(f"  3. Use: python enhanced_monitor.py for monitoring")
        elif status in ["good", "fair"]:
            print(f"  1. Review failed tests above")
            print(f"  2. Try: python quick_run.py --test")
            print(f"  3. Check environment with: conda list")
        else:
            print(f"  1. Re-run: ./clean_install_mld.sh")
            print(f"  2. Check system dependencies")
            print(f"  3. Verify conda environment activation")
        
        return success_rate >= 75

def main():
    """Main validation function"""
    
    parser = argparse.ArgumentParser(
        description='Validate Clean MLD Installation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python validate_clean_install.py              # Full validation
  python validate_clean_install.py --verbose    # Detailed output
  python validate_clean_install.py --quick      # Skip DFT test
        """
    )
    
    parser.add_argument('--verbose', action='store_true',
                       help='Verbose output with details')
    parser.add_argument('--quick', action='store_true',
                       help='Skip time-consuming DFT test')
    
    args = parser.parse_args()
    
    # Create validator
    validator = CleanInstallValidator(verbose=args.verbose)
    
    print(f"{Colors.BOLD}{Colors.BLUE}")
    print("🔍 MLD CLEAN INSTALLATION VALIDATION")
    print("=" * 50)
    print(f"{Colors.END}")
    
    # Run validation tests
    tests = [
        validator.test_environment_setup,
        validator.test_core_packages,
        validator.test_gpaw_setup,
        validator.test_mpi_functionality,
        validator.test_visualization,
        validator.test_structure_files,
        validator.test_progress_tracking
    ]
    
    # Add DFT test unless quick mode
    if not args.quick:
        tests.append(validator.test_basic_dft)
    
    # Run all tests
    for test_func in tests:
        try:
            test_func()
            print()  # Add spacing
        except Exception as e:
            validator.print_status(f"Test {test_func.__name__} failed: {e}", 'fail')
            validator.results[test_func.__name__] = False
    
    # Generate final report
    success = validator.generate_report()
    
    return 0 if success else 1

if __name__ == "__main__":
    exit(main())