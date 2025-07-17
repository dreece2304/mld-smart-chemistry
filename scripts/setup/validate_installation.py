#!/usr/bin/env python3
"""
Comprehensive MLD Installation Validator
Tests all components systematically with detailed reporting and troubleshooting
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

# Add src to path so we can import our modules
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

try:
    from mld_chemistry.utils import (
        get_system_info, 
        validate_installation as utils_validate,
        print_system_summary
    )
    MLD_MODULES_AVAILABLE = True
except ImportError:
    MLD_MODULES_AVAILABLE = False

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

class MLDValidator:
    """Comprehensive MLD installation validator"""
    
    def __init__(self, verbose=True, quick_mode=False):
        self.verbose = verbose
        self.quick_mode = quick_mode
        self.results = {}
        self.start_time = time.time()
        self.validation_log = []
        
    def log(self, message, level="INFO"):
        """Log validation steps"""
        timestamp = datetime.now().strftime("%H:%M:%S")
        log_entry = f"[{timestamp}] {level}: {message}"
        self.validation_log.append(log_entry)
        
        if self.verbose:
            icons = {"INFO": "ℹ️", "SUCCESS": "✅", "WARNING": "⚠️", "ERROR": "❌", "STEP": "🔧"}
            print(f"{icons.get(level, 'ℹ️')} {message}")
    
    def test_system_environment(self):
        """Test basic system environment"""
        self.log("Testing System Environment", "STEP")
        
        tests_passed = 0
        total_tests = 0
        
        # Python version
        total_tests += 1
        py_version = tuple(map(int, platform.python_version().split('.')))
        if py_version >= (3, 8) and py_version <= (3, 11):
            self.log(f"Python {platform.python_version()}: Compatible", "SUCCESS")
            tests_passed += 1
        else:
            self.log(f"Python {platform.python_version()}: May cause issues", "WARNING")
        
        # Platform detection
        total_tests += 1
        if platform.system() == 'Linux':
            self.log(f"Platform: Linux ({platform.platform()})", "SUCCESS")
            tests_passed += 1
        else:
            self.log(f"Platform: {platform.system()} (Linux required)", "ERROR")
        
        # WSL detection and optimization
        total_tests += 1
        if 'Microsoft' in platform.uname().release:
            self.log("Environment: WSL2 detected", "INFO")
            
            # Check WSL optimizations
            optimizations = [
                ('MPLBACKEND', 'Agg'),
                ('QT_QPA_PLATFORM', 'offscreen'),
                ('OMPI_MCA_btl_vader_single_copy_mechanism', 'none')
            ]
            
            wsl_optimized = True
            for env_var, expected in optimizations:
                if os.environ.get(env_var) != expected:
                    wsl_optimized = False
                    break
            
            if wsl_optimized:
                self.log("WSL2 optimizations: Applied", "SUCCESS")
                tests_passed += 1
            else:
                self.log("WSL2 optimizations: Missing (may affect performance)", "WARNING")
        else:
            self.log("Environment: Native Linux", "SUCCESS")
            tests_passed += 1
        
        # Memory check
        total_tests += 1
        try:
            import psutil
            memory_gb = psutil.virtual_memory().total / (1024**3)
            if memory_gb >= 8:
                self.log(f"Memory: {memory_gb:.1f} GB (sufficient)", "SUCCESS")
                tests_passed += 1
            elif memory_gb >= 4:
                self.log(f"Memory: {memory_gb:.1f} GB (limited, may be slow)", "WARNING")
                tests_passed += 1
            else:
                self.log(f"Memory: {memory_gb:.1f} GB (insufficient)", "ERROR")
        except ImportError:
            self.log("Memory check: psutil not available", "WARNING")
        
        self.results['system_environment'] = tests_passed / total_tests
        return tests_passed == total_tests
    
    def test_core_packages(self):
        """Test core scientific packages"""
        self.log("Testing Core Scientific Packages", "STEP")
        
        required_packages = {
            'numpy': '1.21.0',
            'scipy': '1.7.0', 
            'matplotlib': '3.5.0',
            'ase': '3.22.0',
            'gpaw': '22.8.0',
            'mpi4py': '3.1.0',
            'tqdm': '4.0.0'
        }
        
        optional_packages = {
            'plotly': '5.0.0',
            'py3dmol': '2.0.0',
            'nglview': '3.0.0',
            'psutil': '5.8.0'
        }
        
        passed_required = 0
        passed_optional = 0
        
        # Test required packages
        for package, min_version in required_packages.items():
            try:
                module = __import__(package)
                version = getattr(module, '__version__', 'unknown')
                self.log(f"{package}: {version} ✓", "SUCCESS")
                passed_required += 1
            except ImportError:
                self.log(f"{package}: MISSING (required)", "ERROR")
        
        # Test optional packages
        for package, min_version in optional_packages.items():
            try:
                module = __import__(package)
                version = getattr(module, '__version__', 'unknown')
                self.log(f"{package}: {version} (optional)", "SUCCESS")
                passed_optional += 1
            except ImportError:
                self.log(f"{package}: not installed (optional)", "WARNING")
        
        required_success = passed_required == len(required_packages)
        self.results['core_packages'] = passed_required / len(required_packages)
        self.results['optional_packages'] = passed_optional / len(optional_packages)
        
        return required_success
    
    def test_gpaw_installation(self):
        """Test GPAW installation comprehensively"""
        self.log("Testing GPAW Installation", "STEP")
        
        tests_passed = 0
        total_tests = 0
        
        # Basic GPAW import
        total_tests += 1
        try:
            import gpaw
            self.log(f"GPAW import: {gpaw.__version__}", "SUCCESS")
            tests_passed += 1
        except ImportError as e:
            self.log(f"GPAW import failed: {e}", "ERROR")
            self.results['gpaw'] = 0
            return False
        
        # GPAW setup path
        total_tests += 1
        setup_path = os.environ.get('GPAW_SETUP_PATH', '')
        if setup_path and os.path.exists(setup_path):
            self.log(f"GPAW_SETUP_PATH: {setup_path}", "SUCCESS")
            tests_passed += 1
            
            # Count datasets
            dataset_files = list(Path(setup_path).glob('*.gz'))
            self.log(f"PAW datasets: {len(dataset_files)} files", "INFO")
            
            # Check for key elements
            key_elements = ['H.PBE.gz', 'C.PBE.gz', 'O.PBE.gz', 'Al.PBE.gz', 'Si.PBE.gz']
            missing_elements = []
            
            for element in key_elements:
                if not (Path(setup_path) / element).exists():
                    missing_elements.append(element)
            
            if missing_elements:
                self.log(f"Missing key datasets: {missing_elements}", "WARNING")
            else:
                self.log("All key element datasets found", "SUCCESS")
                
        else:
            self.log(f"GPAW_SETUP_PATH not set or invalid: {setup_path}", "ERROR")
        
        # Test GPAW calculator creation
        if not self.quick_mode:
            total_tests += 1
            try:
                from gpaw import GPAW, PW
                # Try to create calculator (this tests compilation)
                calc = GPAW(mode=PW(200), xc='PBE', txt=None)
                self.log("GPAW calculator creation: Working", "SUCCESS")
                tests_passed += 1
            except Exception as e:
                self.log(f"GPAW calculator creation failed: {e}", "ERROR")
        
        self.results['gpaw'] = tests_passed / total_tests
        return tests_passed == total_tests
    
    def test_mpi_functionality(self):
        """Test MPI parallel functionality"""
        self.log("Testing MPI Functionality", "STEP")
        
        tests_passed = 0
        total_tests = 0
        
        # Basic MPI import
        total_tests += 1
        try:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            size = comm.Get_size()
            rank = comm.Get_rank()
            
            self.log(f"MPI: rank {rank}/{size}", "SUCCESS")
            tests_passed += 1
        except ImportError as e:
            self.log(f"MPI import failed: {e}", "ERROR")
            self.results['mpi'] = 0
            return False
        
        # Test MPI execution with subprocess
        if not self.quick_mode:
            total_tests += 1
            try:
                result = subprocess.run([
                    'mpirun', '-np', '2', 'python', '-c',
                    'from mpi4py import MPI; print(f"Process {MPI.COMM_WORLD.Get_rank()}")'
                ], capture_output=True, text=True, timeout=10)
                
                if result.returncode == 0 and len(result.stdout.strip().split('\n')) == 2:
                    self.log("MPI multiprocess execution: Working", "SUCCESS")
                    tests_passed += 1
                else:
                    self.log("MPI multiprocess execution: Issues detected", "WARNING")
                    
            except (subprocess.TimeoutExpired, FileNotFoundError) as e:
                self.log(f"MPI multiprocess test failed: {e}", "WARNING")
        
        self.results['mpi'] = tests_passed / total_tests
        return tests_passed > 0  # At least basic MPI should work
    
    def test_visualization(self):
        """Test visualization capabilities"""
        self.log("Testing Visualization Capabilities", "STEP")
        
        tests_passed = 0
        total_tests = 0
        
        # Test matplotlib headless
        total_tests += 1
        try:
            import matplotlib
            matplotlib.use('Agg')  # Force headless
            import matplotlib.pyplot as plt
            import numpy as np
            
            # Create test plot
            fig, ax = plt.subplots()
            x = np.linspace(0, 10, 100)
            ax.plot(x, np.sin(x))
            
            # Save to temporary file
            with tempfile.NamedTemporaryFile(suffix='.png', delete=True) as tmp:
                fig.savefig(tmp.name)
                plt.close(fig)
                
                if os.path.getsize(tmp.name) > 1000:  # File has content
                    self.log("Matplotlib (headless): Working", "SUCCESS")
                    tests_passed += 1
                else:
                    self.log("Matplotlib (headless): Empty output", "ERROR")
                    
        except Exception as e:
            self.log(f"Matplotlib test failed: {e}", "ERROR")
        
        # Test ASE visualization
        total_tests += 1
        try:
            from ase import Atoms
            from ase.visualize.plot import plot_atoms
            
            # Create simple molecule
            h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
            
            fig, ax = plt.subplots()
            plot_atoms(h2, ax)
            
            with tempfile.NamedTemporaryFile(suffix='.png', delete=True) as tmp:
                fig.savefig(tmp.name)
                plt.close(fig)
                
                if os.path.getsize(tmp.name) > 1000:
                    self.log("ASE visualization: Working", "SUCCESS")
                    tests_passed += 1
                else:
                    self.log("ASE visualization: Empty output", "ERROR")
                    
        except Exception as e:
            self.log(f"ASE visualization test failed: {e}", "ERROR")
        
        self.results['visualization'] = tests_passed / total_tests
        return tests_passed > 0
    
    def test_basic_dft_calculation(self):
        """Test basic DFT calculation"""
        if self.quick_mode:
            self.log("Skipping DFT test (quick mode)", "INFO")
            return True
            
        self.log("Testing Basic DFT Calculation", "STEP")
        
        try:
            from ase import Atoms
            from gpaw import GPAW, PW
            
            # Create simple H2 molecule
            h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
            h2.center(vacuum=6.0)
            
            # Set up minimal GPAW calculator
            with tempfile.TemporaryDirectory() as tmpdir:
                calc = GPAW(
                    mode=PW(200),  # Low cutoff for speed
                    xc='PBE',
                    kpts=(1, 1, 1),
                    txt=os.path.join(tmpdir, 'h2_test.out'),
                    maxiter=20,  # Limit iterations
                    symmetry='off'
                )
                
                h2.calc = calc
                
                # Calculate energy
                start_time = time.time()
                energy = h2.get_potential_energy()
                calc_time = time.time() - start_time
                
                # Sanity check energy
                if -10.0 < energy < -5.0:  # Reasonable H2 energy range
                    self.log(f"H2 DFT calculation: E={energy:.3f} eV ({calc_time:.1f}s)", "SUCCESS")
                    self.results['basic_dft'] = True
                    return True
                else:
                    self.log(f"H2 DFT calculation: Unreasonable energy {energy:.3f} eV", "WARNING")
                    self.results['basic_dft'] = False
                    return False
                    
        except Exception as e:
            self.log(f"Basic DFT calculation failed: {e}", "ERROR")
            self.results['basic_dft'] = False
            return False
    
    def test_mld_modules(self):
        """Test MLD-specific modules"""
        self.log("Testing MLD Framework Modules", "STEP")
        
        if not MLD_MODULES_AVAILABLE:
            self.log("MLD modules not in path - testing basic structure", "WARNING")
            
            # Check if files exist
            src_dir = Path(__file__).parent.parent / 'src' / 'mld_chemistry'
            if src_dir.exists():
                required_files = ['__init__.py', 'smart_optimizer.py', 'visualization.py', 'utils.py']
                missing_files = [f for f in required_files if not (src_dir / f).exists()]
                
                if missing_files:
                    self.log(f"Missing MLD module files: {missing_files}", "ERROR")
                    self.results['mld_modules'] = False
                    return False
                else:
                    self.log("MLD module files: All present", "SUCCESS")
                    self.results['mld_modules'] = True
                    return True
            else:
                self.log("MLD source directory not found", "ERROR")
                self.results['mld_modules'] = False
                return False
        
        try:
            # Test imports
            from mld_chemistry import SmartOptimizer, optimize_molecule
            from mld_chemistry.visualization import visualize_molecule
            from mld_chemistry.utils import get_system_info
            
            self.log("MLD module imports: Working", "SUCCESS")
            
            # Test basic functionality
            system_info = get_system_info()
            if isinstance(system_info, dict) and 'platform' in system_info:
                self.log("MLD utilities: Working", "SUCCESS")
                self.results['mld_modules'] = True
                return True
            else:
                self.log("MLD utilities: Issues detected", "WARNING")
                self.results['mld_modules'] = False
                return False
                
        except Exception as e:
            self.log(f"MLD module test failed: {e}", "ERROR")
            self.results['mld_modules'] = False
            return False
    
    def generate_report(self):
        """Generate comprehensive validation report"""
        elapsed_time = time.time() - self.start_time
        
        print(f"\n{Colors.BOLD}{'='*80}")
        print("MLD INSTALLATION VALIDATION REPORT")
        print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Validation time: {elapsed_time:.1f} seconds")
        print(f"{'='*80}{Colors.END}")
        
        # System information
        if MLD_MODULES_AVAILABLE:
            print(f"\n{Colors.CYAN}SYSTEM INFORMATION:{Colors.END}")
            try:
                print_system_summary()
            except:
                print("System information not available")
        
        # Test results
        print(f"\n{Colors.CYAN}VALIDATION RESULTS:{Colors.END}")
        
        test_results = [
            ('System Environment', self.results.get('system_environment', 0)),
            ('Core Packages', self.results.get('core_packages', 0)),
            ('GPAW Installation', self.results.get('gpaw', 0)),
            ('MPI Functionality', self.results.get('mpi', 0)),
            ('Visualization', self.results.get('visualization', 0)),
            ('MLD Modules', self.results.get('mld_modules', 0))
        ]
        
        if not self.quick_mode:
            test_results.append(('Basic DFT', self.results.get('basic_dft', 0)))
        
        total_score = 0
        max_score = 0
        
        for test_name, score in test_results:
            if isinstance(score, bool):
                score = 1.0 if score else 0.0
            
            total_score += score
            max_score += 1.0
            
            if score >= 0.9:
                status = f"{Colors.GREEN}✅ PASS{Colors.END}"
            elif score >= 0.5:
                status = f"{Colors.YELLOW}⚠️  PARTIAL{Colors.END}"
            else:
                status = f"{Colors.RED}❌ FAIL{Colors.END}"
            
            print(f"  {test_name:<20} {status} ({score*100:.0f}%)")
        
        # Overall assessment
        overall_score = (total_score / max_score) * 100 if max_score > 0 else 0
        
        print(f"\n{Colors.BOLD}OVERALL ASSESSMENT:{Colors.END}")
        print(f"  Score: {overall_score:.1f}%")
        
        if overall_score >= 90:
            print(f"{Colors.GREEN}🎉 EXCELLENT - Ready for production use{Colors.END}")
            status = "excellent"
        elif overall_score >= 75:
            print(f"{Colors.YELLOW}✅ GOOD - Ready for most tasks{Colors.END}")
            status = "good"
        elif overall_score >= 50:
            print(f"{Colors.YELLOW}⚠️  FAIR - Basic functionality works{Colors.END}")
            status = "fair"
        else:
            print(f"{Colors.RED}❌ POOR - Major issues need attention{Colors.END}")
            status = "poor"
        
        # Recommendations
        print(f"\n{Colors.BOLD}RECOMMENDATIONS:{Colors.END}")
        
        if status == "excellent":
            print("  ✨ Installation is complete and optimized!")
            print("  🚀 Run: python scripts/step1_create_h2o.py")
            print("  📊 Try: python scripts/step2_smart_tma.py")
        elif status in ["good", "fair"]:
            print("  🔧 Review failed tests above")
            print("  📖 Check documentation for troubleshooting")
            print("  🧪 Test basic functionality: python scripts/step1_create_h2o.py")
        else:
            print("  🔄 Re-run installer: ./install/setup_complete_environment.sh")
            print("  📋 Check system requirements")
            print("  💬 Seek support if issues persist")
        
        print(f"\n{Colors.BOLD}SUPPORT:{Colors.END}")
        print("  📖 Documentation: docs/troubleshooting.md")
        print("  🐛 Issues: https://github.com/yourusername/mld-smart-chemistry/issues")
        print("  💬 Discussions: https://github.com/yourusername/mld-smart-chemistry/discussions")
        
        return overall_score >= 75

def main():
    """Main validation function"""
    
    parser = argparse.ArgumentParser(
        description='Validate MLD Installation',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--quick', action='store_true',
                       help='Skip time-consuming tests (DFT calculation)')
    parser.add_argument('--quiet', action='store_true',
                       help='Minimal output')
    parser.add_argument('--log-file', type=str,
                       help='Save validation log to file')
    
    args = parser.parse_args()
    
    # Create validator
    validator = MLDValidator(verbose=not args.quiet, quick_mode=args.quick)
    
    if not args.quiet:
        print(f"{Colors.BOLD}{Colors.BLUE}")
        print("🔍 MLD INSTALLATION VALIDATION")
        print("=" * 50)
        print(f"{Colors.END}")
    
    # Run validation tests
    tests = [
        validator.test_system_environment,
        validator.test_core_packages,
        validator.test_gpaw_installation,
        validator.test_mpi_functionality,
        validator.test_visualization,
        validator.test_mld_modules
    ]
    
    if not args.quick:
        tests.append(validator.test_basic_dft_calculation)
    
    # Execute tests
    for test_func in tests:
        try:
            test_func()
            if not args.quiet:
                print()  # Add spacing
        except Exception as e:
            validator.log(f"Test {test_func.__name__} crashed: {e}", "ERROR")
    
    # Generate report
    success = validator.generate_report()
    
    # Save log if requested
    if args.log_file:
        with open(args.log_file, 'w') as f:
            f.write('\n'.join(validator.validation_log))
        print(f"\n📝 Validation log saved to: {args.log_file}")
    
    return 0 if success else 1

if __name__ == "__main__":
    exit(main())