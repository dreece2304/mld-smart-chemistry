#!/usr/bin/env python3
"""
Comprehensive MLD Setup Validation Script
Validates all dependencies, configurations, and scripts before running thermal MLD simulations
"""

import os
import sys
import subprocess
import importlib
import platform
import psutil
from pathlib import Path
from datetime import datetime
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

class MLDValidator:
    """Comprehensive MLD setup validator"""
    
    def __init__(self, fix_issues=False, verbose=False):
        self.fix_issues = fix_issues
        self.verbose = verbose
        self.issues = []
        self.warnings = []
        self.fixes_applied = []
        
        # System info
        self.system_info = {
            'platform': platform.system(),
            'python_version': platform.python_version(),
            'cpu_count': psutil.cpu_count(),
            'memory_gb': psutil.virtual_memory().total / (1024**3),
            'conda_env': os.environ.get('CONDA_DEFAULT_ENV', 'none'),
            'conda_prefix': os.environ.get('CONDA_PREFIX', 'none')
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
        elif status == 'fix':
            print(f"{indent_str}{Colors.MAGENTA}🔧 {message}{Colors.END}")
        else:
            print(f"{indent_str}{message}")
    
    def run_command(self, cmd, capture_output=True):
        """Run shell command safely"""
        try:
            result = subprocess.run(cmd, shell=True, capture_output=capture_output, 
                                  text=True, timeout=30)
            return result.returncode == 0, result.stdout, result.stderr
        except subprocess.TimeoutExpired:
            return False, "", "Command timed out"
        except Exception as e:
            return False, "", str(e)
    
    def check_system_requirements(self):
        """Check basic system requirements"""
        self.print_status("Checking System Requirements", 'info')
        
        # Python version
        py_version = tuple(map(int, platform.python_version().split('.')))
        if py_version >= (3, 8) and py_version <= (3, 11):
            self.print_status(f"Python {platform.python_version()}", 'pass', 1)
        else:
            self.issues.append(f"Python version {platform.python_version()} not recommended (use 3.8-3.11)")
            self.print_status(f"Python {platform.python_version()} (not recommended)", 'warn', 1)
        
        # Memory check
        memory_gb = self.system_info['memory_gb']
        if memory_gb >= 8:
            self.print_status(f"Memory: {memory_gb:.1f} GB", 'pass', 1)
        elif memory_gb >= 4:
            self.warnings.append(f"Low memory: {memory_gb:.1f} GB (recommended: 8+ GB)")
            self.print_status(f"Memory: {memory_gb:.1f} GB (low)", 'warn', 1)
        else:
            self.issues.append(f"Insufficient memory: {memory_gb:.1f} GB (minimum: 4 GB)")
            self.print_status(f"Memory: {memory_gb:.1f} GB (insufficient)", 'fail', 1)
        
        # CPU check
        cpu_count = self.system_info['cpu_count']
        if cpu_count >= 4:
            self.print_status(f"CPUs: {cpu_count}", 'pass', 1)
        else:
            self.warnings.append(f"Low CPU count: {cpu_count} (recommended: 4+)")
            self.print_status(f"CPUs: {cpu_count} (low)", 'warn', 1)
        
        return len(self.issues) == 0
    
    def check_conda_environment(self):
        """Check conda environment setup"""
        self.print_status("Checking Conda Environment", 'info')
        
        # Check if conda is available
        success, stdout, stderr = self.run_command("conda --version")
        if not success:
            self.issues.append("Conda not found or not working")
            self.print_status("Conda not available", 'fail', 1)
            return False
        
        conda_version = stdout.strip()
        self.print_status(f"Conda version: {conda_version}", 'pass', 1)
        
        # Check current environment
        current_env = self.system_info['conda_env']
        if current_env == 'mld_modeling':
            self.print_status("Currently in mld_modeling environment", 'pass', 1)
        elif current_env == 'base':
            self.warnings.append("Currently in base environment (should use mld_modeling)")
            self.print_status("In base environment (recommend mld_modeling)", 'warn', 1)
            
            if self.fix_issues:
                self.print_status("Attempting to activate mld_modeling environment", 'fix', 1)
                # Note: We can't actually change environment from Python, but we can check if it exists
                success, stdout, stderr = self.run_command("conda env list")
                if 'mld_modeling' in stdout:
                    self.print_status("mld_modeling environment exists - please activate manually", 'info', 2)
                    self.fixes_applied.append("Found mld_modeling environment")
                else:
                    self.issues.append("mld_modeling environment not found")
                    self.print_status("mld_modeling environment not found", 'fail', 2)
        else:
            self.warnings.append(f"In unknown environment: {current_env}")
            self.print_status(f"Unknown environment: {current_env}", 'warn', 1)
        
        return True
    
    def check_python_packages(self):
        """Check required Python packages"""
        self.print_status("Checking Python Packages", 'info')
        
        # Core packages with version requirements
        core_packages = {
            'ase': '3.22.0',
            'gpaw': '22.8.0', 
            'numpy': '1.21.0',
            'scipy': '1.7.0',
            'matplotlib': '3.5.0',
            'tqdm': '4.0.0'
        }
        
        # Optional packages
        optional_packages = {
            'pymatgen': '2023.1.0',
            'spglib': '2.0.0',
            'ovito': '3.8.0',
            'mpi4py': '3.1.0'
        }
        
        missing_core = []
        missing_optional = []
        
        # Check core packages
        for package, min_version in core_packages.items():
            try:
                module = importlib.import_module(package)
                version = getattr(module, '__version__', 'unknown')
                self.print_status(f"{package}: {version}", 'pass', 1)
            except ImportError:
                missing_core.append(package)
                self.print_status(f"{package}: not installed", 'fail', 1)
        
        # Check optional packages
        for package, min_version in optional_packages.items():
            try:
                module = importlib.import_module(package)
                version = getattr(module, '__version__', 'unknown')
                self.print_status(f"{package}: {version} (optional)", 'pass', 1)
            except ImportError:
                missing_optional.append(package)
                self.print_status(f"{package}: not installed (optional)", 'warn', 1)
        
        # Handle missing packages
        if missing_core:
            self.issues.append(f"Missing core packages: {', '.join(missing_core)}")
            
            if self.fix_issues:
                self.print_status("Installing missing core packages", 'fix', 1)
                for package in missing_core:
                    self.print_status(f"Installing {package}...", 'fix', 2)
                    success, stdout, stderr = self.run_command(f"pip install {package}")
                    if success:
                        self.fixes_applied.append(f"Installed {package}")
                        self.print_status(f"Successfully installed {package}", 'pass', 2)
                    else:
                        self.print_status(f"Failed to install {package}: {stderr}", 'fail', 2)
        
        if missing_optional:
            self.warnings.append(f"Missing optional packages: {', '.join(missing_optional)}")
        
        return len(missing_core) == 0
    
    def check_gpaw_setup(self):
        """Check GPAW installation and datasets"""
        self.print_status("Checking GPAW Setup", 'info')
        
        # Check GPAW import
        try:
            import gpaw
            self.print_status(f"GPAW version: {gpaw.__version__}", 'pass', 1)
        except ImportError:
            self.issues.append("GPAW not installed")
            self.print_status("GPAW not installed", 'fail', 1)
            return False
        
        # Check GPAW_SETUP_PATH
        gpaw_path = os.environ.get('GPAW_SETUP_PATH')
        expected_path = os.path.expanduser('~/.local/share/gpaw/gpaw-setups-24.11.0')
        
        if gpaw_path:
            self.print_status(f"GPAW_SETUP_PATH: {gpaw_path}", 'pass', 1)
            if gpaw_path != expected_path:
                self.warnings.append(f"GPAW_SETUP_PATH differs from expected: {expected_path}")
        else:
            self.warnings.append("GPAW_SETUP_PATH not set")
            self.print_status("GPAW_SETUP_PATH not set", 'warn', 1)
            
            if self.fix_issues:
                self.print_status("Setting GPAW_SETUP_PATH", 'fix', 1)
                os.environ['GPAW_SETUP_PATH'] = expected_path
                self.fixes_applied.append("Set GPAW_SETUP_PATH environment variable")
        
        # Check datasets directory
        datasets_path = gpaw_path or expected_path
        if os.path.exists(datasets_path):
            # Count dataset files
            dataset_files = list(Path(datasets_path).glob('*.py'))
            self.print_status(f"GPAW datasets: {len(dataset_files)} files found", 'pass', 1)
        else:
            self.issues.append(f"GPAW datasets not found at {datasets_path}")
            self.print_status(f"GPAW datasets not found", 'fail', 1)
            
            if self.fix_issues:
                self.print_status("Attempting to install GPAW datasets", 'fix', 1)
                success, stdout, stderr = self.run_command("gpaw install-data ~/.local/share/gpaw")
                if success:
                    self.fixes_applied.append("Installed GPAW datasets")
                    self.print_status("Successfully installed GPAW datasets", 'pass', 2)
                else:
                    self.print_status(f"Failed to install datasets: {stderr}", 'fail', 2)
        
        return True
    
    def check_structure_files(self):
        """Check structure files integrity"""
        self.print_status("Checking Structure Files", 'info')
        
        required_files = [
            'structures/tma_molecule.xyz',
            'structures/butyne_diol_molecule.xyz', 
            'structures/simple_si_surface.xyz'
        ]
        
        optional_files = [
            'structures/diol_improved.xyz',
            'structures/initial_mld_system.xyz'
        ]
        
        all_good = True
        
        # Check required files
        for file_path in required_files:
            if os.path.exists(file_path):
                # Check file size and content
                size = os.path.getsize(file_path)
                if size > 0:
                    self.print_status(f"{file_path}: {size} bytes", 'pass', 1)
                    
                    # Basic validation: check if it's a valid XYZ file
                    try:
                        with open(file_path, 'r') as f:
                            first_line = f.readline().strip()
                            if first_line.isdigit():
                                self.print_status(f"Valid XYZ format", 'pass', 2)
                            else:
                                self.warnings.append(f"{file_path} may not be valid XYZ format")
                                self.print_status(f"Questionable XYZ format", 'warn', 2)
                    except Exception as e:
                        self.warnings.append(f"Could not validate {file_path}: {e}")
                        self.print_status(f"Could not validate format", 'warn', 2)
                else:
                    self.issues.append(f"{file_path} is empty")
                    self.print_status(f"{file_path}: empty file", 'fail', 1)
                    all_good = False
            else:
                self.issues.append(f"Missing required file: {file_path}")
                self.print_status(f"{file_path}: missing", 'fail', 1)
                all_good = False
        
        # Check optional files
        for file_path in optional_files:
            if os.path.exists(file_path):
                size = os.path.getsize(file_path)
                self.print_status(f"{file_path}: {size} bytes (optional)", 'pass', 1)
            else:
                self.print_status(f"{file_path}: not found (optional)", 'warn', 1)
        
        return all_good
    
    def check_script_integrity(self):
        """Check Python scripts for syntax errors"""
        self.print_status("Checking Script Integrity", 'info')
        
        scripts = [
            'quick_run.py',
            'thermal_mld_simulation.py',
            'progress_optimization.py',
            'run_mld_with_progress.py',
            'fix_diol_structure.py'
        ]
        
        all_good = True
        
        for script in scripts:
            if os.path.exists(script):
                try:
                    # Try to compile the script
                    with open(script, 'r') as f:
                        source = f.read()
                    compile(source, script, 'exec')
                    self.print_status(f"{script}: syntax OK", 'pass', 1)
                except SyntaxError as e:
                    self.issues.append(f"Syntax error in {script}: {e}")
                    self.print_status(f"{script}: syntax error", 'fail', 1)
                    all_good = False
                except Exception as e:
                    self.warnings.append(f"Could not check {script}: {e}")
                    self.print_status(f"{script}: could not check", 'warn', 1)
            else:
                self.issues.append(f"Missing script: {script}")
                self.print_status(f"{script}: missing", 'fail', 1)
                all_good = False
        
        return all_good
    
    def run_basic_functionality_test(self):
        """Run basic functionality tests"""
        self.print_status("Running Basic Functionality Tests", 'info')
        
        tests_passed = 0
        total_tests = 3
        
        # Test 1: Import test
        try:
            import ase
            import gpaw
            import numpy as np
            from ase import Atoms
            from gpaw import GPAW, PW
            self.print_status("Package imports: OK", 'pass', 1)
            tests_passed += 1
        except Exception as e:
            self.issues.append(f"Import test failed: {e}")
            self.print_status(f"Package imports: FAILED", 'fail', 1)
        
        # Test 2: ASE atoms creation
        try:
            from ase import Atoms
            h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
            h2.center(vacuum=5.0)
            self.print_status("ASE Atoms creation: OK", 'pass', 1)
            tests_passed += 1
        except Exception as e:
            self.issues.append(f"ASE test failed: {e}")
            self.print_status(f"ASE Atoms creation: FAILED", 'fail', 1)
        
        # Test 3: GPAW calculator creation (without running)
        try:
            from gpaw import GPAW, PW
            calc = GPAW(mode=PW(200), xc='PBE', kpts=(1,1,1), txt=None)
            self.print_status("GPAW calculator creation: OK", 'pass', 1)
            tests_passed += 1
        except Exception as e:
            self.issues.append(f"GPAW test failed: {e}")
            self.print_status(f"GPAW calculator creation: FAILED", 'fail', 1)
        
        self.print_status(f"Basic tests: {tests_passed}/{total_tests} passed", 
                         'pass' if tests_passed == total_tests else 'warn', 1)
        
        return tests_passed == total_tests
    
    def generate_report(self):
        """Generate comprehensive validation report"""
        print(f"\n{Colors.BOLD}{'='*80}")
        print(f"MLD SETUP VALIDATION REPORT")
        print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"{'='*80}{Colors.END}")
        
        # System summary
        print(f"\n{Colors.CYAN}SYSTEM INFORMATION:{Colors.END}")
        print(f"  Platform: {self.system_info['platform']}")
        print(f"  Python: {self.system_info['python_version']}")
        print(f"  CPUs: {self.system_info['cpu_count']}")
        print(f"  Memory: {self.system_info['memory_gb']:.1f} GB")
        print(f"  Conda Environment: {self.system_info['conda_env']}")
        
        # Issues summary
        if self.issues:
            print(f"\n{Colors.RED}CRITICAL ISSUES ({len(self.issues)}):{Colors.END}")
            for i, issue in enumerate(self.issues, 1):
                print(f"  {i}. {issue}")
        
        if self.warnings:
            print(f"\n{Colors.YELLOW}WARNINGS ({len(self.warnings)}):{Colors.END}")
            for i, warning in enumerate(self.warnings, 1):
                print(f"  {i}. {warning}")
        
        if self.fixes_applied:
            print(f"\n{Colors.MAGENTA}FIXES APPLIED ({len(self.fixes_applied)}):{Colors.END}")
            for i, fix in enumerate(self.fixes_applied, 1):
                print(f"  {i}. {fix}")
        
        # Overall status
        print(f"\n{Colors.BOLD}OVERALL STATUS:{Colors.END}")
        if not self.issues:
            print(f"{Colors.GREEN}✅ SETUP VALIDATION PASSED{Colors.END}")
            print(f"   Ready to run MLD simulations!")
        elif len(self.issues) <= 2:
            print(f"{Colors.YELLOW}⚠️  SETUP HAS MINOR ISSUES{Colors.END}")
            print(f"   May still work but recommend fixing issues")
        else:
            print(f"{Colors.RED}❌ SETUP VALIDATION FAILED{Colors.END}")
            print(f"   Must fix issues before running simulations")
        
        if self.warnings:
            print(f"{Colors.YELLOW}⚠️  {len(self.warnings)} warnings detected{Colors.END}")
        
        print(f"\n{Colors.BOLD}NEXT STEPS:{Colors.END}")
        if not self.issues:
            print(f"  1. Run: python quick_run.py --test")
            print(f"  2. Run: python test_suite.py --level 1-3")
            print(f"  3. Start thermal MLD simulations")
        else:
            print(f"  1. Fix critical issues listed above")
            print(f"  2. Re-run: python validate_setup.py --fix-issues")
            print(f"  3. Run test suite when validation passes")
        
        return len(self.issues) == 0

def main():
    """Main validation function"""
    
    parser = argparse.ArgumentParser(
        description='Comprehensive MLD Setup Validation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python validate_setup.py                    # Basic validation
  python validate_setup.py --fix-issues       # Validate and fix issues
  python validate_setup.py --verbose          # Detailed output
  python validate_setup.py --full --fix-issues # Complete validation with fixes
        """
    )
    
    parser.add_argument('--fix-issues', action='store_true',
                       help='Attempt to fix issues automatically')
    parser.add_argument('--verbose', action='store_true',
                       help='Verbose output')
    parser.add_argument('--full', action='store_true',
                       help='Run all validation checks')
    
    args = parser.parse_args()
    
    # Create validator
    validator = MLDValidator(fix_issues=args.fix_issues, verbose=args.verbose)
    
    print(f"{Colors.BOLD}{Colors.BLUE}")
    print("🔍 MLD SETUP VALIDATION")
    print("=" * 50)
    print(f"{Colors.END}")
    
    # Run validation checks
    checks = [
        validator.check_system_requirements,
        validator.check_conda_environment, 
        validator.check_python_packages,
        validator.check_gpaw_setup,
        validator.check_structure_files,
        validator.check_script_integrity
    ]
    
    if args.full:
        checks.append(validator.run_basic_functionality_test)
    
    passed_checks = 0
    for check in checks:
        try:
            if check():
                passed_checks += 1
        except Exception as e:
            validator.issues.append(f"Validation check failed: {e}")
            validator.print_status(f"Check failed: {e}", 'fail')
    
    # Generate report
    success = validator.generate_report()
    
    print(f"\n{Colors.BOLD}Validation completed: {passed_checks}/{len(checks)} checks passed{Colors.END}")
    
    return 0 if success else 1

if __name__ == "__main__":
    exit(main())