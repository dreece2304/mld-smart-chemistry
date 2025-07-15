#!/usr/bin/env python3
"""
DFT Troubleshooting Tool for Surface SCF Convergence Issues
Addresses common GPAW convergence problems with surfaces

Common issues and solutions:
1. SCF convergence failure → Better mixing, level shifting
2. Surface systems harder than molecules → Smearing, different functionals
3. Conservative parameters too restrictive → Adaptive parameter adjustment
"""

import numpy as np
from ase import Atoms
from ase.io import read, write
import os
import time

try:
    from gpaw import GPAW, PW
    from gpaw.mixer import Mixer, MixerSum
    GPAW_AVAILABLE = True
except ImportError:
    GPAW_AVAILABLE = False

class DFTTroubleshooter:
    """Diagnose and fix DFT convergence issues"""
    
    def __init__(self):
        """Initialize with progressive parameter sets"""
        
        # Parameter sets from most conservative to most aggressive
        self.parameter_sets = {
            'ultra_conservative': {
                'mode': PW(200),
                'xc': 'PBE',
                'kpts': (1, 1, 1),
                'convergence': {'energy': 1e-3, 'density': 1e-2, 'eigenstates': 1e-2},
                'mixer': {'beta': 0.01, 'nmaxold': 3},
                'maxiter': 500,
                'occupations': {'name': 'fermi-dirac', 'width': 0.1},
                'txt': 'dft_ultra_conservative.txt'
            },
            
            'conservative': {
                'mode': PW(250),
                'xc': 'PBE',
                'kpts': (1, 1, 1),
                'convergence': {'energy': 1e-4, 'density': 1e-3, 'eigenstates': 1e-3},
                'mixer': {'beta': 0.05, 'nmaxold': 5},
                'maxiter': 400,
                'occupations': {'name': 'fermi-dirac', 'width': 0.05},
                'txt': 'dft_conservative.txt'
            },
            
            'standard': {
                'mode': PW(300),
                'xc': 'PBE', 
                'kpts': (1, 1, 1),
                'convergence': {'energy': 1e-4, 'density': 1e-3, 'eigenstates': 1e-3},
                'mixer': {'beta': 0.1, 'nmaxold': 8},
                'maxiter': 300,
                'occupations': {'name': 'fermi-dirac', 'width': 0.02},
                'txt': 'dft_standard.txt'
            },
            
            'aggressive': {
                'mode': PW(350),
                'xc': 'PBE',
                'kpts': (1, 1, 1),
                'convergence': {'energy': 1e-5, 'density': 1e-4, 'eigenstates': 1e-4},
                'mixer': {'beta': 0.2, 'nmaxold': 10},
                'maxiter': 200,
                'occupations': {'name': 'fermi-dirac', 'width': 0.01},
                'txt': 'dft_aggressive.txt'
            }
        }
        
        print("🔧 DFT Troubleshooter")
        print("   Progressive parameter adjustment for convergence")
        print("   Specialized for surface systems")
        
    def analyze_convergence_log(self, log_file: str) -> dict:
        """Analyze GPAW log file to diagnose convergence issues"""
        
        if not os.path.exists(log_file):
            return {'error': 'Log file not found'}
        
        analysis = {
            'total_iterations': 0,
            'final_energy_diff': None,
            'final_density_diff': None,
            'converged': False,
            'likely_issues': []
        }
        
        try:
            with open(log_file, 'r') as f:
                lines = f.readlines()
            
            # Look for convergence information
            for i, line in enumerate(lines):
                if 'iter:' in line and 'energy' in line:
                    analysis['total_iterations'] += 1
                    
                    # Extract energy and density differences
                    parts = line.split()
                    try:
                        energy_idx = parts.index('energy') + 1
                        if energy_idx < len(parts):
                            analysis['final_energy_diff'] = float(parts[energy_idx])
                    except (ValueError, IndexError):
                        pass
                
                if 'Converged' in line:
                    analysis['converged'] = True
                elif 'Did not converge' in line:
                    analysis['converged'] = False
        
        # Diagnose likely issues
        if analysis['total_iterations'] >= 300:
            analysis['likely_issues'].append("Too many SCF iterations - try better mixing")
        
        if analysis['final_energy_diff'] and analysis['final_energy_diff'] > 1e-3:
            analysis['likely_issues'].append("Energy not converging - try smaller mixing")
            
        if not analysis['converged']:
            analysis['likely_issues'].append("SCF failure - surface systems need special treatment")
        
        return analysis
    
    def test_parameter_set(self, atoms: Atoms, param_name: str) -> tuple:
        """Test a specific parameter set"""
        
        print(f"\n🧪 Testing {param_name} parameters...")
        
        if not GPAW_AVAILABLE:
            return False, 0.0, "GPAW not available"
        
        # Copy atoms and center
        test_atoms = atoms.copy()
        test_atoms.center(vacuum=10.0)
        
        # Create calculator
        params = self.parameter_sets[param_name].copy()
        
        try:
            calc = GPAW(**params)
            test_atoms.calc = calc
            
            print(f"   Parameters: {param_name}")
            print(f"   Cutoff: {params['mode']}")
            print(f"   Mixing: β={params['mixer']['beta']}")
            print(f"   Smearing: {params['occupations']['width']} eV")
            
            # Attempt single point calculation
            start_time = time.time()
            energy = test_atoms.get_potential_energy()
            elapsed = time.time() - start_time
            
            print(f"   ✅ SUCCESS: {energy:.4f} eV ({elapsed:.1f}s)")
            return True, energy, "Converged"
            
        except Exception as e:
            elapsed = time.time() - start_time
            error_msg = str(e)
            print(f"   ❌ FAILED: {error_msg[:50]}... ({elapsed:.1f}s)")
            
            # Analyze the log file
            log_analysis = self.analyze_convergence_log(params['txt'])
            return False, 0.0, error_msg
    
    def progressive_parameter_test(self, atoms: Atoms, name: str) -> dict:
        """Test progressively aggressive parameters until one works"""
        
        print(f"\n{'='*60}")
        print(f"🔧 PROGRESSIVE PARAMETER TEST: {name}")
        print(f"{'='*60}")
        print(f"Atoms: {len(atoms)}")
        print(f"Elements: {set(atoms.get_chemical_symbols())}")
        
        results = {
            'system': name,
            'atoms': len(atoms),
            'working_params': None,
            'energy': None,
            'all_results': {}
        }
        
        # Test each parameter set
        for param_name in ['ultra_conservative', 'conservative', 'standard', 'aggressive']:
            success, energy, message = self.test_parameter_set(atoms, param_name)
            
            results['all_results'][param_name] = {
                'success': success,
                'energy': energy,
                'message': message
            }
            
            if success:
                results['working_params'] = param_name
                results['energy'] = energy
                print(f"   🎉 Found working parameters: {param_name}")
                break
        
        return results
    
    def create_robust_calculator(self, atoms: Atoms, label: str) -> GPAW:
        """Create calculator with parameters that work for this system"""
        
        print(f"\n🔧 Creating robust calculator for {label}...")
        
        # Test what works for this system
        test_results = self.progressive_parameter_test(atoms, label)
        
        if test_results['working_params']:
            param_name = test_results['working_params']
            params = self.parameter_sets[param_name].copy()
            params['txt'] = f'{label}_robust.txt'
            
            print(f"   ✅ Using {param_name} parameters")
            return GPAW(**params)
        else:
            print(f"   ❌ No parameters worked - using ultra-conservative fallback")
            params = self.parameter_sets['ultra_conservative'].copy()
            params['txt'] = f'{label}_fallback.txt'
            return GPAW(**params)
    
    def diagnose_surface_issues(self) -> dict:
        """Diagnose specific issues with surface calculations"""
        
        print(f"\n🔍 DIAGNOSING SURFACE CALCULATION ISSUES")
        print("="*60)
        
        diagnostics = {
            'cluster_files': [],
            'surface_files': [],
            'recommendations': []
        }
        
        # Check for cluster files
        cluster_files = [
            'cluster_si4o6h6.xyz',
            'cluster_surface_patch_small.xyz',
            'test_surface_cluster_dft.xyz'
        ]
        
        for f in cluster_files:
            if os.path.exists(f):
                diagnostics['cluster_files'].append(f)
        
        # Check for surface files  
        surface_files = [
            'small_si100_2x2_4L_OH.xyz',
            'small_si100_2x2_6L_OH.xyz'
        ]
        
        for f in surface_files:
            if os.path.exists(f):
                diagnostics['surface_files'].append(f)
        
        # Generate recommendations
        if not diagnostics['cluster_files']:
            diagnostics['recommendations'].append("Run cluster_models.py to create small test systems")
        
        if not diagnostics['surface_files']:
            diagnostics['recommendations'].append("Run create_small_surfaces.py to create DFT-ready surfaces")
        
        diagnostics['recommendations'].extend([
            "Surface systems are harder to converge than molecules",
            "Try ultra-conservative parameters first",
            "Use smearing for metallic/surface systems",
            "Increase maxiter to 500+ for difficult systems"
        ])
        
        return diagnostics
    
    def fix_surface_convergence(self):
        """Attempt to fix surface convergence issues"""
        
        print(f"\n🔧 SURFACE CONVERGENCE TROUBLESHOOTING")
        print("="*60)
        
        # Diagnose issues
        diagnostics = self.diagnose_surface_issues()
        
        print(f"Available systems:")
        for f in diagnostics['cluster_files'] + diagnostics['surface_files']:
            print(f"  📄 {f}")
        
        print(f"\nRecommendations:")
        for rec in diagnostics['recommendations']:
            print(f"  • {rec}")
        
        # Test available systems
        test_files = diagnostics['cluster_files'][:2]  # Test first 2 clusters
        
        for test_file in test_files:
            if os.path.exists(test_file):
                print(f"\n{'─'*50}")
                print(f"Testing: {test_file}")
                print("─"*50)
                
                try:
                    atoms = read(test_file)
                    results = self.progressive_parameter_test(atoms, test_file)
                    
                    if results['working_params']:
                        print(f"✅ {test_file} works with {results['working_params']}")
                        
                        # Save working structure
                        working_file = test_file.replace('.xyz', '_working.xyz')
                        write(working_file, atoms)
                        
                    else:
                        print(f"❌ {test_file} failed with all parameter sets")
                        
                except Exception as e:
                    print(f"❌ Error testing {test_file}: {e}")

def main():
    """Run DFT troubleshooting"""
    
    troubleshooter = DFTTroubleshooter()
    
    print("🔧 DFT Troubleshooting for Surface Convergence")
    print("Addresses SCF convergence failures in surface systems")
    print("="*70)
    
    # Run comprehensive troubleshooting
    troubleshooter.fix_surface_convergence()
    
    print(f"\n💡 Next steps:")
    print(f"  1. Use working parameters identified above")
    print(f"  2. Start with ultra-conservative settings")
    print(f"  3. Increase smearing width for difficult systems")
    print(f"  4. Use structures marked as '_working.xyz'")
    
    print(f"\n🎯 For production runs:")
    print(f"  • Use troubleshooter.create_robust_calculator()")
    print(f"  • Test small systems before large ones")
    print(f"  • Monitor convergence in log files")
    
    return 0

if __name__ == "__main__":
    exit(main())