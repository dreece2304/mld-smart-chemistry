#!/usr/bin/env python3
"""
DFT Parameter Optimizer
Automatically optimizes DFT parameters based on system characteristics

Features:
- System size-based parameter selection
- Convergence difficulty detection
- Automatic parameter adjustment
- Smart restart logic
- Material-specific optimizations
"""

import numpy as np
from ase import Atoms
from ase.io import read, write
import time
import os
import json
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Any

try:
    from gpaw import GPAW, PW
    GPAW_AVAILABLE = True
except ImportError:
    GPAW_AVAILABLE = False

class DFTParameterOptimizer:
    """Automatically optimize DFT parameters for different systems"""
    
    def __init__(self):
        """Initialize with parameter rules and optimization strategies"""
        
        # System size categories (atoms)
        self.SIZE_CATEGORIES = {
            'tiny': (1, 20),        # Molecules
            'small': (21, 50),      # Small clusters
            'medium': (51, 150),    # Medium surfaces
            'large': (151, 300),    # Large surfaces
            'huge': (301, 1000)     # Very large systems
        }
        
        # Base parameter sets for different system sizes
        self.BASE_PARAMETERS = {
            'tiny': {
                'mode': PW(400),           # High accuracy for molecules
                'kpts': (1, 1, 1),
                'convergence': {'energy': 1e-6, 'density': 1e-5, 'eigenstates': 1e-5},
                'mixer': {'beta': 0.25},
                'maxiter': 200
            },
            'small': {
                'mode': PW(350),
                'kpts': (1, 1, 1),
                'convergence': {'energy': 1e-5, 'density': 1e-4, 'eigenstates': 1e-4},
                'mixer': {'beta': 0.2},
                'maxiter': 250
            },
            'medium': {
                'mode': PW(300),
                'kpts': (1, 1, 1),
                'convergence': {'energy': 1e-5, 'density': 1e-4, 'eigenstates': 1e-4},
                'mixer': {'beta': 0.15},
                'maxiter': 300
            },
            'large': {
                'mode': PW(250),          # Lower cutoff for speed
                'kpts': (1, 1, 1),
                'convergence': {'energy': 2e-5, 'density': 2e-4, 'eigenstates': 2e-4},
                'mixer': {'beta': 0.1},
                'maxiter': 400
            },
            'huge': {
                'mode': PW(200),          # Much lower cutoff
                'kpts': (1, 1, 1),
                'convergence': {'energy': 5e-5, 'density': 5e-4, 'eigenstates': 5e-4},
                'mixer': {'beta': 0.05},
                'maxiter': 500
            }
        }
        
        # Material-specific adjustments
        self.MATERIAL_ADJUSTMENTS = {
            'silicon': {
                'mixer_beta_factor': 0.8,     # Si needs gentler mixing
                'convergence_factor': 1.5,    # Looser convergence ok
            },
            'metal': {
                'mixer_beta_factor': 0.6,     # Metals harder to converge
                'convergence_factor': 2.0,
            },
            'organic': {
                'mixer_beta_factor': 1.2,     # Organics converge easier
                'convergence_factor': 0.8,
            }
        }
        
        # Convergence difficulty indicators
        self.DIFFICULTY_INDICATORS = {
            'oscillating_energy': 5,      # Energy oscillates for N steps
            'slow_density_conv': 1e-2,    # Density convergence slower than this
            'high_iteration_count': 50,   # More than N iterations without convergence
        }
        
        print("🎯 DFT Parameter Optimizer")
        print("   Automatically optimizes parameters based on system characteristics")
        
    def classify_system(self, atoms: Atoms) -> str:
        """Classify system by size"""
        n_atoms = len(atoms)
        
        for category, (min_atoms, max_atoms) in self.SIZE_CATEGORIES.items():
            if min_atoms <= n_atoms <= max_atoms:
                return category
        
        return 'huge'  # Default for very large systems
    
    def detect_material_type(self, atoms: Atoms) -> str:
        """Detect material type from atomic composition"""
        symbols = atoms.get_chemical_symbols()
        unique_symbols = set(symbols)
        
        # Silicon-based systems
        if 'Si' in unique_symbols:
            return 'silicon'
        
        # Metal systems
        metals = {'Al', 'Fe', 'Cu', 'Au', 'Ag', 'Pt', 'Pd', 'Ni', 'Ti', 'Zr'}
        if any(symbol in metals for symbol in unique_symbols):
            return 'metal'
        
        # Organic systems
        if 'C' in unique_symbols and 'H' in unique_symbols:
            return 'organic'
        
        return 'general'
    
    def get_optimized_parameters(self, atoms: Atoms, calculation_type: str = 'optimization') -> Dict[str, Any]:
        """
        Get optimized DFT parameters for a system
        
        Args:
            atoms: ASE Atoms object
            calculation_type: 'optimization', 'single_point', 'frequency'
        
        Returns:
            Optimized parameter dictionary
        """
        # Classify system
        size_category = self.classify_system(atoms)
        material_type = self.detect_material_type(atoms)
        n_atoms = len(atoms)
        
        print(f"\n🎯 Parameter Optimization")
        print(f"   System size: {n_atoms} atoms ({size_category})")
        print(f"   Material type: {material_type}")
        print(f"   Calculation: {calculation_type}")
        
        # Get base parameters
        params = self.BASE_PARAMETERS[size_category].copy()
        
        # Apply material-specific adjustments
        if material_type in self.MATERIAL_ADJUSTMENTS:
            adjustments = self.MATERIAL_ADJUSTMENTS[material_type]
            
            # Adjust mixer beta
            if 'mixer_beta_factor' in adjustments:
                params['mixer']['beta'] *= adjustments['mixer_beta_factor']
            
            # Adjust convergence criteria
            if 'convergence_factor' in adjustments:
                factor = adjustments['convergence_factor']
                for key in params['convergence']:
                    params['convergence'][key] *= factor
        
        # Calculation-type specific adjustments
        if calculation_type == 'single_point':
            # Tighter convergence for single point calculations
            for key in params['convergence']:
                params['convergence'][key] *= 0.5
        elif calculation_type == 'frequency':
            # Much tighter convergence for frequencies
            for key in params['convergence']:
                params['convergence'][key] *= 0.1
        
        # Add common parameters
        params.update({
            'xc': 'PBE',
            'symmetry': 'off',  # Safer for surfaces
            'nbands': 'nao'
        })
        
        print(f"   Cutoff: {params['mode'].ecut} eV")
        print(f"   Mixer β: {params['mixer']['beta']:.3f}")
        print(f"   Max iter: {params['maxiter']}")
        
        return params
    
    def analyze_convergence_progress(self, dft_output_file: str) -> Dict[str, Any]:
        """
        Analyze DFT convergence progress from output file
        
        Returns:
            Dictionary with convergence analysis
        """
        if not os.path.exists(dft_output_file):
            return {'status': 'no_file'}
        
        analysis = {
            'status': 'unknown',
            'iterations': 0,
            'energy_history': [],
            'density_conv_history': [],
            'eigst_conv_history': [],
            'is_oscillating': False,
            'is_slow_converging': False,
            'needs_parameter_adjustment': False
        }
        
        try:
            with open(dft_output_file, 'r') as f:
                lines = f.readlines()
            
            # Parse SCF iterations
            for line in lines:
                # Look for GPAW iteration lines: "iter:   1 21:21:08  -660.451065   -1.71  -0.82"
                if line.startswith('iter:'):
                    parts = line.split()
                    if len(parts) >= 6:
                        iter_num = int(parts[1])
                        energy = float(parts[3])
                        
                        analysis['iterations'] = iter_num
                        analysis['energy_history'].append(energy)
                        
                        # Parse convergence indicators if present
                        if len(parts) >= 5:
                            eigst_conv = parts[4].replace('c', '')  # Remove 'c' if converged
                            analysis['eigst_conv_history'].append(float(eigst_conv))
                        
                        if len(parts) >= 6:
                            dens_conv = parts[5].replace('c', '')
                            analysis['density_conv_history'].append(float(dens_conv))
            
            # Analyze patterns
            if len(analysis['energy_history']) >= 5:
                # Check for energy oscillations
                recent_energies = analysis['energy_history'][-5:]
                energy_range = max(recent_energies) - min(recent_energies)
                mean_energy = np.mean(recent_energies)
                if abs(energy_range) > 0.01 * abs(mean_energy):  # 1% oscillation
                    analysis['is_oscillating'] = True
            
            # Check for slow density convergence
            if analysis['density_conv_history']:
                latest_dens_conv = analysis['density_conv_history'][-1]
                if latest_dens_conv > self.DIFFICULTY_INDICATORS['slow_density_conv']:
                    analysis['is_slow_converging'] = True
            
            # Check iteration count
            if analysis['iterations'] > self.DIFFICULTY_INDICATORS['high_iteration_count']:
                analysis['is_slow_converging'] = True
            
            # Determine if parameter adjustment needed
            analysis['needs_parameter_adjustment'] = (
                analysis['is_oscillating'] or 
                analysis['is_slow_converging']
            )
            
            # Determine status
            if 'converged' in ''.join(lines[-10:]).lower():
                analysis['status'] = 'converged'
            elif analysis['needs_parameter_adjustment']:
                analysis['status'] = 'needs_adjustment'
            else:
                analysis['status'] = 'progressing'
                
        except Exception as e:
            analysis['error'] = str(e)
        
        return analysis
    
    def suggest_parameter_improvements(self, atoms: Atoms, convergence_analysis: Dict) -> Dict[str, Any]:
        """
        Suggest parameter improvements based on convergence analysis
        """
        current_params = self.get_optimized_parameters(atoms)
        improved_params = current_params.copy()
        
        suggestions = []
        
        if convergence_analysis['is_oscillating']:
            # Reduce mixer beta for oscillating calculations
            improved_params['mixer']['beta'] *= 0.5
            suggestions.append("Reduced mixer β for oscillation control")
            
            # Increase maxiter
            improved_params['maxiter'] = int(improved_params['maxiter'] * 1.5)
            suggestions.append("Increased max iterations")
        
        if convergence_analysis['is_slow_converging']:
            # Loosen convergence criteria slightly
            for key in improved_params['convergence']:
                improved_params['convergence'][key] *= 2.0
            suggestions.append("Loosened convergence criteria")
            
            # Reduce cutoff slightly for speed
            if hasattr(improved_params['mode'], 'ecut'):
                new_cutoff = max(200, improved_params['mode'].ecut - 50)
                improved_params['mode'] = PW(new_cutoff)
                suggestions.append(f"Reduced cutoff to {new_cutoff} eV")
        
        return {
            'improved_parameters': improved_params,
            'suggestions': suggestions,
            'original_parameters': current_params
        }
    
    def create_optimized_calculator(self, atoms: Atoms, label: str, 
                                  calculation_type: str = 'optimization') -> Optional[GPAW]:
        """
        Create an optimized GPAW calculator
        """
        if not GPAW_AVAILABLE:
            print("❌ GPAW not available")
            return None
        
        params = self.get_optimized_parameters(atoms, calculation_type)
        params['txt'] = f'{label}_dft.txt'
        
        try:
            calc = GPAW(**params)
            print(f"   ✅ Optimized calculator created for {label}")
            return calc
        except Exception as e:
            print(f"   ❌ Calculator creation failed: {e}")
            return None
    
    def save_optimization_report(self, atoms: Atoms, label: str, 
                               params: Dict, analysis: Dict = None):
        """Save parameter optimization report"""
        
        report = {
            'timestamp': datetime.now().isoformat(),
            'system_info': {
                'n_atoms': len(atoms),
                'size_category': self.classify_system(atoms),
                'material_type': self.detect_material_type(atoms),
                'chemical_formula': atoms.get_chemical_formula()
            },
            'optimized_parameters': params,
            'convergence_analysis': analysis or {}
        }
        
        report_file = f'{label}_param_optimization.json'
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2, default=str)
        
        print(f"   📊 Optimization report saved: {report_file}")

def main():
    """Example usage"""
    optimizer = DFTParameterOptimizer()
    
    # Example: optimize parameters for current surface
    surface_files = [
        'small_si100_2x2_4L.xyz',
        'small_si100_2x2_4L_OH.xyz'
    ]
    
    for surf_file in surface_files:
        if os.path.exists(surf_file):
            print(f"\n🔍 Analyzing: {surf_file}")
            atoms = read(surf_file)
            
            # Get optimized parameters
            params = optimizer.get_optimized_parameters(atoms, 'optimization')
            
            # Check if there's an ongoing calculation to analyze
            dft_file = f"{os.path.splitext(surf_file)[0]}_surf_opt.txt"
            if os.path.exists(dft_file):
                analysis = optimizer.analyze_convergence_progress(dft_file)
                print(f"   Convergence status: {analysis['status']}")
                
                if analysis['needs_parameter_adjustment']:
                    improvements = optimizer.suggest_parameter_improvements(atoms, analysis)
                    print("   💡 Suggested improvements:")
                    for suggestion in improvements['suggestions']:
                        print(f"      - {suggestion}")
            
            # Save report
            optimizer.save_optimization_report(atoms, os.path.splitext(surf_file)[0], params)

if __name__ == "__main__":
    main()