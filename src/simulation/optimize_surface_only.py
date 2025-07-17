#!/usr/bin/env python3
"""
Standalone Surface Optimizer for DFT
Pre-optimize surface structures before MLD simulations

This script:
- Optimizes only the surface structure
- Saves the optimized structure and energy
- Can be run independently before MLD simulations
"""

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS, LBFGS
from ase.constraints import FixAtoms
import time
import os
import json
from datetime import datetime
from typing import Optional, Dict
import threading
from tqdm import tqdm
import sys

try:
    from gpaw import GPAW, PW
    GPAW_AVAILABLE = True
except ImportError:
    GPAW_AVAILABLE = False
    print("❌ GPAW not available! Install with: conda install -c conda-forge gpaw")

class DFTProgressMonitor:
    """Monitor DFT calculation progress"""
    
    def __init__(self, max_iter=300):
        self.max_iter = max_iter
        self.current_iter = 0
        self.pbar = None
        self.converged = False
        self.energy_history = []
        self.start_time = time.time()
        
    def start(self, desc="DFT SCF"):
        """Start progress bar"""
        self.pbar = tqdm(total=self.max_iter, desc=desc, unit="iter", 
                        bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]')
        
    def update(self, calc):
        """Update progress from calculator"""
        try:
            # Get iteration info from calculator
            if hasattr(calc, 'iter'):
                new_iter = calc.iter
                if new_iter > self.current_iter:
                    self.pbar.update(new_iter - self.current_iter)
                    self.current_iter = new_iter
                    
            # Get convergence info
            if hasattr(calc, 'hamiltonian'):
                ham = calc.hamiltonian
                if hasattr(ham, 'e_entropy'):
                    energy = ham.e_entropy
                    self.energy_history.append(energy)
                    
                    # Show energy change
                    if len(self.energy_history) > 1:
                        de = abs(self.energy_history[-1] - self.energy_history[-2])
                        self.pbar.set_postfix({'ΔE': f'{de:.2e} eV'})
                        
        except Exception:
            pass  # Silently ignore monitoring errors
            
    def close(self):
        """Close progress bar"""
        if self.pbar:
            self.pbar.close()

class SurfaceOptimizer:
    """Standalone surface optimizer for DFT calculations"""
    
    def __init__(self):
        """Initialize with optimized parameters for surfaces"""
        
        # DFT parameters optimized for surface calculations
        self.dft_params = {
            'mode': PW(300),           # Moderate cutoff for surfaces
            'xc': 'PBE',              # Standard for surfaces
            'kpts': (1, 1, 1),        # Gamma point for slabs
            'symmetry': 'off',        # Important for surfaces
            'convergence': {
                'energy': 1e-5,
                'density': 1e-4,
                'eigenstates': 1e-4
            },
            'mixer': {'beta': 0.15},   # Slightly higher for surfaces
            'maxiter': 300,
            'nbands': 'nao'
        }
        
        # Optimization parameters
        self.opt_params = {
            'fmax': 0.05,              # Reasonable for surfaces
            'max_steps': 100,
            'optimizer': LBFGS
        }
        
        print("🔬 Surface Optimizer for DFT")
        print("   Optimizes surface structures independently")
        print("   Saves results for use in MLD simulations")
        
    def optimize_surface(self, surface_file: str, output_prefix: str = None, restart: bool = False) -> Dict:
        """
        Optimize a surface structure
        
        Args:
            surface_file: Path to surface structure file
            output_prefix: Prefix for output files (default: use input filename)
            restart: Whether to restart from previous calculation
            
        Returns:
            Dictionary with optimization results
        """
        
        if not os.path.exists(surface_file):
            print(f"❌ Surface file not found: {surface_file}")
            return {'success': False, 'error': 'File not found'}
        
        # Load surface
        print(f"\n📂 Loading surface: {surface_file}")
        surface = read(surface_file)
        n_atoms = len(surface)
        print(f"   Atoms: {n_atoms}")
        
        # Determine output prefix
        if output_prefix is None:
            output_prefix = os.path.splitext(os.path.basename(surface_file))[0]
        
        # Check for restart
        gpw_file = f'{output_prefix}_surf_opt.gpw'
        partial_structure = f'{output_prefix}_surf_partial.xyz'
        optimized_structure = f'{output_prefix}_surf_optimized.xyz'
        
        if restart:
            # Try to load from previous calculation
            if os.path.exists(optimized_structure):
                print(f"   ✅ Found completed optimization: {optimized_structure}")
                surface = read(optimized_structure)
                print(f"   Loading completed structure...")
            elif os.path.exists(partial_structure):
                print(f"   📂 Found partial structure: {partial_structure}")
                surface = read(partial_structure)
                print(f"   Restarting from partial optimization...")
            elif os.path.exists(gpw_file):
                print(f"   📂 Found previous calculation: {gpw_file}")
                print(f"   Will restart from saved state...")
            else:
                print(f"   ⚠️  No previous calculation found, starting fresh...")
                restart = False
        
        # Check size
        if n_atoms > 200:
            print(f"   ⚠️  WARNING: {n_atoms} atoms is large for DFT!")
            response = input("   Continue anyway? (y/n): ")
            if response.lower() != 'y':
                return {'success': False, 'error': 'User cancelled due to size'}
        
        # Set up constraints (fix bottom layers)
        positions = surface.get_positions()
        z_coords = positions[:, 2]
        z_min = z_coords.min()
        z_cutoff = z_min + 3.0  # Fix atoms within 3Å of bottom
        
        fixed_indices = [i for i, z in enumerate(z_coords) if z < z_cutoff]
        if fixed_indices:
            constraint = FixAtoms(indices=fixed_indices)
            surface.set_constraint(constraint)
            print(f"   Fixed {len(fixed_indices)} bottom atoms")
        
        # Create calculator
        print("\n⚙️  Setting up DFT calculator...")
        calc_params = self.dft_params.copy()
        calc_params['txt'] = f'{output_prefix}_surf_opt.txt'
        
        # Enable verbose output if requested
        if hasattr(self, 'monitor_log') and self.monitor_log:
            calc_params['txt'] = '-'  # Output to stdout
        
        try:
            if restart and os.path.exists(gpw_file):
                # Restart from saved calculation
                print(f"   📂 Restarting from: {gpw_file}")
                calc = GPAW(gpw_file, txt=calc_params['txt'])
                surface.calc = calc
                print("   ✅ Calculator restarted from previous state")
            else:
                # Fresh calculation
                calc = GPAW(**calc_params)
                surface.calc = calc
                print("   ✅ Calculator ready")
        except Exception as e:
            print(f"   ❌ Calculator setup failed: {e}")
            return {'success': False, 'error': str(e)}
        
        # Run optimization
        print(f"\n🏃 Starting optimization...")
        print(f"   Target fmax: {self.opt_params['fmax']} eV/Å")
        print(f"   Max steps: {self.opt_params['max_steps']}")
        print(f"   Log file: {output_prefix}_surf_opt.log")
        
        opt_log = f'{output_prefix}_surf_opt.log'
        opt = LBFGS(surface, logfile=opt_log)
        
        start_time = time.time()
        
        # Progress tracking variables
        step_progress = None
        dft_monitor = DFTProgressMonitor(max_iter=self.dft_params['maxiter'])
        optimization_data = {
            'steps': [],
            'energies': [],
            'forces': [],
            'times': []
        }
        
        try:
            # Monitor DFT iterations during initial energy calculation
            print("   Calculating initial energy...")
            
            # Start monitoring thread for DFT progress
            monitor_active = True
            def monitor_dft():
                dft_monitor.start("Initial DFT")
                while monitor_active:
                    dft_monitor.update(calc)
                    time.sleep(0.1)
                dft_monitor.close()
                
            monitor_thread = threading.Thread(target=monitor_dft)
            monitor_thread.start()
            
            initial_energy = surface.get_potential_energy()
            monitor_active = False
            monitor_thread.join()
            
            print(f"   Initial energy: {initial_energy:.4f} eV")
            
            # Save initial state
            calc.write(gpw_file, mode='all')
            
            # Run optimization with comprehensive progress tracking
            print("\n📊 Starting optimization with live progress tracking...")
            
            # Create optimization progress bar
            opt_pbar = tqdm(total=self.opt_params['max_steps'], 
                           desc="Optimization", 
                           unit="steps",
                           bar_format='{desc}: {percentage:3.0f}%|{bar}| Step {n_fmt}/{total_fmt} [{elapsed}<{remaining}]')
            
            step_count = 0
            last_energy = initial_energy
            
            def optimization_callback():
                nonlocal step_count, last_energy, monitor_active
                step_count += 1
                
                # Update optimization progress
                opt_pbar.update(1)
                
                # Get current state
                energy = surface.get_potential_energy()
                forces = surface.get_forces()
                max_force = np.max(np.abs(forces[~surface.constraints[0].index]))
                
                # Calculate changes
                energy_change = energy - last_energy
                last_energy = energy
                
                # Update progress bar info
                opt_pbar.set_postfix({
                    'E': f'{energy:.3f} eV',
                    'ΔE': f'{energy_change:.2e}',
                    'Fmax': f'{max_force:.3f} eV/Å'
                })
                
                # Store optimization data
                optimization_data['steps'].append(step_count)
                optimization_data['energies'].append(float(energy))
                optimization_data['forces'].append(float(max_force))
                optimization_data['times'].append(time.time() - start_time)
                
                # Save progress every 5 steps
                if step_count % 5 == 0:
                    calc.write(gpw_file, mode='all')
                    write(partial_structure, surface)
                    
                    # Save optimization data
                    opt_data_file = f'{output_prefix}_optimization_progress.json'
                    with open(opt_data_file, 'w') as f:
                        json.dump(optimization_data, f, indent=2)
                
                # Check convergence
                if max_force < self.opt_params['fmax']:
                    opt_pbar.set_description("✅ Converged!")
                    return True
                    
            # Attach callback
            opt.attach(optimization_callback)
            
            # Custom optimization loop with DFT monitoring
            converged = False
            for i in range(self.opt_params['max_steps']):
                if converged:
                    break
                    
                # Monitor DFT iterations for this step
                dft_monitor = DFTProgressMonitor(max_iter=self.dft_params['maxiter'])
                monitor_active = True
                
                def monitor_step_dft():
                    dft_monitor.start(f"DFT Step {i+1}")
                    while monitor_active:
                        dft_monitor.update(calc)
                        time.sleep(0.1)
                    dft_monitor.close()
                
                monitor_thread = threading.Thread(target=monitor_step_dft)
                monitor_thread.start()
                
                # Run one optimization step
                converged = opt.step()
                
                # Stop DFT monitoring
                monitor_active = False
                monitor_thread.join()
                
                # Check forces
                forces = surface.get_forces()
                max_force = np.max(np.abs(forces[~surface.constraints[0].index]))
                if max_force < self.opt_params['fmax']:
                    converged = True
                    break
                    
            opt_pbar.close()
            
            # Get final results
            final_energy = surface.get_potential_energy()
            forces = surface.get_forces()
            max_force = np.max(np.abs(forces[~surface.constraints[0].index]))  # Max force on unfixed atoms
            elapsed = time.time() - start_time
            
            print(f"\n✅ Optimization complete!")
            print(f"   Initial energy: {initial_energy:.4f} eV")
            print(f"   Final energy: {final_energy:.4f} eV")
            print(f"   Energy change: {final_energy - initial_energy:.4f} eV")
            print(f"   Max force: {max_force:.4f} eV/Å")
            print(f"   Converged: {converged}")
            print(f"   Time: {elapsed/60:.1f} minutes")
            
            # Save optimized structure
            opt_structure_file = f'{output_prefix}_surf_optimized.xyz'
            write(opt_structure_file, surface)
            print(f"\n💾 Saved optimized structure: {opt_structure_file}")
            
            # Save results to JSON
            results = {
                'success': True,
                'surface_file': surface_file,
                'n_atoms': n_atoms,
                'initial_energy': float(initial_energy),
                'final_energy': float(final_energy),
                'energy_change': float(final_energy - initial_energy),
                'max_force': float(max_force),
                'converged': bool(converged),
                'optimization_time': elapsed,
                'fixed_atoms': len(fixed_indices) if fixed_indices else 0,
                'timestamp': datetime.now().isoformat()
            }
            
            results_file = f'{output_prefix}_surf_results.json'
            with open(results_file, 'w') as f:
                json.dump(results, f, indent=2)
            print(f"💾 Saved results: {results_file}")
            
            # Update pre_calculated_energies.json for dft_fast_mld.py
            energy_cache_file = "config/pre_calculated_energies.json"
            energy_cache = {}
            if os.path.exists(energy_cache_file):
                try:
                    with open(energy_cache_file, 'r') as f:
                        energy_cache = json.load(f)
                except:
                    pass
            
            energy_cache['surface'] = float(final_energy)
            energy_cache[f'{output_prefix}_surface'] = float(final_energy)
            
            with open(energy_cache_file, 'w') as f:
                json.dump(energy_cache, f, indent=2)
            print(f"💾 Updated energy cache: {energy_cache_file}")
            
            # Generate optimization plot
            self._generate_optimization_plot(optimization_data, output_prefix)
            
            return results
            
        except Exception as e:
            elapsed = time.time() - start_time
            print(f"\n❌ Optimization failed: {e}")
            print(f"   Time before failure: {elapsed/60:.1f} minutes")
            
            # Try to save partial results
            try:
                write(f'{output_prefix}_surf_partial.xyz', surface)
                print(f"   💾 Saved partial structure: {output_prefix}_surf_partial.xyz")
            except:
                pass
            
            return {
                'success': False,
                'error': str(e),
                'time': elapsed
            }
    
    def _generate_optimization_plot(self, opt_data: Dict, prefix: str):
        """Generate optimization progress plot"""
        try:
            import matplotlib.pyplot as plt
            
            if len(opt_data['steps']) < 2:
                return
                
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
            
            # Energy plot
            ax1.plot(opt_data['steps'], opt_data['energies'], 'b-', linewidth=2)
            ax1.set_ylabel('Energy (eV)', fontsize=12)
            ax1.set_title('Optimization Progress', fontsize=14)
            ax1.grid(True, alpha=0.3)
            
            # Force plot
            ax2.plot(opt_data['steps'], opt_data['forces'], 'r-', linewidth=2)
            ax2.axhline(y=self.opt_params['fmax'], color='g', linestyle='--', 
                       label=f'Target: {self.opt_params["fmax"]} eV/Å')
            ax2.set_ylabel('Max Force (eV/Å)', fontsize=12)
            ax2.set_xlabel('Optimization Step', fontsize=12)
            ax2.set_yscale('log')
            ax2.grid(True, alpha=0.3)
            ax2.legend()
            
            plt.tight_layout()
            plot_file = f'{prefix}_optimization_progress.png'
            plt.savefig(plot_file, dpi=150)
            plt.close()
            
            print(f"📈 Saved optimization plot: {plot_file}")
        except ImportError:
            pass  # Matplotlib not available

def main():
    """Main execution"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Optimize surface structures for DFT MLD simulations'
    )
    parser.add_argument(
        'surface_file',
        help='Path to surface structure file (xyz, traj, etc.)'
    )
    parser.add_argument(
        '--output-prefix',
        help='Prefix for output files (default: use input filename)'
    )
    parser.add_argument(
        '--fmax',
        type=float,
        default=0.05,
        help='Force convergence criterion in eV/Å (default: 0.05)'
    )
    parser.add_argument(
        '--steps',
        type=int,
        default=100,
        help='Maximum optimization steps (default: 100)'
    )
    parser.add_argument(
        '--restart',
        action='store_true',
        help='Restart from previous calculation if available'
    )
    parser.add_argument(
        '--monitor-log',
        action='store_true',
        help='Show detailed DFT output in real-time (useful for debugging)'
    )
    
    args = parser.parse_args()
    
    # Create optimizer
    optimizer = SurfaceOptimizer()
    
    # Override parameters if specified
    if args.fmax:
        optimizer.opt_params['fmax'] = args.fmax
    if args.steps:
        optimizer.opt_params['max_steps'] = args.steps
    if args.monitor_log:
        optimizer.monitor_log = True
    
    # Run optimization
    results = optimizer.optimize_surface(
        args.surface_file,
        args.output_prefix,
        args.restart
    )
    
    if results['success']:
        print("\n🎉 Surface optimization successful!")
        print(f"   Use '{results.get('output_prefix', 'surface')}_surf_optimized.xyz' in your MLD simulations")
    else:
        print("\n❌ Surface optimization failed")
        print(f"   Error: {results.get('error', 'Unknown error')}")

if __name__ == "__main__":
    main()