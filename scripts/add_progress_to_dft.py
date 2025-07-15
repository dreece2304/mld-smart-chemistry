#!/usr/bin/env python3
"""
Progress Bar Integration for DFT Scripts
Adds real-time progress tracking to existing DFT calculations

This module can be imported by other scripts to add progress bars
"""

import sys
import time
import numpy as np
from datetime import datetime, timedelta
from typing import Callable, Optional

# Try importing tqdm
try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False

class DFTProgress:
    """Add progress tracking to DFT calculations"""
    
    def __init__(self, use_color: bool = True):
        self.use_color = use_color
        self.colors = {
            'green': '\033[92m',
            'yellow': '\033[93m',
            'red': '\033[91m',
            'blue': '\033[94m',
            'bold': '\033[1m',
            'reset': '\033[0m'
        } if use_color else {k: '' for k in ['green', 'yellow', 'red', 'blue', 'bold', 'reset']}
    
    def print_status(self, message: str, status: str = 'info'):
        """Print colored status message"""
        colors = {
            'info': self.colors['blue'],
            'success': self.colors['green'],
            'warning': self.colors['yellow'],
            'error': self.colors['red']
        }
        color = colors.get(status, '')
        print(f"{color}{message}{self.colors['reset']}")
    
    def create_scf_monitor(self, max_iter: int = 300):
        """Create SCF convergence monitor"""
        
        class SCFMonitor:
            def __init__(self, max_iter: int):
                self.max_iter = max_iter
                self.iteration = 0
                self.start_time = time.time()
                self.energies = []
                self.density_errors = []
                
                # Progress bar
                if TQDM_AVAILABLE:
                    self.pbar = tqdm(total=max_iter, desc="SCF Iterations", 
                                   bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]')
                else:
                    self.pbar = None
                    print(f"SCF Progress: 0/{max_iter}")
            
            def update(self, energy: Optional[float] = None, 
                      density_error: Optional[float] = None):
                """Update SCF progress"""
                self.iteration += 1
                
                if energy is not None:
                    self.energies.append(energy)
                if density_error is not None:
                    self.density_errors.append(density_error)
                
                # Update progress bar
                if self.pbar:
                    self.pbar.update(1)
                    
                    # Add convergence info
                    if self.energies and len(self.energies) > 1:
                        de = abs(self.energies[-1] - self.energies[-2])
                        self.pbar.set_postfix({
                            'ΔE': f'{de:.2e}',
                            'ρ_err': f'{density_error:.2e}' if density_error else '---'
                        })
                else:
                    # Simple text progress
                    if self.iteration % 10 == 0 or self.iteration == 1:
                        elapsed = time.time() - self.start_time
                        rate = self.iteration / elapsed if elapsed > 0 else 0
                        
                        if self.energies and len(self.energies) > 1:
                            de = abs(self.energies[-1] - self.energies[-2])
                            print(f"SCF Progress: {self.iteration}/{self.max_iter} "
                                  f"(ΔE: {de:.2e}, {rate:.1f} iter/s)")
                        else:
                            print(f"SCF Progress: {self.iteration}/{self.max_iter} "
                                  f"({rate:.1f} iter/s)")
            
            def close(self):
                """Close progress bar"""
                if self.pbar:
                    self.pbar.close()
                
                # Summary
                elapsed = time.time() - self.start_time
                print(f"\nSCF Summary:")
                print(f"  Iterations: {self.iteration}")
                print(f"  Time: {elapsed:.1f} seconds")
                
                if self.energies and len(self.energies) > 1:
                    final_de = abs(self.energies[-1] - self.energies[-2])
                    print(f"  Final ΔE: {final_de:.2e} eV")
        
        return SCFMonitor(max_iter)
    
    def create_optimization_monitor(self, max_steps: int = 100, fmax_target: float = 0.05):
        """Create geometry optimization monitor"""
        
        class OptMonitor:
            def __init__(self, max_steps: int, fmax_target: float):
                self.max_steps = max_steps
                self.fmax_target = fmax_target
                self.step = 0
                self.start_time = time.time()
                self.energies = []
                self.forces = []
                
                # Progress bar
                if TQDM_AVAILABLE:
                    self.pbar = tqdm(total=max_steps, desc="Optimization Steps",
                                   bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [Fmax: {postfix}]')
                else:
                    self.pbar = None
                    print(f"Optimization Progress: 0/{max_steps} (Target Fmax: {fmax_target} eV/Å)")
            
            def update(self, energy: float, fmax: float):
                """Update optimization progress"""
                self.step += 1
                self.energies.append(energy)
                self.forces.append(fmax)
                
                # Check if converged
                converged = fmax < self.fmax_target
                
                # Update progress
                if self.pbar:
                    self.pbar.update(1)
                    
                    # Color code fmax based on convergence
                    if converged:
                        fmax_str = f'\033[92m{fmax:.4f}\033[0m'  # Green
                    elif fmax < self.fmax_target * 2:
                        fmax_str = f'\033[93m{fmax:.4f}\033[0m'  # Yellow
                    else:
                        fmax_str = f'\033[91m{fmax:.4f}\033[0m'  # Red
                    
                    self.pbar.set_postfix_str(f'{fmax_str} eV/Å')
                    
                    if converged:
                        self.pbar.set_description("Optimization CONVERGED")
                else:
                    # Simple text progress
                    if self.step % 5 == 0 or self.step == 1 or converged:
                        elapsed = time.time() - self.start_time
                        
                        # Energy change
                        if len(self.energies) > 1:
                            de = self.energies[-1] - self.energies[-2]
                            print(f"Step {self.step}: E = {energy:.4f} eV (ΔE = {de:+.4f}), "
                                  f"Fmax = {fmax:.4f} eV/Å {'✓' if converged else ''}")
                        else:
                            print(f"Step {self.step}: E = {energy:.4f} eV, "
                                  f"Fmax = {fmax:.4f} eV/Å")
                
                return converged
            
            def close(self):
                """Close progress bar and show summary"""
                if self.pbar:
                    self.pbar.close()
                
                # Summary
                elapsed = time.time() - self.start_time
                final_fmax = self.forces[-1] if self.forces else 999.9
                converged = final_fmax < self.fmax_target
                
                print(f"\nOptimization Summary:")
                print(f"  Steps: {self.step}")
                print(f"  Converged: {'✅ Yes' if converged else '❌ No'}")
                print(f"  Final Fmax: {final_fmax:.4f} eV/Å")
                print(f"  Time: {elapsed/60:.1f} minutes")
                
                # Show convergence trend
                if len(self.forces) > 3:
                    self.plot_convergence()
            
            def plot_convergence(self):
                """Simple ASCII convergence plot"""
                if len(self.forces) < 2:
                    return
                
                print("\n📈 Force Convergence:")
                
                # Get last 20 points or all if less
                plot_forces = self.forces[-20:]
                max_f = max(plot_forces)
                
                for i, f in enumerate(plot_forces):
                    # Normalize to 30 character width
                    bar_length = int(30 * (f / max_f)) if max_f > 0 else 0
                    bar = "█" * bar_length
                    
                    # Color based on convergence
                    if f < self.fmax_target:
                        color = '\033[92m'  # Green
                    elif f < self.fmax_target * 2:
                        color = '\033[93m'  # Yellow  
                    else:
                        color = '\033[91m'  # Red
                    
                    print(f"  Step {len(self.forces)-len(plot_forces)+i+1:3d}: "
                          f"{color}{bar:<30}\033[0m {f:.4f}")
                
                # Show target line
                print(f"  {'─'*10} Target: {self.fmax_target:.4f} eV/Å")
        
        return OptMonitor(max_steps, fmax_target)
    
    def create_mld_cycle_monitor(self, n_cycles: int):
        """Create MLD cycle progress monitor"""
        
        class MLDMonitor:
            def __init__(self, n_cycles: int):
                self.n_cycles = n_cycles
                self.current_cycle = 0
                self.start_time = time.time()
                self.cycle_times = []
                self.tma_reactions = []
                self.h2o_reactions = []
                
                print(f"\n🔄 MLD Process: {n_cycles} cycles planned")
                print("─" * 60)
            
            def start_cycle(self, cycle_num: int):
                """Start a new cycle"""
                self.current_cycle = cycle_num
                self.cycle_start = time.time()
                
                print(f"\n{'='*60}")
                print(f"📍 MLD Cycle {cycle_num}/{self.n_cycles}")
                print(f"{'='*60}")
                
                # Estimate remaining time
                if self.cycle_times:
                    avg_cycle_time = np.mean(self.cycle_times)
                    remaining_cycles = self.n_cycles - cycle_num + 1
                    eta = datetime.now() + timedelta(seconds=avg_cycle_time * remaining_cycles)
                    print(f"⏱️  Estimated completion: {eta.strftime('%H:%M:%S')}")
            
            def update_tma_pulse(self, n_reactions: int):
                """Update TMA pulse results"""
                self.tma_reactions.append(n_reactions)
                print(f"✅ TMA pulse complete: {n_reactions} reactions")
            
            def update_h2o_pulse(self, n_reactions: int):
                """Update H2O pulse results"""
                self.h2o_reactions.append(n_reactions)
                print(f"✅ H2O pulse complete: {n_reactions} reactions")
            
            def end_cycle(self):
                """End current cycle"""
                cycle_time = time.time() - self.cycle_start
                self.cycle_times.append(cycle_time)
                
                print(f"\n⏱️  Cycle {self.current_cycle} time: {cycle_time/60:.1f} minutes")
                
                # Progress bar for overall process
                progress = self.current_cycle / self.n_cycles
                bar_length = 40
                filled = int(bar_length * progress)
                bar = "█" * filled + "░" * (bar_length - filled)
                
                print(f"\nOverall Progress: [{bar}] {progress*100:.0f}%")
            
            def close(self):
                """Final summary"""
                total_time = time.time() - self.start_time
                
                print(f"\n{'='*60}")
                print(f"✅ MLD PROCESS COMPLETE")
                print(f"{'='*60}")
                print(f"Total cycles: {len(self.cycle_times)}")
                print(f"Total time: {total_time/3600:.2f} hours")
                
                if self.cycle_times:
                    print(f"Average cycle time: {np.mean(self.cycle_times)/60:.1f} minutes")
                
                if self.tma_reactions:
                    print(f"Total TMA reactions: {sum(self.tma_reactions)}")
                    print(f"Average per cycle: {np.mean(self.tma_reactions):.1f}")
        
        return MLDMonitor(n_cycles)

# Example usage functions
def example_scf_progress():
    """Example of SCF progress monitoring"""
    progress = DFTProgress()
    scf_monitor = progress.create_scf_monitor(max_iter=100)
    
    # Simulate SCF iterations
    energy = -10.0
    for i in range(75):
        energy += np.random.normal(0, 0.1) * np.exp(-i/20)
        density_error = 0.1 * np.exp(-i/15)
        
        scf_monitor.update(energy=energy, density_error=density_error)
        time.sleep(0.01)
    
    scf_monitor.close()

def example_optimization_progress():
    """Example of optimization progress monitoring"""
    progress = DFTProgress()
    opt_monitor = progress.create_optimization_monitor(max_steps=50, fmax_target=0.05)
    
    # Simulate optimization
    energy = -10.0
    fmax = 0.5
    
    for i in range(35):
        energy -= 0.1 * np.exp(-i/10)
        fmax = 0.5 * np.exp(-i/8) + 0.02
        
        converged = opt_monitor.update(energy, fmax)
        time.sleep(0.05)
        
        if converged:
            break
    
    opt_monitor.close()

def main():
    """Demo progress tracking"""
    print("🎯 DFT Progress Tracking Demo")
    print("="*60)
    
    print("\n1. SCF Convergence Progress:")
    example_scf_progress()
    
    print("\n2. Geometry Optimization Progress:")
    example_optimization_progress()
    
    print("\n✅ Progress tracking demo complete!")
    print("\n💡 To use in your scripts:")
    print("   from add_progress_to_dft import DFTProgress")
    print("   progress = DFTProgress()")
    print("   scf_monitor = progress.create_scf_monitor()")
    
    return 0

if __name__ == "__main__":
    exit(main())