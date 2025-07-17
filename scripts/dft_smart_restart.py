#!/usr/bin/env python3
"""
DFT Smart Restart System
Automatically detects failed/stalled calculations and restarts with improved parameters

Features:
- Automatic failure detection
- Parameter adjustment based on failure type
- Checkpoint recovery
- Progress preservation
- Intelligent retry strategies
"""

import os
import time
import shutil
import json
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Tuple
from ase import Atoms
from ase.io import read, write

try:
    from gpaw import GPAW, PW
    GPAW_AVAILABLE = True
except ImportError:
    GPAW_AVAILABLE = False

from dft_parameter_optimizer import DFTParameterOptimizer
from dft_convergence_monitor import DFTConvergenceMonitor

class DFTSmartRestart:
    """Smart restart system for failed DFT calculations"""
    
    def __init__(self):
        """Initialize restart system"""
        
        self.optimizer = DFTParameterOptimizer()
        self.monitor = DFTConvergenceMonitor()
        
        # Failure detection criteria
        self.FAILURE_CRITERIA = {
            'file_age_hours': 2.0,          # File not updated in 2 hours
            'stall_iterations': 50,         # No progress for 50 iterations
            'energy_explosion': 1000.0,     # Energy > 1000 eV indicates failure
            'oscillation_cycles': 10,       # Oscillating for 10+ cycles
            'memory_error_keywords': ['out of memory', 'memory', 'malloc', 'killed'],
            'convergence_error_keywords': ['convergence', 'failed', 'error', 'exception']
        }
        
        # Restart strategies
        self.RESTART_STRATEGIES = {
            'stalled_scf': {
                'mixer_reduction': 0.5,     # Reduce mixer beta
                'cutoff_reduction': 50,     # Reduce cutoff by 50 eV
                'convergence_relaxation': 2.0  # Relax convergence by factor of 2
            },
            'memory_error': {
                'cutoff_reduction': 100,    # Aggressive cutoff reduction
                'convergence_relaxation': 5.0,
                'block_size_reduction': True
            },
            'oscillating': {
                'mixer_reduction': 0.3,     # Strong mixer reduction
                'density_mixing_change': True,
                'max_iter_increase': 1.5
            },
            'slow_convergence': {
                'cutoff_reduction': 25,
                'convergence_relaxation': 1.5,
                'mixer_adjustment': 0.8
            }
        }
        
        self.restart_history = []
        
        print("🔄 DFT Smart Restart System")
        print("   Automatic failure detection and intelligent restart")
    
    def detect_calculation_status(self, calculation_dir: str, 
                                calculation_name: str) -> Dict[str, any]:
        """
        Detect the status of a DFT calculation
        
        Returns:
            Dictionary with status information
        """
        
        status = {
            'calculation_name': calculation_name,
            'status': 'unknown',
            'failure_type': None,
            'needs_restart': False,
            'last_activity': None,
            'file_info': {},
            'convergence_info': {}
        }
        
        # Look for relevant files
        dft_file = os.path.join(calculation_dir, f"{calculation_name}_dft.txt")
        opt_log = os.path.join(calculation_dir, f"{calculation_name}_opt.log")
        structure_file = os.path.join(calculation_dir, f"{calculation_name}_optimized.xyz")
        
        status['file_info'] = {
            'dft_file': dft_file if os.path.exists(dft_file) else None,
            'opt_log': opt_log if os.path.exists(opt_log) else None,
            'structure_file': structure_file if os.path.exists(structure_file) else None
        }
        
        # Check if calculation completed successfully
        if os.path.exists(structure_file):
            status['status'] = 'completed'
            return status
        
        # Analyze main DFT file
        if os.path.exists(dft_file):
            status.update(self._analyze_dft_file(dft_file))
        else:
            status['status'] = 'not_started'
            return status
        
        # Check for specific failure patterns
        status['failure_type'] = self._detect_failure_type(status)
        
        # Determine if restart is needed
        status['needs_restart'] = self._should_restart(status)
        
        return status
    
    def _analyze_dft_file(self, dft_file: str) -> Dict:
        """Analyze DFT output file for status"""
        
        analysis = {
            'file_size': os.path.getsize(dft_file),
            'last_modified': datetime.fromtimestamp(os.path.getmtime(dft_file)),
            'age_hours': (datetime.now() - datetime.fromtimestamp(os.path.getmtime(dft_file))).total_seconds() / 3600
        }
        
        # Parse convergence information
        convergence_data = self.monitor.parse_dft_file(dft_file)
        analysis['convergence_info'] = convergence_data
        
        # Determine status from file content
        if convergence_data['is_converged']:
            analysis['status'] = 'converged'
        elif convergence_data['current_status'] == 'error':
            analysis['status'] = 'failed'
        elif analysis['age_hours'] > self.FAILURE_CRITERIA['file_age_hours']:
            analysis['status'] = 'stalled'
        elif convergence_data['scf_iterations']:
            analysis['status'] = 'running'
        else:
            analysis['status'] = 'starting'
        
        return analysis
    
    def _detect_failure_type(self, status: Dict) -> Optional[str]:
        """Detect specific type of failure"""
        
        convergence_info = status.get('convergence_info', {})
        scf_data = convergence_info.get('scf_iterations', [])
        
        # Check for stalled calculation
        if status.get('age_hours', 0) > self.FAILURE_CRITERIA['file_age_hours']:
            return 'stalled'
        
        # Check for too many iterations
        if scf_data and len(scf_data) > self.FAILURE_CRITERIA['stall_iterations']:
            return 'slow_convergence'
        
        # Check for energy explosion
        if scf_data:
            latest_energy = scf_data[-1].get('energy', 0)
            if abs(latest_energy) > self.FAILURE_CRITERIA['energy_explosion']:
                return 'energy_explosion'
        
        # Check for oscillations
        if len(scf_data) >= self.FAILURE_CRITERIA['oscillation_cycles']:
            recent_energies = [step['energy'] for step in scf_data[-self.FAILURE_CRITERIA['oscillation_cycles']:]]
            energy_range = max(recent_energies) - min(recent_energies)
            mean_energy = abs(sum(recent_energies) / len(recent_energies))
            if energy_range > 0.1 * mean_energy:  # 10% oscillation
                return 'oscillating'
        
        # Check for error messages in file
        dft_file = status['file_info'].get('dft_file')
        if dft_file and os.path.exists(dft_file):
            try:
                with open(dft_file, 'r') as f:
                    content = f.read().lower()
                
                for keyword in self.FAILURE_CRITERIA['memory_error_keywords']:
                    if keyword in content:
                        return 'memory_error'
                
                for keyword in self.FAILURE_CRITERIA['convergence_error_keywords']:
                    if keyword in content:
                        return 'convergence_error'
                        
            except:
                pass
        
        return None
    
    def _should_restart(self, status: Dict) -> bool:
        """Determine if calculation should be restarted"""
        
        # Don't restart if completed or running normally
        if status['status'] in ['completed', 'converged']:
            return False
        
        # Restart if we detected a failure
        if status['failure_type'] is not None:
            return True
        
        # Restart if stalled
        if status['status'] == 'stalled':
            return True
        
        return False
    
    def create_restart_parameters(self, original_atoms: Atoms, 
                                failure_type: str, 
                                restart_attempt: int = 1) -> Dict:
        """
        Create improved parameters for restart based on failure type
        """
        
        # Get base optimized parameters
        base_params = self.optimizer.get_optimized_parameters(original_atoms)
        restart_params = base_params.copy()
        
        # Apply failure-specific strategies
        if failure_type in self.RESTART_STRATEGIES:
            strategy = self.RESTART_STRATEGIES[failure_type]
            
            # Adjust mixer beta
            if 'mixer_reduction' in strategy:
                restart_params['mixer']['beta'] *= strategy['mixer_reduction']
            elif 'mixer_adjustment' in strategy:
                restart_params['mixer']['beta'] *= strategy['mixer_adjustment']
            
            # Adjust cutoff energy
            if 'cutoff_reduction' in strategy:
                current_cutoff = restart_params['mode'].ecut
                new_cutoff = max(150, current_cutoff - strategy['cutoff_reduction'])
                restart_params['mode'] = PW(new_cutoff)
            
            # Relax convergence criteria
            if 'convergence_relaxation' in strategy:
                factor = strategy['convergence_relaxation']
                for key in restart_params['convergence']:
                    restart_params['convergence'][key] *= factor
            
            # Increase max iterations
            if 'max_iter_increase' in strategy:
                restart_params['maxiter'] = int(restart_params['maxiter'] * strategy['max_iter_increase'])
        
        # Apply progressive adjustments for multiple restart attempts
        if restart_attempt > 1:
            # Get more aggressive with each restart
            restart_params['mixer']['beta'] *= (0.8 ** restart_attempt)
            
            # Reduce cutoff further
            current_cutoff = restart_params['mode'].ecut
            new_cutoff = max(100, current_cutoff - (25 * restart_attempt))
            restart_params['mode'] = PW(new_cutoff)
            
            # Relax convergence more
            factor = 1.5 ** restart_attempt
            for key in restart_params['convergence']:
                restart_params['convergence'][key] *= factor
        
        return restart_params
    
    def backup_calculation_files(self, calculation_dir: str, 
                               calculation_name: str) -> str:
        """Backup existing calculation files before restart"""
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_dir = os.path.join(calculation_dir, f"backup_{calculation_name}_{timestamp}")
        
        os.makedirs(backup_dir, exist_ok=True)
        
        # Files to backup
        files_to_backup = [
            f"{calculation_name}_dft.txt",
            f"{calculation_name}_opt.log",
            f"{calculation_name}.xyz",
            f"{calculation_name}_partial.xyz"
        ]
        
        backed_up_files = []
        for filename in files_to_backup:
            source_path = os.path.join(calculation_dir, filename)
            if os.path.exists(source_path):
                backup_path = os.path.join(backup_dir, filename)
                shutil.copy2(source_path, backup_path)
                backed_up_files.append(filename)
        
        # Save restart info
        restart_info = {
            'timestamp': datetime.now().isoformat(),
            'backed_up_files': backed_up_files,
            'backup_dir': backup_dir
        }
        
        with open(os.path.join(backup_dir, 'restart_info.json'), 'w') as f:
            json.dump(restart_info, f, indent=2)
        
        print(f"   📁 Backed up {len(backed_up_files)} files to {backup_dir}")
        return backup_dir
    
    def execute_restart(self, calculation_dir: str, calculation_name: str, 
                       atoms: Atoms, failure_type: str) -> bool:
        """
        Execute restart with improved parameters
        """
        
        print(f"\n🔄 Restarting calculation: {calculation_name}")
        print(f"   Failure type: {failure_type}")
        
        # Check restart history
        restart_count = len([r for r in self.restart_history 
                           if r['calculation_name'] == calculation_name])
        
        if restart_count >= 3:
            print(f"   ❌ Too many restart attempts ({restart_count}). Manual intervention needed.")
            return False
        
        # Backup existing files
        backup_dir = self.backup_calculation_files(calculation_dir, calculation_name)
        
        # Create improved parameters
        restart_params = self.create_restart_parameters(atoms, failure_type, restart_count + 1)
        
        print(f"   🎯 New parameters:")
        print(f"      Cutoff: {restart_params['mode'].ecut} eV")
        print(f"      Mixer β: {restart_params['mixer']['beta']:.4f}")
        print(f"      Max iter: {restart_params['maxiter']}")
        
        # Record restart attempt
        restart_record = {
            'timestamp': datetime.now().isoformat(),
            'calculation_name': calculation_name,
            'failure_type': failure_type,
            'restart_attempt': restart_count + 1,
            'backup_dir': backup_dir,
            'new_parameters': restart_params
        }
        
        self.restart_history.append(restart_record)
        
        # Save restart history
        with open(os.path.join(calculation_dir, 'restart_history.json'), 'w') as f:
            json.dump(self.restart_history, f, indent=2, default=str)
        
        print(f"   ✅ Restart prepared. Use improved parameters in your calculation script.")
        return True
    
    def monitor_and_restart(self, calculation_dir: str, calculation_name: str, 
                          atoms_file: str, check_interval: int = 300):
        """
        Monitor calculation and automatically restart if needed
        
        Args:
            calculation_dir: Directory containing calculation
            calculation_name: Base name of calculation
            atoms_file: Structure file for restart
            check_interval: Check interval in seconds (default 5 minutes)
        """
        
        print(f"\n🔍 Monitoring calculation: {calculation_name}")
        print(f"   Check interval: {check_interval/60:.1f} minutes")
        print("   Will auto-restart on failure")
        
        if not os.path.exists(atoms_file):
            print(f"❌ Structure file not found: {atoms_file}")
            return
        
        atoms = read(atoms_file)
        
        try:
            while True:
                status = self.detect_calculation_status(calculation_dir, calculation_name)
                
                print(f"   Status: {status['status']}", end="")
                if status['failure_type']:
                    print(f" | Failure: {status['failure_type']}")
                else:
                    print()
                
                # Check if restart needed
                if status['needs_restart']:
                    success = self.execute_restart(
                        calculation_dir, calculation_name, atoms, status['failure_type']
                    )
                    
                    if not success:
                        print("❌ Restart failed. Stopping monitoring.")
                        break
                
                # Check if completed
                if status['status'] in ['completed', 'converged']:
                    print("🎉 Calculation completed successfully!")
                    break
                
                time.sleep(check_interval)
                
        except KeyboardInterrupt:
            print("\n⏹️  Monitoring stopped by user")

def main():
    """Example usage"""
    import argparse
    
    parser = argparse.ArgumentParser(description='DFT Smart Restart System')
    parser.add_argument('--check', help='Check status of calculation')
    parser.add_argument('--monitor', help='Monitor calculation name')
    parser.add_argument('--atoms', help='Structure file for restart')
    parser.add_argument('--dir', default='.', help='Calculation directory')
    
    args = parser.parse_args()
    
    restart_system = DFTSmartRestart()
    
    if args.check:
        status = restart_system.detect_calculation_status(args.dir, args.check)
        print(json.dumps(status, indent=2, default=str))
    
    elif args.monitor and args.atoms:
        restart_system.monitor_and_restart(args.dir, args.monitor, args.atoms)
    
    else:
        print("Usage examples:")
        print("  python dft_smart_restart.py --check surface_opt")
        print("  python dft_smart_restart.py --monitor surface_opt --atoms surface.xyz")

if __name__ == "__main__":
    main()