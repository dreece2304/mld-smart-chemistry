#!/usr/bin/env python3
"""
Real-time monitoring script for MLD calculations
Provides live updates on running calculations
"""

import os
import sys
import time
import glob
import argparse
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np

try:
    from tqdm.auto import tqdm
except ImportError:
    os.system("pip install tqdm")
    from tqdm.auto import tqdm

class CalculationMonitor:
    """Monitor running GPAW calculations"""
    
    def __init__(self, watch_dir='.', update_interval=5):
        self.watch_dir = watch_dir
        self.update_interval = update_interval
        self.start_time = time.time()
        
        # Storage for monitoring data
        self.calculation_data = {}
        self.file_sizes = {}
        
        print(f"🔍 Monitoring calculations in: {os.path.abspath(watch_dir)}")
        print(f"🔄 Update interval: {update_interval} seconds")
        print(f"⏰ Started: {datetime.now().strftime('%H:%M:%S')}")
        print("=" * 60)
    
    def find_active_calculations(self):
        """Find currently running calculations"""
        
        # Look for GPAW output files
        gpaw_files = glob.glob(os.path.join(self.watch_dir, '**/*_gpaw.out'), recursive=True)
        
        # Look for optimization log files
        opt_files = glob.glob(os.path.join(self.watch_dir, '**/*_opt.log'), recursive=True)
        
        # Look for trajectory files
        traj_files = glob.glob(os.path.join(self.watch_dir, '**/*.traj'), recursive=True)
        
        active_files = set(gpaw_files + opt_files + traj_files)
        
        # Filter for recently modified files (active calculations)
        active_calcs = []
        current_time = time.time()
        
        for file_path in active_files:
            if os.path.exists(file_path):
                mtime = os.path.getmtime(file_path)
                if current_time - mtime < 300:  # Modified in last 5 minutes
                    active_calcs.append(file_path)
        
        return active_calcs
    
    def parse_gpaw_output(self, file_path):
        """Parse GPAW output file for progress information"""
        
        if not os.path.exists(file_path):
            return None
        
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
        except:
            return None
        
        data = {
            'energies': [],
            'forces': [],
            'scf_iterations': [],
            'current_step': 0,
            'converged': False,
            'last_update': datetime.now()
        }
        
        # Parse SCF iterations and energies
        for line in lines:
            if 'iter:' in line and 'eV' in line:
                try:
                    parts = line.split()
                    if len(parts) >= 4:
                        energy = float(parts[3])
                        data['energies'].append(energy)
                except:
                    continue
            
            elif 'FORCES:' in line or 'max|F|' in line:
                try:
                    # Extract force information
                    force_val = float(line.split()[-1])
                    data['forces'].append(force_val)
                except:
                    continue
            
            elif 'Converged' in line:
                data['converged'] = True
        
        data['current_step'] = len(data['energies'])
        
        return data
    
    def parse_optimization_log(self, file_path):
        """Parse ASE optimization log file"""
        
        if not os.path.exists(file_path):
            return None
        
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
        except:
            return None
        
        data = {
            'step': 0,
            'energy': None,
            'fmax': None,
            'last_update': datetime.now()
        }
        
        # Parse optimization steps
        for line in reversed(lines[-20:]):  # Check last 20 lines
            if line.strip() and not line.startswith('#'):
                try:
                    parts = line.split()
                    if len(parts) >= 3:
                        data['step'] = int(parts[0])
                        data['energy'] = float(parts[1])
                        data['fmax'] = float(parts[2])
                        break
                except:
                    continue
        
        return data
    
    def get_file_size_info(self, file_path):
        """Get file size and growth rate"""
        
        if not os.path.exists(file_path):
            return None
        
        current_size = os.path.getsize(file_path)
        current_time = time.time()
        
        if file_path not in self.file_sizes:
            self.file_sizes[file_path] = [(current_time, current_size)]
            return {
                'size_mb': current_size / (1024 * 1024),
                'growth_rate': 0
            }
        
        # Keep last 10 measurements
        self.file_sizes[file_path].append((current_time, current_size))
        self.file_sizes[file_path] = self.file_sizes[file_path][-10:]
        
        # Calculate growth rate
        if len(self.file_sizes[file_path]) > 1:
            old_time, old_size = self.file_sizes[file_path][0]
            time_diff = current_time - old_time
            size_diff = current_size - old_size
            growth_rate = size_diff / time_diff if time_diff > 0 else 0
        else:
            growth_rate = 0
        
        return {
            'size_mb': current_size / (1024 * 1024),
            'growth_rate': growth_rate / (1024 * 1024)  # MB/s
        }
    
    def display_status(self, active_calcs):
        """Display current status of all calculations"""
        
        # Clear screen
        os.system('clear' if os.name == 'posix' else 'cls')
        
        print(f"🔍 MLD CALCULATION MONITOR")
        print(f"⏰ {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | "
              f"Runtime: {str(timedelta(seconds=int(time.time() - self.start_time)))}")
        print("=" * 80)
        
        if not active_calcs:
            print("😴 No active calculations detected")
            print("\nLooking for files modified in the last 5 minutes...")
            print("If calculations are running, they should appear here.")
            return
        
        for i, file_path in enumerate(active_calcs):
            rel_path = os.path.relpath(file_path, self.watch_dir)
            file_info = self.get_file_size_info(file_path)
            
            print(f"\n📊 Calculation {i+1}: {rel_path}")
            print("-" * 40)
            
            if file_info:
                print(f"💾 File size: {file_info['size_mb']:.2f} MB | "
                      f"Growth: {file_info['growth_rate']:.3f} MB/s")
            
            # Parse specific file types
            if file_path.endswith('_gpaw.out'):
                data = self.parse_gpaw_output(file_path)
                if data:
                    print(f"🔧 GPAW calculation:")
                    print(f"   SCF steps: {data['current_step']}")
                    if data['energies']:
                        print(f"   Last energy: {data['energies'][-1]:.6f} eV")
                    if data['forces']:
                        print(f"   Last max force: {data['forces'][-1]:.6f} eV/Å")
                    print(f"   Status: {'✅ Converged' if data['converged'] else '🔄 Running'}")
            
            elif file_path.endswith('_opt.log'):
                data = self.parse_optimization_log(file_path)
                if data:
                    print(f"📈 Optimization:")
                    print(f"   Step: {data['step']}")
                    if data['energy']:
                        print(f"   Energy: {data['energy']:.6f} eV")
                    if data['fmax']:
                        print(f"   Max force: {data['fmax']:.6f} eV/Å")
            
            # Show last modification time
            mtime = os.path.getmtime(file_path)
            last_mod = datetime.fromtimestamp(mtime)
            time_diff = datetime.now() - last_mod
            print(f"⏱️  Last updated: {time_diff.seconds} seconds ago")
        
        print("\n" + "=" * 80)
        print(f"Press Ctrl+C to stop monitoring | Next update in {self.update_interval}s")
    
    def plot_convergence(self, gpaw_files):
        """Plot convergence for GPAW calculations"""
        
        plt.figure(figsize=(12, 8))
        
        for i, file_path in enumerate(gpaw_files):
            data = self.parse_gpaw_output(file_path)
            if data and data['energies']:
                plt.subplot(2, 2, 1)
                plt.plot(data['energies'], label=f'Calc {i+1}')
                plt.xlabel('SCF Iteration')
                plt.ylabel('Energy (eV)')
                plt.title('Energy Convergence')
                plt.legend()
                
                if data['forces']:
                    plt.subplot(2, 2, 2)
                    plt.plot(data['forces'], label=f'Calc {i+1}')
                    plt.xlabel('Step')
                    plt.ylabel('Max Force (eV/Å)')
                    plt.title('Force Convergence')
                    plt.legend()
        
        plt.tight_layout()
        plt.savefig('convergence_plot.png', dpi=150)
        plt.close()
        
        print("📈 Convergence plot saved: convergence_plot.png")
    
    def start_monitoring(self, plot_interval=60):
        """Start the monitoring loop"""
        
        plot_counter = 0
        
        try:
            while True:
                active_calcs = self.find_active_calculations()
                self.display_status(active_calcs)
                
                # Generate plots periodically
                plot_counter += self.update_interval
                if plot_counter >= plot_interval:
                    gpaw_files = [f for f in active_calcs if f.endswith('_gpaw.out')]
                    if gpaw_files:
                        try:
                            self.plot_convergence(gpaw_files)
                        except Exception as e:
                            print(f"⚠️  Plot generation failed: {e}")
                    plot_counter = 0
                
                time.sleep(self.update_interval)
                
        except KeyboardInterrupt:
            print("\n\n⚠️  Monitoring stopped by user")
        except Exception as e:
            print(f"\n\n❌ Monitoring error: {e}")

def main():
    """Main function"""
    
    parser = argparse.ArgumentParser(description='Monitor MLD calculations in real-time')
    parser.add_argument('--dir', default='.', 
                       help='Directory to monitor (default: current)')
    parser.add_argument('--interval', type=int, default=5,
                       help='Update interval in seconds (default: 5)')
    parser.add_argument('--plot-interval', type=int, default=60,
                       help='Plot generation interval in seconds (default: 60)')
    parser.add_argument('--once', action='store_true',
                       help='Run once and exit (no continuous monitoring)')
    
    args = parser.parse_args()
    
    monitor = CalculationMonitor(args.dir, args.interval)
    
    if args.once:
        # Single check
        active_calcs = monitor.find_active_calculations()
        monitor.display_status(active_calcs)
    else:
        # Continuous monitoring
        monitor.start_monitoring(args.plot_interval)

if __name__ == "__main__":
    main()