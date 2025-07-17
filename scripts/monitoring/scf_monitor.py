#!/usr/bin/env python3
"""
Simple SCF Iteration Monitor
Shows real-time SCF convergence progress with timing information
"""

import os
import time
import re
from datetime import datetime
import sys

class SCFMonitor:
    """Monitor SCF iterations in real-time"""
    
    def __init__(self):
        self.scf_pattern = re.compile(r'iter:\s+(\d+)\s+(\d+:\d+:\d+)\s+([-\d.]+)')
        self.converged_pattern = re.compile(r'converged')
        
    def monitor_scf(self, filepath, refresh_rate=1.0):
        """Monitor SCF iterations with live updates"""
        if not os.path.exists(filepath):
            print(f"❌ File not found: {filepath}")
            return
            
        print(f"🔍 Monitoring SCF progress in: {filepath}")
        print("Press Ctrl+C to stop monitoring\n")
        
        last_size = 0
        scf_history = []
        start_time = None
        
        while True:
            try:
                # Check if file has been updated
                current_size = os.path.getsize(filepath)
                if current_size > last_size:
                    
                    # Read new content
                    with open(filepath, 'r') as f:
                        lines = f.readlines()
                    
                    # Parse SCF iterations
                    for line in lines:
                        scf_match = self.scf_pattern.search(line)
                        if scf_match:
                            iter_num = int(scf_match.group(1))
                            time_str = scf_match.group(2)
                            energy = float(scf_match.group(3))
                            
                            # Track timing
                            if start_time is None:
                                start_time = datetime.now()
                            
                            # Update history
                            scf_data = {
                                'iter': iter_num,
                                'time': time_str,
                                'energy': energy
                            }
                            
                            # Only add new iterations
                            if not scf_history or iter_num > scf_history[-1]['iter']:
                                scf_history.append(scf_data)
                                
                                # Calculate progress metrics
                                if len(scf_history) > 1:
                                    prev_energy = scf_history[-2]['energy']
                                    energy_change = energy - prev_energy
                                    
                                    # Estimate time per iteration
                                    elapsed = (datetime.now() - start_time).total_seconds()
                                    time_per_iter = elapsed / iter_num if iter_num > 0 else 0
                                    
                                    # Progress bar
                                    progress = min(iter_num / 300 * 100, 100)
                                    bar_length = 30
                                    filled = int(bar_length * progress / 100)
                                    bar = "█" * filled + "░" * (bar_length - filled)
                                    
                                    print(f"\r🔄 SCF [{bar}] {progress:5.1f}% | "
                                          f"Iter: {iter_num:3d}/300 | "
                                          f"E: {energy:12.6f} eV | "
                                          f"ΔE: {energy_change:8.2e} | "
                                          f"Time: {time_str} | "
                                          f"Rate: {time_per_iter:.1f}s/iter", end="")
                                else:
                                    print(f"\r🔄 SCF Starting... | Iter: {iter_num:3d}/300 | "
                                          f"E: {energy:12.6f} eV | Time: {time_str}", end="")
                        
                        # Check for convergence
                        if self.converged_pattern.search(line):
                            print(f"\n✅ SCF Converged in {len(scf_history)} iterations!")
                            if scf_history:
                                final_energy = scf_history[-1]['energy']
                                elapsed = (datetime.now() - start_time).total_seconds()
                                print(f"   Final Energy: {final_energy:.6f} eV")
                                print(f"   Total Time: {elapsed:.1f} seconds")
                                print(f"   Average: {elapsed/len(scf_history):.1f} s/iteration")
                            return
                    
                    last_size = current_size
                
                time.sleep(refresh_rate)
                
            except KeyboardInterrupt:
                print(f"\n\n🛑 Monitoring stopped by user")
                if scf_history:
                    print(f"   Last iteration: {scf_history[-1]['iter']}")
                    print(f"   Last energy: {scf_history[-1]['energy']:.6f} eV")
                break
            except Exception as e:
                print(f"\n❌ Error: {e}")
                break

def main():
    """Main execution"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Monitor SCF convergence progress')
    parser.add_argument('file', help='Path to DFT output file')
    parser.add_argument('--refresh', type=float, default=1.0, help='Refresh rate in seconds')
    
    args = parser.parse_args()
    
    monitor = SCFMonitor()
    monitor.monitor_scf(args.file, args.refresh)

if __name__ == "__main__":
    main()