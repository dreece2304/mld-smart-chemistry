#!/usr/bin/env python3
"""
Real-time DFT Progress Monitor with Progress Bars
Monitors GPAW output files and displays progress bars for:
- SCF convergence
- Geometry optimization
- Overall simulation progress
"""

import os
import time
import re
from datetime import datetime
from typing import Dict, List, Optional, Tuple
import sys

try:
    from rich.console import Console
    from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeRemainingColumn
    from rich.table import Table
    from rich.live import Live
    from rich.layout import Layout
    from rich.panel import Panel
    RICH_AVAILABLE = True
except ImportError:
    RICH_AVAILABLE = False
    print("Installing rich for better progress display...")
    os.system("pip install rich")
    try:
        from rich.console import Console
        from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeRemainingColumn
        from rich.table import Table
        from rich.live import Live
        from rich.layout import Layout
        from rich.panel import Panel
        RICH_AVAILABLE = True
    except:
        print("Failed to install rich. Using basic progress display.")

class DFTProgressMonitor:
    """Monitor DFT calculation progress in real-time"""
    
    def __init__(self):
        self.console = Console() if RICH_AVAILABLE else None
        self.scf_pattern = re.compile(r'iter:\s+(\d+)\s+.+')
        self.energy_pattern = re.compile(r'Free energy:\s+([-\d.]+)')
        self.force_pattern = re.compile(r'Maximum force:\s+([\d.]+)')
        self.opt_step_pattern = re.compile(r'Step=\s*(\d+)')
        self.gpaw_scf_pattern = re.compile(r'(\d+)\s+([-\d.]+)\s+\d+\s+([-\d.]+)')
        self.lbfgs_pattern = re.compile(r'LBFGS:\s+(\d+)\s+[\d:]+\s+([-\d.]+)\s+([\d.]+)')
        
    def parse_dft_file(self, filepath: str) -> Dict:
        """Parse DFT output file for progress information"""
        if not os.path.exists(filepath):
            return {'exists': False}
        
        stats = {
            'exists': True,
            'size': os.path.getsize(filepath),
            'modified': datetime.fromtimestamp(os.path.getmtime(filepath)),
            'scf_iter': 0,
            'scf_energy': None,
            'opt_step': 0,
            'max_force': None,
            'converged': False,
            'scf_history': [],
            'opt_history': []
        }
        
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
                
            for line in lines:
                # Check for SCF iterations (GPAW style)
                gpaw_match = self.gpaw_scf_pattern.search(line)
                if gpaw_match:
                    stats['scf_iter'] = int(gpaw_match.group(1))
                    stats['scf_energy'] = float(gpaw_match.group(2))
                    stats['scf_history'].append({
                        'iter': stats['scf_iter'],
                        'energy': stats['scf_energy']
                    })
                
                # Check for optimization steps (LBFGS)
                lbfgs_match = self.lbfgs_pattern.search(line)
                if lbfgs_match:
                    stats['opt_step'] = int(lbfgs_match.group(1))
                    stats['opt_energy'] = float(lbfgs_match.group(2))
                    stats['max_force'] = float(lbfgs_match.group(3))
                    stats['opt_history'].append({
                        'step': stats['opt_step'],
                        'energy': stats['opt_energy'],
                        'force': stats['max_force']
                    })
                
                # Check for convergence
                if 'converged' in line.lower() and 'not' not in line.lower():
                    stats['converged'] = True
                    
        except Exception as e:
            stats['error'] = str(e)
            
        return stats
    
    def parse_opt_log(self, filepath: str) -> Dict:
        """Parse optimization log file"""
        if not os.path.exists(filepath):
            return {'exists': False}
        
        stats = {
            'exists': True,
            'steps': [],
            'current_step': 0,
            'current_fmax': None,
            'current_energy': None
        }
        
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    if 'Step' in line and 'Time' in line and 'Energy' in line:
                        continue  # Header line
                    
                    parts = line.strip().split()
                    if len(parts) >= 4 and parts[0].isdigit():
                        step_data = {
                            'step': int(parts[0]),
                            'time': parts[1],
                            'energy': float(parts[2]),
                            'fmax': float(parts[3])
                        }
                        stats['steps'].append(step_data)
                        stats['current_step'] = step_data['step']
                        stats['current_energy'] = step_data['energy']
                        stats['current_fmax'] = step_data['fmax']
                        
        except Exception as e:
            stats['error'] = str(e)
            
        return stats
    
    def create_progress_display(self, stats: Dict, opt_stats: Dict) -> str:
        """Create a formatted progress display"""
        if not RICH_AVAILABLE:
            return self.create_simple_display(stats, opt_stats)
        
        # Create layout
        layout = Layout()
        
        # File info
        if stats['exists']:
            file_info = Table(show_header=False)
            file_info.add_row("File Size", f"{stats['size'] / 1024:.1f} KB")
            file_info.add_row("Last Modified", stats['modified'].strftime("%H:%M:%S"))
            
            # SCF Progress
            scf_progress = Progress(
                TextColumn("[bold blue]SCF Iteration"),
                BarColumn(),
                TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                TextColumn("iter {task.completed}/{task.total}")
            )
            
            max_scf = 300  # Typical max iterations
            current_scf = stats.get('scf_iter', 0)
            scf_task = scf_progress.add_task("SCF", total=max_scf, completed=current_scf)
            
            # Optimization Progress
            if opt_stats['exists'] and opt_stats['steps']:
                opt_progress = Progress(
                    TextColumn("[bold green]Optimization"),
                    BarColumn(),
                    TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                    TextColumn("step {task.completed}/{task.total}")
                )
                
                max_opt = 100  # Typical max steps
                current_opt = opt_stats['current_step']
                opt_task = opt_progress.add_task("Opt", total=max_opt, completed=current_opt)
                
                # Force convergence progress
                force_progress = Progress(
                    TextColumn("[bold yellow]Force Conv"),
                    BarColumn(),
                    TextColumn("{task.description}")
                )
                
                target_force = 0.03  # Target fmax
                current_force = opt_stats['current_fmax'] or 999
                force_percentage = max(0, min(100, (1 - current_force/5.0) * 100))
                force_task = force_progress.add_task(
                    f"Force", 
                    total=100, 
                    completed=force_percentage,
                    description=f"fmax={current_force:.4f} eV/Å (target: {target_force})"
                )
            
            # Energy info
            energy_table = Table(title="Current Values")
            energy_table.add_column("Property", style="cyan")
            energy_table.add_column("Value", style="green")
            
            if stats.get('scf_energy'):
                energy_table.add_row("SCF Energy", f"{stats['scf_energy']:.6f} eV")
            if opt_stats['exists'] and opt_stats['current_energy']:
                energy_table.add_row("Total Energy", f"{opt_stats['current_energy']:.6f} eV")
            if opt_stats['exists'] and opt_stats['current_fmax']:
                energy_table.add_row("Max Force", f"{opt_stats['current_fmax']:.6f} eV/Å")
            
            # Combine into panels
            panels = [
                Panel(file_info, title="File Info"),
                Panel(scf_progress, title="SCF Convergence"),
            ]
            
            if opt_stats['exists'] and opt_stats['steps']:
                panels.extend([
                    Panel(opt_progress, title="Geometry Optimization"),
                    Panel(force_progress, title="Force Convergence"),
                ])
            
            panels.append(Panel(energy_table, title="Current State"))
            
            return panels
        else:
            return [Panel("Waiting for DFT output file...", title="Status")]
    
    def create_simple_display(self, stats: Dict, opt_stats: Dict) -> str:
        """Create simple text progress display"""
        output = ["\n" + "="*60]
        output.append(f"DFT Progress Monitor - {datetime.now().strftime('%H:%M:%S')}")
        output.append("="*60)
        
        if stats['exists']:
            output.append(f"File size: {stats['size']/1024:.1f} KB")
            output.append(f"Last modified: {stats['modified'].strftime('%H:%M:%S')}")
            output.append("")
            
            if stats.get('scf_iter'):
                scf_bar = self.create_text_bar(stats['scf_iter'], 300, 30)
                output.append(f"SCF: {scf_bar} {stats['scf_iter']}/300")
                
            if stats.get('scf_energy'):
                output.append(f"Energy: {stats['scf_energy']:.6f} eV")
                
            if opt_stats['exists'] and opt_stats['steps']:
                opt_bar = self.create_text_bar(opt_stats['current_step'], 100, 30)
                output.append(f"Opt: {opt_bar} {opt_stats['current_step']}/100")
                output.append(f"Force: {opt_stats['current_fmax']:.4f} eV/Å")
                
        else:
            output.append("Waiting for DFT output file...")
            
        return "\n".join(output)
    
    def create_text_bar(self, current: int, total: int, width: int = 30) -> str:
        """Create a simple text progress bar"""
        percentage = current / total
        filled = int(width * percentage)
        bar = "█" * filled + "░" * (width - filled)
        return f"[{bar}] {percentage*100:.0f}%"
    
    def monitor_calculation(self, dft_file: str, opt_log: str = None, refresh_rate: float = 1.0):
        """Monitor DFT calculation with live updates"""
        print(f"Monitoring: {dft_file}")
        if opt_log:
            print(f"Optimization log: {opt_log}")
        print("Press Ctrl+C to stop monitoring\n")
        
        if RICH_AVAILABLE:
            with Live(console=self.console, refresh_per_second=1/refresh_rate) as live:
                while True:
                    try:
                        stats = self.parse_dft_file(dft_file)
                        opt_stats = self.parse_opt_log(opt_log) if opt_log else {'exists': False}
                        
                        display = self.create_progress_display(stats, opt_stats)
                        
                        # Update live display
                        if isinstance(display, list):
                            from rich.columns import Columns
                            live.update(Columns(display))
                        else:
                            live.update(display)
                        
                        # Check if converged
                        if stats.get('converged'):
                            live.console.print("\n[bold green]✓ Calculation converged![/bold green]")
                            break
                            
                        time.sleep(refresh_rate)
                        
                    except KeyboardInterrupt:
                        live.console.print("\n[yellow]Monitoring stopped by user[/yellow]")
                        break
                    except Exception as e:
                        live.console.print(f"\n[red]Error: {e}[/red]")
                        break
        else:
            # Simple monitoring without rich
            while True:
                try:
                    os.system('clear' if os.name == 'posix' else 'cls')
                    stats = self.parse_dft_file(dft_file)
                    opt_stats = self.parse_opt_log(opt_log) if opt_log else {'exists': False}
                    print(self.create_simple_display(stats, opt_stats))
                    
                    if stats.get('converged'):
                        print("\n✓ Calculation converged!")
                        break
                        
                    time.sleep(refresh_rate)
                    
                except KeyboardInterrupt:
                    print("\nMonitoring stopped by user")
                    break

def main():
    """Main execution"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Monitor DFT calculation progress')
    parser.add_argument('dft_file', help='Path to DFT output file (e.g., tma_dft.txt)')
    parser.add_argument('--opt-log', help='Path to optimization log file (e.g., cycle1_tma_opt.log)')
    parser.add_argument('--refresh', type=float, default=1.0, help='Refresh rate in seconds')
    
    args = parser.parse_args()
    
    monitor = DFTProgressMonitor()
    monitor.monitor_calculation(args.dft_file, args.opt_log, args.refresh)

if __name__ == "__main__":
    main()