#!/usr/bin/env python3
"""
DFT Convergence Detection and Analysis System
Monitors DFT calculations in real-time and provides intelligent analysis

Features:
- Real-time convergence monitoring
- Automatic problem detection
- Performance analysis
- Convergence prediction
- Smart alerts for intervention
"""

import os
import time
import re
import numpy as np
from datetime import datetime, timedelta
from typing import Dict, List, Tuple, Optional
import json

class DFTConvergenceMonitor:
    """Real-time DFT convergence monitoring and analysis"""
    
    def __init__(self):
        """Initialize convergence monitoring"""
        
        # Convergence patterns for different codes
        self.PATTERNS = {
            'gpaw_scf': re.compile(r'iter:\s+(\d+)\s+[\d:]+\s+([-\d.]+)\s+([-\d.c]+)\s+([-\d.c]+)'),
            'gpaw_opt': re.compile(r'LBFGS:\s+(\d+)\s+[\d:]+\s+([-\d.]+)\s+([\d.]+)'),
            'converged': re.compile(r'converged|Converged|CONVERGED')
        }
        
        # Convergence thresholds
        self.THRESHOLDS = {
            'energy_oscillation': 0.01,    # eV
            'slow_progress_steps': 10,      # iterations without improvement
            'density_warning': -2.0,       # log scale
            'eigenstates_warning': -2.0,   # log scale
            'max_reasonable_iter': 200,    # iterations before concern
        }
        
        # Performance tracking
        self.performance_history = []
        
        print("📊 DFT Convergence Monitor")
        print("   Real-time analysis and problem detection")
    
    def parse_dft_file(self, filepath: str) -> Dict:
        """Parse DFT output file for convergence data"""
        
        if not os.path.exists(filepath):
            return {'status': 'file_not_found', 'data': []}
        
        # Get file info
        file_size = os.path.getsize(filepath)
        modified_time = datetime.fromtimestamp(os.path.getmtime(filepath))
        
        convergence_data = {
            'file_info': {
                'size_kb': file_size / 1024,
                'modified': modified_time,
                'age_minutes': (datetime.now() - modified_time).total_seconds() / 60
            },
            'scf_iterations': [],
            'optimization_steps': [],
            'is_converged': False,
            'current_status': 'unknown'
        }
        
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            # Parse SCF iterations
            for line in lines:
                scf_match = self.PATTERNS['gpaw_scf'].search(line)
                if scf_match:
                    iter_num = int(scf_match.group(1))
                    energy = float(scf_match.group(2))
                    eigst = scf_match.group(3).replace('c', '')  # Remove 'c' flag
                    dens = scf_match.group(4).replace('c', '')
                    
                    convergence_data['scf_iterations'].append({
                        'iteration': iter_num,
                        'energy': energy,
                        'eigenstates_conv': float(eigst),
                        'density_conv': float(dens),
                        'eigst_converged': 'c' in scf_match.group(3),
                        'dens_converged': 'c' in scf_match.group(4)
                    })
                
                # Parse optimization steps
                opt_match = self.PATTERNS['gpaw_opt'].search(line)
                if opt_match:
                    step = int(opt_match.group(1))
                    energy = float(opt_match.group(2))
                    fmax = float(opt_match.group(3))
                    
                    convergence_data['optimization_steps'].append({
                        'step': step,
                        'energy': energy,
                        'max_force': fmax
                    })
                
                # Check for convergence
                if self.PATTERNS['converged'].search(line):
                    convergence_data['is_converged'] = True
            
            # Determine current status
            if convergence_data['is_converged']:
                convergence_data['current_status'] = 'converged'
            elif convergence_data['scf_iterations']:
                convergence_data['current_status'] = 'scf_running'
            elif convergence_data['optimization_steps']:
                convergence_data['current_status'] = 'optimization_running'
            else:
                convergence_data['current_status'] = 'starting'
                
        except Exception as e:
            convergence_data['error'] = str(e)
            convergence_data['current_status'] = 'error'
        
        return convergence_data
    
    def analyze_convergence_health(self, convergence_data: Dict) -> Dict:
        """Analyze convergence health and detect problems"""
        
        analysis = {
            'overall_health': 'good',
            'warnings': [],
            'problems': [],
            'recommendations': [],
            'performance': {},
            'predictions': {}
        }
        
        scf_data = convergence_data.get('scf_iterations', [])
        opt_data = convergence_data.get('optimization_steps', [])
        
        if not scf_data and not opt_data:
            analysis['overall_health'] = 'no_data'
            return analysis
        
        # Analyze SCF convergence
        if scf_data:
            self._analyze_scf_health(scf_data, analysis)
        
        # Analyze optimization convergence
        if opt_data:
            self._analyze_optimization_health(opt_data, analysis)
        
        # Performance analysis
        self._analyze_performance(convergence_data, analysis)
        
        # Overall health assessment
        if analysis['problems']:
            analysis['overall_health'] = 'poor'
        elif analysis['warnings']:
            analysis['overall_health'] = 'concerning'
        else:
            analysis['overall_health'] = 'good'
        
        return analysis
    
    def _analyze_scf_health(self, scf_data: List[Dict], analysis: Dict):
        """Analyze SCF convergence health"""
        
        if len(scf_data) < 3:
            return
        
        latest = scf_data[-1]
        recent = scf_data[-5:] if len(scf_data) >= 5 else scf_data
        
        # Check energy oscillations
        energies = [step['energy'] for step in recent]
        if len(energies) >= 3:
            energy_range = max(energies) - min(energies)
            mean_energy = abs(np.mean(energies))
            if energy_range > self.THRESHOLDS['energy_oscillation'] * mean_energy:
                analysis['problems'].append("Energy oscillating - reduce mixer beta")
        
        # Check convergence progress
        if latest['density_conv'] > self.THRESHOLDS['density_warning']:
            analysis['warnings'].append(f"Slow density convergence: {latest['density_conv']:.2f}")
        
        if latest['eigenstates_conv'] > self.THRESHOLDS['eigenstates_warning']:
            analysis['warnings'].append(f"Slow eigenstate convergence: {latest['eigenstates_conv']:.2f}")
        
        # Check iteration count
        if latest['iteration'] > self.THRESHOLDS['max_reasonable_iter']:
            analysis['problems'].append(f"High iteration count: {latest['iteration']}")
            analysis['recommendations'].append("Consider looser convergence criteria")
        
        # Predict convergence
        if len(scf_data) >= 10:
            self._predict_scf_convergence(scf_data, analysis)
    
    def _analyze_optimization_health(self, opt_data: List[Dict], analysis: Dict):
        """Analyze geometry optimization health"""
        
        if len(opt_data) < 2:
            return
        
        latest = opt_data[-1]
        
        # Check force convergence progress
        if len(opt_data) >= 5:
            recent_forces = [step['max_force'] for step in opt_data[-5:]]
            force_improvement = recent_forces[0] - recent_forces[-1]
            
            if force_improvement <= 0:
                analysis['warnings'].append("Forces not decreasing")
                analysis['recommendations'].append("Check for convergence issues")
        
        # Check energy conservation
        if len(opt_data) >= 3:
            recent_energies = [step['energy'] for step in opt_data[-3:]]
            energy_increase = max(recent_energies) - min(recent_energies)
            
            if energy_increase > 0.1:  # 0.1 eV increase
                analysis['warnings'].append("Energy increasing during optimization")
    
    def _analyze_performance(self, convergence_data: Dict, analysis: Dict):
        """Analyze calculation performance"""
        
        scf_data = convergence_data.get('scf_iterations', [])
        file_info = convergence_data.get('file_info', {})
        
        if len(scf_data) >= 2:
            # Estimate iteration time (rough approximation)
            total_iterations = len(scf_data)
            file_age_minutes = file_info.get('age_minutes', 1)
            
            avg_time_per_iter = file_age_minutes / total_iterations
            analysis['performance'] = {
                'avg_time_per_iteration': avg_time_per_iter,
                'total_iterations': total_iterations,
                'runtime_minutes': file_age_minutes
            }
            
            # Performance warnings
            if avg_time_per_iter > 2.0:  # More than 2 minutes per iteration
                analysis['warnings'].append(f"Slow iterations: {avg_time_per_iter:.1f} min/iter")
                analysis['recommendations'].append("Consider reducing cutoff energy")
    
    def _predict_scf_convergence(self, scf_data: List[Dict], analysis: Dict):
        """Predict SCF convergence time"""
        
        try:
            # Look at density convergence trend
            recent_dens = [step['density_conv'] for step in scf_data[-10:]]
            
            if len(recent_dens) >= 5:
                # Simple linear fit to predict convergence
                x = np.arange(len(recent_dens))
                coeffs = np.polyfit(x, recent_dens, 1)
                slope = coeffs[0]
                
                if slope < -0.1:  # Converging
                    # Predict when it reaches -4.0 (typical target)
                    current_conv = recent_dens[-1]
                    if current_conv > -4.0:
                        steps_to_convergence = int((-4.0 - current_conv) / slope)
                        analysis['predictions']['scf_convergence'] = {
                            'estimated_steps': steps_to_convergence,
                            'confidence': 'medium' if abs(slope) > 0.2 else 'low'
                        }
                else:
                    analysis['warnings'].append("SCF convergence has stalled")
                    
        except Exception:
            pass  # Prediction failed, no big deal
    
    def generate_status_report(self, filepath: str) -> Dict:
        """Generate comprehensive status report"""
        
        convergence_data = self.parse_dft_file(filepath)
        health_analysis = self.analyze_convergence_health(convergence_data)
        
        report = {
            'timestamp': datetime.now().isoformat(),
            'file_path': filepath,
            'convergence_data': convergence_data,
            'health_analysis': health_analysis,
            'summary': self._create_summary(convergence_data, health_analysis)
        }
        
        return report
    
    def _create_summary(self, convergence_data: Dict, health_analysis: Dict) -> Dict:
        """Create human-readable summary"""
        
        summary = {
            'status': convergence_data.get('current_status', 'unknown'),
            'health': health_analysis.get('overall_health', 'unknown'),
            'progress': {},
            'time_info': {},
            'key_metrics': {}
        }
        
        # Progress summary
        scf_data = convergence_data.get('scf_iterations', [])
        opt_data = convergence_data.get('optimization_steps', [])
        
        if scf_data:
            latest_scf = scf_data[-1]
            summary['progress']['scf'] = {
                'iterations': latest_scf['iteration'],
                'energy': latest_scf['energy'],
                'density_conv': latest_scf['density_conv'],
                'eigst_converged': latest_scf['eigst_converged'],
                'dens_converged': latest_scf['dens_converged']
            }
        
        if opt_data:
            latest_opt = opt_data[-1]
            summary['progress']['optimization'] = {
                'steps': latest_opt['step'],
                'energy': latest_opt['energy'],
                'max_force': latest_opt['max_force']
            }
        
        # Time info
        file_info = convergence_data.get('file_info', {})
        summary['time_info'] = {
            'runtime_minutes': file_info.get('age_minutes', 0),
            'last_update': file_info.get('modified')
        }
        
        # Performance metrics
        if 'performance' in health_analysis:
            summary['key_metrics'] = health_analysis['performance']
        
        return summary
    
    def monitor_calculation(self, filepath: str, check_interval: int = 30, 
                          save_reports: bool = True) -> None:
        """Monitor calculation with periodic reporting"""
        
        print(f"\n📊 Monitoring: {filepath}")
        print(f"   Check interval: {check_interval} seconds")
        print("   Press Ctrl+C to stop monitoring\n")
        
        try:
            while True:
                report = self.generate_status_report(filepath)
                
                # Display summary
                self._display_monitoring_summary(report)
                
                # Save report if requested
                if save_reports:
                    report_file = f"{os.path.splitext(filepath)[0]}_convergence_report.json"
                    with open(report_file, 'w') as f:
                        json.dump(report, f, indent=2, default=str)
                
                # Check if converged
                if report['convergence_data']['is_converged']:
                    print("\n🎉 Calculation converged!")
                    break
                
                time.sleep(check_interval)
                
        except KeyboardInterrupt:
            print("\n⏹️  Monitoring stopped by user")
    
    def _display_monitoring_summary(self, report: Dict):
        """Display concise monitoring summary"""
        
        summary = report['summary']
        health = report['health_analysis']
        
        # Clear screen and show status
        print(f"\r📊 Status: {summary['status']} | Health: {summary['health']} | Runtime: {summary['time_info']['runtime_minutes']:.1f} min", end="")
        
        # Show warnings if any
        if health['warnings']:
            print(f"\n⚠️  Warnings: {', '.join(health['warnings'])}")
        
        # Show problems if any
        if health['problems']:
            print(f"\n🚨 Problems: {', '.join(health['problems'])}")

def main():
    """Example usage and testing"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Monitor DFT convergence')
    parser.add_argument('dft_file', help='DFT output file to monitor')
    parser.add_argument('--interval', type=int, default=30, help='Check interval (seconds)')
    parser.add_argument('--report', action='store_true', help='Generate detailed report')
    
    args = parser.parse_args()
    
    monitor = DFTConvergenceMonitor()
    
    if args.report:
        # Generate single report
        report = monitor.generate_status_report(args.dft_file)
        print(json.dumps(report, indent=2, default=str))
    else:
        # Start monitoring
        monitor.monitor_calculation(args.dft_file, args.interval)

if __name__ == "__main__":
    main()