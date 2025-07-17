#!/usr/bin/env python3
"""
MLD Data Analyzer - Terminal-based analysis tool
Provides detailed statistics and analysis without requiring GUI
"""

import os
import glob
import numpy as np
from ase.io import read
from collections import Counter, defaultdict
import json
from datetime import datetime

class MLDAnalyzer:
    """Comprehensive MLD data analysis"""
    
    def __init__(self):
        self.cycle_data = {}
        self.load_mld_data()
    
    def load_mld_data(self):
        """Load all MLD simulation data"""
        print("📊 Loading MLD data...")
        
        # Find all MLD files
        cycle_files = sorted(glob.glob('mld_cycle*_after_*.xyz'))
        
        for cycle_file in cycle_files:
            try:
                atoms = read(cycle_file)
                filename = os.path.basename(cycle_file)
                
                # Parse filename: mld_cycle1_after_tma.xyz
                parts = filename.replace('.xyz', '').split('_')
                cycle_num = int(parts[1].replace('cycle', ''))
                phase = parts[3]  # 'tma' or 'h2o'
                
                if cycle_num not in self.cycle_data:
                    self.cycle_data[cycle_num] = {}
                
                self.cycle_data[cycle_num][phase] = {
                    'file': cycle_file,
                    'atoms': atoms,
                    'composition': Counter(atoms.get_chemical_symbols()),
                    'total_atoms': len(atoms),
                    'positions': atoms.get_positions()
                }
                
            except Exception as e:
                print(f"Warning: Could not read {cycle_file}: {e}")
        
        print(f"✅ Loaded data for {len(self.cycle_data)} cycles")
    
    def analyze_growth_kinetics(self):
        """Analyze MLD growth kinetics"""
        print("\n" + "="*60)
        print("🔬 MLD GROWTH KINETICS ANALYSIS")
        print("="*60)
        
        if not self.cycle_data:
            print("❌ No cycle data available")
            return
        
        cycles = sorted(self.cycle_data.keys())
        
        # Growth per cycle analysis
        print("\n📈 Growth per Cycle:")
        print("-" * 50)
        print("Cycle │ After TMA │ After H2O │ Net Growth │ Efficiency")
        print("──────┼───────────┼───────────┼────────────┼───────────")
        
        total_growth = 0
        previous_atoms = None
        
        for cycle in cycles:
            tma_atoms = self.cycle_data[cycle].get('tma', {}).get('total_atoms', 0)
            h2o_atoms = self.cycle_data[cycle].get('h2o', {}).get('total_atoms', 0)
            
            if previous_atoms is not None:
                net_growth = h2o_atoms - previous_atoms
                total_growth += net_growth
                
                # Calculate efficiency (atoms added vs theoretical maximum)
                # Theoretical: each TMA adds ~8-10 atoms (Al + 2×CH3 groups)
                theoretical_max = 10  # Conservative estimate
                efficiency = (net_growth / theoretical_max) * 100 if theoretical_max > 0 else 0
                
                print(f"{cycle:5} │ {tma_atoms:9} │ {h2o_atoms:9} │ {net_growth:10} │ {efficiency:8.1f}%")
            else:
                print(f"{cycle:5} │ {tma_atoms:9} │ {h2o_atoms:9} │      -     │     -")
            
            previous_atoms = h2o_atoms
        
        print("-" * 50)
        print(f"Total growth: {total_growth} atoms")
        
        # Growth rate analysis
        if len(cycles) > 1:
            avg_growth = total_growth / (len(cycles) - 1)
            print(f"Average growth per cycle: {avg_growth:.1f} atoms")
    
    def analyze_composition_evolution(self):
        """Track how composition changes over cycles"""
        print("\n" + "="*60)
        print("🧪 COMPOSITION EVOLUTION ANALYSIS")
        print("="*60)
        
        cycles = sorted(self.cycle_data.keys())
        elements = set()
        
        # Collect all elements
        for cycle in cycles:
            for phase in ['tma', 'h2o']:
                if phase in self.cycle_data[cycle]:
                    comp = self.cycle_data[cycle][phase]['composition']
                    elements.update(comp.keys())
        
        elements = sorted(elements)
        
        print(f"\n📊 Element Evolution (After H2O phase):")
        print("-" * (15 + len(elements) * 8))
        
        # Header
        header = "Cycle │"
        for el in elements:
            header += f" {el:>6} │"
        print(header)
        print("─" * (len(header) - len("│")))
        
        # Data rows
        for cycle in cycles:
            if 'h2o' in self.cycle_data[cycle]:
                comp = self.cycle_data[cycle]['h2o']['composition']
                row = f"{cycle:5} │"
                for el in elements:
                    count = comp.get(el, 0)
                    row += f" {count:6} │"
                print(row)
        
        # Calculate growth rates for each element
        print(f"\n📈 Element Growth Rates:")
        print("-" * 40)
        
        if len(cycles) >= 2:
            first_cycle = cycles[0]
            last_cycle = cycles[-1]
            
            if 'h2o' in self.cycle_data[first_cycle] and 'h2o' in self.cycle_data[last_cycle]:
                first_comp = self.cycle_data[first_cycle]['h2o']['composition']
                last_comp = self.cycle_data[last_cycle]['h2o']['composition']
                
                for el in elements:
                    initial = first_comp.get(el, 0)
                    final = last_comp.get(el, 0)
                    growth = final - initial
                    
                    if initial > 0:
                        growth_percent = (growth / initial) * 100
                        print(f"{el:2}: {initial:4} → {final:4} (+{growth:3}) [{growth_percent:+6.1f}%]")
                    else:
                        print(f"{el:2}: {initial:4} → {final:4} (+{growth:3}) [new element]")
    
    def analyze_surface_coverage(self):
        """Analyze surface coverage and uniformity"""
        print("\n" + "="*60)
        print("🎯 SURFACE COVERAGE ANALYSIS")
        print("="*60)
        
        # Load initial surface for comparison
        initial_files = ['uv_ozone_si_surface.xyz', 'surface_si100_hydroxylated.xyz', 'simple_surface.xyz']
        initial_surface = None
        
        for init_file in initial_files:
            if os.path.exists(init_file):
                initial_surface = read(init_file)
                break
        
        if not initial_surface:
            print("❌ No initial surface found for comparison")
            return
        
        initial_atoms = len(initial_surface)
        initial_oh = initial_surface.get_chemical_symbols().count('H')  # OH groups
        
        print(f"\n📐 Initial Surface:")
        print(f"  Total atoms: {initial_atoms}")
        print(f"  OH groups: {initial_oh}")
        print(f"  Surface area: {self.calculate_surface_area(initial_surface):.1f} Ų")
        
        # Analyze coverage progression
        print(f"\n📊 Coverage Evolution:")
        print("-" * 60)
        print("Cycle │ Total Atoms │ Added │ OH Density │ Coverage")
        print("──────┼─────────────┼───────┼────────────┼─────────")
        
        cycles = sorted(self.cycle_data.keys())
        
        for cycle in cycles:
            if 'h2o' in self.cycle_data[cycle]:
                atoms = self.cycle_data[cycle]['h2o']['atoms']
                total_atoms = len(atoms)
                added_atoms = total_atoms - initial_atoms
                
                # Estimate OH density
                oh_count = atoms.get_chemical_symbols().count('H')
                surface_area = self.calculate_surface_area(atoms)
                oh_density = oh_count / surface_area if surface_area > 0 else 0
                
                # Estimate coverage (rough approximation)
                coverage = min(100, (added_atoms / initial_atoms) * 100)
                
                print(f"{cycle:5} │ {total_atoms:11} │ {added_atoms:5} │ {oh_density:8.2f}   │ {coverage:6.1f}%")
    
    def calculate_surface_area(self, atoms):
        """Rough surface area calculation"""
        cell = atoms.get_cell()
        if cell is not None:
            # For surfaces, use XY area
            area = np.linalg.norm(np.cross(cell[0], cell[1]))
            return area * 1e-2  # Convert Ų to nm²
        return 1.0  # Fallback
    
    def analyze_reaction_energetics(self):
        """Analyze energetics if available"""
        print("\n" + "="*60)
        print("⚡ REACTION ENERGETICS ANALYSIS")
        print("="*60)
        
        # Look for DFT output files
        dft_files = glob.glob('dft_*.xyz')
        gpaw_files = glob.glob('*_gpaw.txt')
        
        if not dft_files and not gpaw_files:
            print("❌ No DFT calculation results found")
            print("Run step5_true_dft_mld.py to get energetic data")
            return
        
        print(f"✅ Found DFT calculation files:")
        for f in dft_files[:5]:  # Show first 5
            print(f"  📄 {f}")
        
        if len(dft_files) > 5:
            print(f"  ... and {len(dft_files) - 5} more")
        
        # Try to extract energies from GPAW output
        energies = []
        for gpaw_file in gpaw_files:
            try:
                with open(gpaw_file, 'r') as f:
                    content = f.read()
                    # Look for final energy
                    if 'converged after' in content:
                        lines = content.split('\n')
                        for line in lines:
                            if 'Extrapolated:' in line and 'eV' in line:
                                energy_str = line.split()
                                for i, word in enumerate(energy_str):
                                    if word == 'eV' and i > 0:
                                        energy = float(energy_str[i-1])
                                        energies.append(energy)
                                        break
            except Exception as e:
                continue
        
        if energies:
            print(f"\n📊 Energy Analysis:")
            print(f"  Calculations found: {len(energies)}")
            print(f"  Energy range: {min(energies):.3f} to {max(energies):.3f} eV")
            print(f"  Average energy: {np.mean(energies):.3f} eV")
        else:
            print("⚠️  Could not extract energy data from GPAW files")
    
    def generate_summary_report(self):
        """Generate comprehensive summary report"""
        print("\n" + "="*80)
        print("📋 MLD SIMULATION SUMMARY REPORT")
        print("="*80)
        print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        
        # Basic statistics
        cycles = sorted(self.cycle_data.keys())
        if cycles:
            first_cycle = cycles[0]
            last_cycle = cycles[-1]
            
            if 'h2o' in self.cycle_data[last_cycle]:
                final_atoms = self.cycle_data[last_cycle]['h2o']['total_atoms']
                final_comp = self.cycle_data[last_cycle]['h2o']['composition']
                
                print(f"\n🎯 Key Results:")
                print(f"  Cycles completed: {len(cycles)}")
                print(f"  Final structure: {final_atoms} atoms")
                print(f"  Elements present: {', '.join(sorted(final_comp.keys()))}")
                
                # Growth efficiency
                if len(cycles) > 1 and 'h2o' in self.cycle_data[first_cycle]:
                    initial_atoms = self.cycle_data[first_cycle]['h2o']['total_atoms']
                    total_growth = final_atoms - initial_atoms
                    avg_growth = total_growth / (len(cycles) - 1)
                    
                    print(f"  Total growth: {total_growth} atoms")
                    print(f"  Average per cycle: {avg_growth:.1f} atoms")
        
        # File summary
        all_files = glob.glob('mld_*.xyz') + glob.glob('dft_*.xyz')
        print(f"\n📁 Output Files: {len(all_files)} files generated")
        
        # Recommendations
        print(f"\n💡 Recommendations:")
        
        if len(cycles) < 5:
            print("  • Consider running more cycles for better statistics")
        
        dft_files = glob.glob('dft_*.xyz')
        if not dft_files:
            print("  • Run DFT calculations for accurate energetics")
        
        print("  • Use terminal_visualizer.py for interactive analysis")
        print("  • Export data with save_analysis_data() for external tools")
    
    def save_analysis_data(self, filename='mld_analysis_data.json'):
        """Save analysis data to JSON for external tools"""
        print(f"\n💾 Saving analysis data to {filename}...")
        
        # Prepare data for JSON export
        export_data = {
            'metadata': {
                'generated': datetime.now().isoformat(),
                'total_cycles': len(self.cycle_data),
                'analysis_version': '1.0'
            },
            'cycles': {}
        }
        
        for cycle, data in self.cycle_data.items():
            export_data['cycles'][cycle] = {}
            
            for phase, phase_data in data.items():
                export_data['cycles'][cycle][phase] = {
                    'total_atoms': phase_data['total_atoms'],
                    'composition': dict(phase_data['composition']),
                    'file': phase_data['file']
                }
        
        with open(filename, 'w') as f:
            json.dump(export_data, f, indent=2)
        
        print(f"✅ Analysis data saved to {filename}")
        print(f"📊 Data can be imported into Python, R, or other analysis tools")

def main():
    """Main analysis function"""
    analyzer = MLDAnalyzer()
    
    if not analyzer.cycle_data:
        print("❌ No MLD simulation data found!")
        print("Please run step5_realistic_mld.py or step5_true_dft_mld.py first")
        return 1
    
    print(f"🔍 MLD Data Analyzer")
    print(f"Found data for {len(analyzer.cycle_data)} cycles")
    
    # Run all analyses
    analyzer.analyze_growth_kinetics()
    analyzer.analyze_composition_evolution()
    analyzer.analyze_surface_coverage()
    analyzer.analyze_reaction_energetics()
    analyzer.generate_summary_report()
    
    # Save data
    analyzer.save_analysis_data()
    
    print(f"\n✅ Analysis complete!")
    print(f"Use terminal_visualizer.py for interactive exploration")
    
    return 0

if __name__ == "__main__":
    exit(main())