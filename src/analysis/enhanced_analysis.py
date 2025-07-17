#!/usr/bin/env python3
"""
Enhanced analysis tools with progress tracking
Integrates bulk structure analysis and radiation modeling with progress bars
"""

import os
import sys
import time
import numpy as np
from datetime import datetime

try:
    from tqdm.auto import tqdm
except ImportError:
    os.system("pip install tqdm")
    from tqdm.auto import tqdm

class ProgressAnalysis:
    """Enhanced analysis with progress tracking"""
    
    def __init__(self, structure_file, name="Analysis"):
        self.structure_file = structure_file
        self.name = name
        self.start_time = time.time()
        
        # Check if file exists
        if not os.path.exists(structure_file):
            raise FileNotFoundError(f"Structure file not found: {structure_file}")
        
        print(f"🔬 Enhanced Analysis: {name}")
        print(f"📁 Structure: {structure_file}")
        print(f"🕒 Started: {datetime.now().strftime('%H:%M:%S')}")
        print("=" * 60)
    
    def run_bulk_analysis_with_progress(self):
        """Run bulk structure analysis with progress tracking"""
        
        try:
            from analysis.bulk_structure_analysis import BulkMLDAnalyzer
        except ImportError:
            print("❌ Bulk analysis module not available")
            return None
        
        print("\n📊 Bulk Structure Analysis")
        print("-" * 40)
        
        # Create progress bar for analysis steps
        analysis_steps = [
            "Loading structure",
            "Analyzing crystallinity", 
            "Calculating density",
            "Analyzing layer structure",
            "Computing mechanical properties",
            "Analyzing bonding",
            "Calculating XRD pattern",
            "Generating report"
        ]
        
        with tqdm(total=len(analysis_steps), desc="🔍 Bulk Analysis", 
                 unit="steps", colour='blue') as pbar:
            
            # Step 1: Load structure
            pbar.set_description("🔍 Loading structure")
            analyzer = BulkMLDAnalyzer(self.structure_file)
            pbar.update(1)
            time.sleep(0.5)  # Small delay for visual feedback
            
            # Step 2: Crystallinity
            pbar.set_description("🔍 Analyzing crystallinity")
            crystallinity = analyzer.analyze_crystallinity()
            pbar.update(1)
            time.sleep(0.5)
            
            # Step 3: Density
            pbar.set_description("🔍 Calculating density")
            density = analyzer.calculate_density()
            pbar.update(1)
            time.sleep(0.5)
            
            # Step 4: Layer structure
            pbar.set_description("🔍 Analyzing layers")
            layer_info = analyzer.analyze_layer_structure()
            pbar.update(1)
            time.sleep(0.5)
            
            # Step 5: Mechanical properties
            pbar.set_description("🔍 Computing mechanics")
            mechanical = analyzer.calculate_mechanical_properties()
            pbar.update(1)
            time.sleep(0.5)
            
            # Step 6: Bonding
            pbar.set_description("🔍 Analyzing bonding")
            bonding = analyzer.analyze_bonding()
            pbar.update(1)
            time.sleep(0.5)
            
            # Step 7: XRD
            pbar.set_description("🔍 Calculating XRD")
            try:
                xrd = analyzer.calculate_xrd_pattern()
            except:
                xrd = None  # XRD might fail for some structures
            pbar.update(1)
            time.sleep(0.5)
            
            # Step 8: Report
            pbar.set_description("🔍 Generating report")
            output_file = f"{self.name.lower()}_bulk_analysis.txt"
            analyzer.generate_report(output_file)
            pbar.update(1)
            
            pbar.set_description("✅ Bulk analysis complete")
        
        results = {
            'crystallinity': crystallinity,
            'density': density,
            'layer_info': layer_info,
            'mechanical': mechanical,
            'bonding': bonding,
            'xrd': xrd,
            'report_file': output_file
        }
        
        print(f"📄 Bulk analysis report: {output_file}")
        return results
    
    def run_radiation_analysis_with_progress(self, exposures=None):
        """Run radiation damage analysis with progress tracking"""
        
        try:
            from calculations.radiation.radiation_damage import RadiationDamageSimulator
            from calculations.radiation.uv_ebeam_exposure import UVEbeamExposureModel
        except ImportError:
            print("❌ Radiation analysis modules not available")
            return None
        
        print("\n☢️  Radiation Damage Analysis")
        print("-" * 40)
        
        if exposures is None:
            exposures = [
                {'type': 'uv', 'wavelength': 185, 'fluence': 1e15},
                {'type': 'uv', 'wavelength': 254, 'fluence': 1e15},
                {'type': 'ebeam', 'energy_keV': 10, 'dose': 1e16},
                {'type': 'ebeam', 'energy_keV': 50, 'dose': 1e15}
            ]
        
        # Progress bar for radiation simulations
        with tqdm(total=len(exposures), desc="☢️  Radiation Tests", 
                 unit="tests", colour='red') as pbar:
            
            results = []
            
            for i, exposure in enumerate(exposures):
                if exposure['type'] == 'uv':
                    pbar.set_description(f"☢️  UV {exposure['wavelength']}nm")
                    
                    sim = RadiationDamageSimulator(self.structure_file)
                    damage = sim.uv_photolysis_simulation(
                        wavelength=exposure['wavelength'],
                        fluence=exposure['fluence']
                    )
                    
                    result = {
                        'type': 'uv',
                        'wavelength': exposure['wavelength'],
                        'fluence': exposure['fluence'],
                        'damage_events': len(damage),
                        'simulator': sim
                    }
                    
                elif exposure['type'] == 'ebeam':
                    pbar.set_description(f"☢️  E-beam {exposure['energy_keV']}keV")
                    
                    sim = RadiationDamageSimulator(self.structure_file)
                    damage = sim.electron_beam_simulation(
                        energy_keV=exposure['energy_keV'],
                        dose=exposure['dose']
                    )
                    
                    result = {
                        'type': 'ebeam',
                        'energy_keV': exposure['energy_keV'],
                        'dose': exposure['dose'],
                        'damage_events': damage['damage_events'],
                        'simulator': sim
                    }
                
                results.append(result)
                pbar.update(1)
                time.sleep(0.3)  # Small delay
        
        # Generate comprehensive radiation report
        print("\n📄 Generating radiation damage report...")
        self._generate_radiation_report(results)
        
        return results
    
    def _generate_radiation_report(self, results):
        """Generate comprehensive radiation damage report"""
        
        report_file = f"{self.name.lower()}_radiation_analysis.txt"
        
        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("        RADIATION DAMAGE ANALYSIS REPORT\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"Structure: {self.structure_file}\n")
            f.write(f"Analysis: {self.name}\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("RADIATION EXPOSURE SUMMARY:\n")
            f.write("-" * 40 + "\n")
            
            for i, result in enumerate(results):
                f.write(f"\nExposure {i+1}:\n")
                if result['type'] == 'uv':
                    f.write(f"  Type: UV Photolysis\n")
                    f.write(f"  Wavelength: {result['wavelength']} nm\n")
                    f.write(f"  Fluence: {result['fluence']:.2e} photons/cm²\n")
                    f.write(f"  Damage Events: {result['damage_events']}\n")
                elif result['type'] == 'ebeam':
                    f.write(f"  Type: Electron Beam\n")
                    f.write(f"  Energy: {result['energy_keV']} keV\n")
                    f.write(f"  Dose: {result['dose']:.2e} e⁻/cm²\n")
                    f.write(f"  Damage Events: {result['damage_events']}\n")
            
            f.write(f"\n\nDetailed simulation reports generated by individual simulators.\n")
            f.write("=" * 80 + "\n")
        
        print(f"📄 Radiation analysis report: {report_file}")
    
    def run_complete_analysis(self, include_radiation=True):
        """Run complete analysis suite with progress tracking"""
        
        print(f"\n🚀 Complete Analysis Suite: {self.name}")
        print("=" * 60)
        
        # Overall progress tracker
        total_steps = 2 if include_radiation else 1
        
        with tqdm(total=100, desc="🔬 Complete Analysis", 
                 bar_format='{l_bar}{bar}| {n:.0f}% [{elapsed}<{remaining}]',
                 colour='green') as overall_pbar:
            
            results = {}
            
            # Bulk analysis
            overall_pbar.set_description("🔬 Bulk structure analysis")
            bulk_results = self.run_bulk_analysis_with_progress()
            results['bulk'] = bulk_results
            overall_pbar.update(50 if include_radiation else 100)
            
            # Radiation analysis (optional)
            if include_radiation:
                overall_pbar.set_description("🔬 Radiation damage analysis")
                radiation_results = self.run_radiation_analysis_with_progress()
                results['radiation'] = radiation_results
                overall_pbar.update(50)
            
            overall_pbar.set_description("✅ Analysis complete")
        
        # Generate summary report
        self._generate_summary_report(results)
        
        total_time = time.time() - self.start_time
        print(f"\n✅ Complete analysis finished in {total_time:.2f} seconds")
        
        return results
    
    def _generate_summary_report(self, results):
        """Generate summary report combining all analyses"""
        
        summary_file = f"{self.name.lower()}_complete_analysis.txt"
        
        with open(summary_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("        COMPLETE MLD ANALYSIS SUMMARY\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"Structure: {self.structure_file}\n")
            f.write(f"Analysis: {self.name}\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Runtime: {time.time() - self.start_time:.2f} seconds\n\n")
            
            # Bulk analysis summary
            if 'bulk' in results and results['bulk']:
                bulk = results['bulk']
                f.write("BULK STRUCTURE PROPERTIES:\n")
                f.write("-" * 40 + "\n")
                
                if 'density' in bulk:
                    density = bulk['density']
                    f.write(f"Density: {density.get('density', 'N/A'):.3f} g/cm³\n")
                    f.write(f"Porosity: {density.get('porosity', 'N/A'):.1%}\n")
                
                if 'layer_info' in bulk:
                    layer = bulk['layer_info']
                    f.write(f"Total thickness: {layer.get('total_thickness', 'N/A'):.2f} Å\n")
                    f.write(f"Number of layers: {layer.get('num_layers', 'N/A')}\n")
                
                f.write(f"Detailed report: {bulk.get('report_file', 'N/A')}\n\n")
            
            # Radiation analysis summary
            if 'radiation' in results and results['radiation']:
                radiation = results['radiation']
                f.write("RADIATION STABILITY:\n")
                f.write("-" * 40 + "\n")
                
                total_damage = sum(r['damage_events'] for r in radiation)
                f.write(f"Total radiation tests: {len(radiation)}\n")
                f.write(f"Total damage events: {total_damage}\n")
                
                uv_tests = [r for r in radiation if r['type'] == 'uv']
                ebeam_tests = [r for r in radiation if r['type'] == 'ebeam']
                
                if uv_tests:
                    f.write(f"UV damage events: {sum(r['damage_events'] for r in uv_tests)}\n")
                if ebeam_tests:
                    f.write(f"E-beam damage events: {sum(r['damage_events'] for r in ebeam_tests)}\n")
                
                f.write(f"Detailed report: {self.name.lower()}_radiation_analysis.txt\n\n")
            
            f.write("FILES GENERATED:\n")
            f.write("-" * 40 + "\n")
            f.write(f"• {summary_file} (this file)\n")
            if 'bulk' in results and results['bulk']:
                f.write(f"• {results['bulk'].get('report_file', 'bulk_analysis.txt')}\n")
            if 'radiation' in results:
                f.write(f"• {self.name.lower()}_radiation_analysis.txt\n")
            
            f.write("\n" + "=" * 80 + "\n")
        
        print(f"📄 Complete analysis summary: {summary_file}")

def analyze_structure_with_progress(structure_file, name=None, include_radiation=True):
    """Convenience function to analyze a structure with progress tracking"""
    
    if name is None:
        name = os.path.splitext(os.path.basename(structure_file))[0]
    
    analyzer = ProgressAnalysis(structure_file, name)
    return analyzer.run_complete_analysis(include_radiation)

def main():
    """Main function for command line usage"""
    
    import argparse
    
    parser = argparse.ArgumentParser(description='Enhanced MLD Analysis with Progress Tracking')
    parser.add_argument('structure', help='Structure file to analyze')
    parser.add_argument('--name', help='Analysis name (default: filename)')
    parser.add_argument('--no-radiation', action='store_true', 
                       help='Skip radiation damage analysis')
    parser.add_argument('--bulk-only', action='store_true',
                       help='Run only bulk structure analysis')
    parser.add_argument('--radiation-only', action='store_true',
                       help='Run only radiation damage analysis')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.structure):
        print(f"❌ Structure file not found: {args.structure}")
        return
    
    name = args.name or os.path.splitext(os.path.basename(args.structure))[0]
    
    try:
        analyzer = ProgressAnalysis(args.structure, name)
        
        if args.bulk_only:
            analyzer.run_bulk_analysis_with_progress()
        elif args.radiation_only:
            analyzer.run_radiation_analysis_with_progress()
        else:
            include_radiation = not args.no_radiation
            analyzer.run_complete_analysis(include_radiation)
            
    except Exception as e:
        print(f"❌ Analysis failed: {e}")

if __name__ == "__main__":
    main()