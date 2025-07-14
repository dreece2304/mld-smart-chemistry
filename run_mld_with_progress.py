#!/usr/bin/env python3
"""
Complete MLD simulation workflow with comprehensive progress tracking
Runs TMA + 2-butyne-1,4-diol deposition with real-time monitoring
"""

import os
import sys
import time
import argparse
import numpy as np
from datetime import datetime, timedelta
from ase.io import read, write
from ase.build import surface
from ase.constraints import FixAtoms
from gpaw import GPAW, PW
from ase.parallel import world

# Import our progress tracking module
from progress_optimization import optimize_structure, GPAWOptimizationMonitor

try:
    from tqdm.auto import tqdm
except ImportError:
    os.system("pip install tqdm")
    from tqdm.auto import tqdm

class MLDWorkflowManager:
    """Manages complete MLD workflow with progress tracking"""
    
    def __init__(self, precision='fast', output_dir='mld_results'):
        self.precision = precision
        self.output_dir = output_dir
        self.start_time = time.time()
        self.stage_times = {}
        self.structures = {}
        self.energies = {}
        
        # Create output directory
        if world.rank == 0:
            os.makedirs(output_dir, exist_ok=True)
            os.chdir(output_dir)
        
        # Set calculation parameters based on precision
        self.calc_params = self._get_calc_params()
        
        # Initialize overall progress tracking
        if world.rank == 0:
            self.overall_progress = tqdm(
                total=100,
                desc="🧪 MLD Workflow Progress",
                bar_format='{l_bar}{bar}| {n:.0f}% [{elapsed}<{remaining}]',
                position=0,
                leave=True,
                colour='magenta'
            )
            
            self.stage_progress = tqdm(
                total=6,  # Number of main stages
                desc="📋 Current Stage",
                unit="stages",
                position=1,
                leave=True,
                colour='cyan'
            )
            
            print(f"\n🚀 MLD Workflow Manager Initialized")
            print(f"📊 Precision: {precision}")
            print(f"📁 Output Directory: {output_dir}")
            print(f"🕒 Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            print(f"{'='*60}\n")
    
    def _get_calc_params(self):
        """Get calculation parameters based on precision level"""
        
        if self.precision == 'fast':
            return {
                'mode': PW(200),
                'xc': 'PBE',
                'kpts': (1, 1, 1),
                'convergence': {'energy': 0.01},
                'maxiter': 100
            }
        elif self.precision == 'medium':
            return {
                'mode': PW(300),
                'xc': 'PBE',
                'kpts': (2, 2, 1),
                'convergence': {'energy': 0.005},
                'maxiter': 200
            }
        elif self.precision == 'production':
            return {
                'mode': PW(400),
                'xc': 'PBE',
                'kpts': (4, 4, 1),
                'convergence': {'energy': 0.0005},
                'maxiter': 300
            }
        else:
            raise ValueError(f"Unknown precision: {self.precision}")
    
    def _update_overall_progress(self, percentage, stage_name):
        """Update overall workflow progress"""
        if world.rank == 0:
            self.overall_progress.n = percentage
            self.overall_progress.set_description(f"🧪 MLD Workflow: {stage_name}")
            self.overall_progress.refresh()
    
    def _advance_stage(self, stage_name):
        """Advance to next stage"""
        if world.rank == 0:
            self.stage_progress.update(1)
            self.stage_progress.set_postfix_str(stage_name)
    
    def stage_1_preparation(self):
        """Stage 1: Structure preparation and loading"""
        stage_start = time.time()
        
        if world.rank == 0:
            print("📋 Stage 1: Structure Preparation")
            print("-" * 40)
        
        self._update_overall_progress(0, "Structure Preparation")
        
        # Load or create molecular structures
        try:
            self.tma = read('../structures/tma_molecule.xyz')
            if world.rank == 0:
                print(f"✅ Loaded TMA: {len(self.tma)} atoms")
        except:
            if world.rank == 0:
                print("⚠️  Creating simple TMA structure")
            self.tma = self._create_simple_tma()
        
        try:
            self.diol = read('../structures/butyne_diol_molecule.xyz')
            if world.rank == 0:
                print(f"✅ Loaded Diol: {len(self.diol)} atoms")
        except:
            if world.rank == 0:
                print("⚠️  Creating simple diol structure")
            self.diol = self._create_simple_diol()
        
        # Create surface
        if world.rank == 0:
            print("🏗️  Creating hydroxylated surface...")
        
        self.surface = self._create_hydroxylated_surface()
        
        if world.rank == 0:
            print(f"✅ Surface created: {len(self.surface)} atoms")
        
        self.stage_times['preparation'] = time.time() - stage_start
        self._advance_stage("Preparation Complete")
        self._update_overall_progress(10, "Structures Ready")
        
        if world.rank == 0:
            print(f"⏱️  Stage 1 completed in {self.stage_times['preparation']:.1f}s\n")
    
    def stage_2_optimize_molecules(self):
        """Stage 2: Optimize individual molecules"""
        stage_start = time.time()
        
        if world.rank == 0:
            print("🔧 Stage 2: Molecule Optimization")
            print("-" * 40)
        
        self._update_overall_progress(20, "Optimizing TMA")
        
        # Optimize TMA
        self.tma.center(vacuum=8.0)
        
        calc_params_mol = self.calc_params.copy()
        calc_params_mol['kpts'] = (1, 1, 1)  # Gamma point for molecules
        
        self.tma_opt, tma_stats = optimize_structure(
            self.tma, calc_params_mol,
            optimizer='BFGS', fmax=0.05, max_steps=100,
            name="TMA_Molecule"
        )
        
        self.structures['tma'] = self.tma_opt
        self.energies['tma'] = tma_stats['final_energy']
        
        self._update_overall_progress(35, "Optimizing Diol")
        
        # Optimize Diol
        self.diol.center(vacuum=8.0)
        
        self.diol_opt, diol_stats = optimize_structure(
            self.diol, calc_params_mol,
            optimizer='BFGS', fmax=0.05, max_steps=100,
            name="Diol_Molecule"
        )
        
        self.structures['diol'] = self.diol_opt
        self.energies['diol'] = diol_stats['final_energy']
        
        self.stage_times['molecules'] = time.time() - stage_start
        self._advance_stage("Molecules Optimized")
        self._update_overall_progress(50, "Molecules Ready")
        
        if world.rank == 0:
            print(f"⏱️  Stage 2 completed in {self.stage_times['molecules']:.1f}s\n")
    
    def stage_3_optimize_surface(self):
        """Stage 3: Optimize surface structure"""
        stage_start = time.time()
        
        if world.rank == 0:
            print("🏔️  Stage 3: Surface Optimization")
            print("-" * 40)
        
        self._update_overall_progress(55, "Optimizing Surface")
        
        # Fix bottom layers
        positions = self.surface.get_positions()
        z_min = positions[:, 2].min()
        
        fixed_indices = []
        for i, pos in enumerate(positions):
            if pos[2] < z_min + 3.0:  # Fix bottom 3 Å
                fixed_indices.append(i)
        
        self.surface.set_constraint(FixAtoms(indices=fixed_indices))
        
        # Optimize surface
        calc_params_surf = self.calc_params.copy()
        
        self.surface_opt, surf_stats = optimize_structure(
            self.surface, calc_params_surf,
            optimizer='BFGS', fmax=0.05, max_steps=150,
            name="OH_Surface"
        )
        
        self.structures['surface'] = self.surface_opt
        self.energies['surface'] = surf_stats['final_energy']
        
        self.stage_times['surface'] = time.time() - stage_start
        self._advance_stage("Surface Optimized")
        self._update_overall_progress(65, "Surface Ready")
        
        if world.rank == 0:
            print(f"⏱️  Stage 3 completed in {self.stage_times['surface']:.1f}s\n")
    
    def stage_4_tma_adsorption(self):
        """Stage 4: Simulate TMA adsorption"""
        stage_start = time.time()
        
        if world.rank == 0:
            print("🔗 Stage 4: TMA Adsorption Simulation")
            print("-" * 40)
        
        self._update_overall_progress(70, "TMA Adsorption")
        
        # Create TMA + surface system
        combined = self.surface_opt.copy()
        
        # Position TMA above surface
        surf_positions = combined.get_positions()
        surf_center = [surf_positions[:, 0].mean(), surf_positions[:, 1].mean()]
        surf_top = surf_positions[:, 2].max()
        
        tma_positioned = self.tma_opt.copy()
        tma_positioned.translate([surf_center[0], surf_center[1], surf_top + 4.0] - 
                               tma_positioned.get_center_of_mass())
        
        # Combine structures
        combined.extend(tma_positioned)
        
        # Remove surface constraints and add new ones
        combined.set_constraint(FixAtoms(indices=list(range(len(self.surface_opt)//2))))
        
        # Optimize TMA + surface
        self.tma_surface, tma_surf_stats = optimize_structure(
            combined, self.calc_params,
            optimizer='BFGS', fmax=0.05, max_steps=100,
            name="TMA_Adsorption"
        )
        
        self.structures['tma_surface'] = self.tma_surface
        self.energies['tma_surface'] = tma_surf_stats['final_energy']
        
        # Calculate adsorption energy
        self.adsorption_energy = (self.energies['tma_surface'] - 
                                self.energies['surface'] - self.energies['tma'])
        
        if world.rank == 0:
            print(f"🧮 TMA Adsorption Energy: {self.adsorption_energy:.6f} eV")
        
        self.stage_times['tma_adsorption'] = time.time() - stage_start
        self._advance_stage("TMA Adsorbed")
        self._update_overall_progress(80, "TMA Adsorption Complete")
        
        if world.rank == 0:
            print(f"⏱️  Stage 4 completed in {self.stage_times['tma_adsorption']:.1f}s\n")
    
    def stage_5_diol_addition(self):
        """Stage 5: Simulate diol addition"""
        stage_start = time.time()
        
        if world.rank == 0:
            print("🔬 Stage 5: Diol Addition Simulation")
            print("-" * 40)
        
        self._update_overall_progress(85, "Diol Addition")
        
        # Create diol + TMA-surface system
        combined = self.tma_surface.copy()
        
        # Position diol above TMA
        surf_positions = combined.get_positions()
        surf_center = [surf_positions[:, 0].mean(), surf_positions[:, 1].mean()]
        surf_top = surf_positions[:, 2].max()
        
        diol_positioned = self.diol_opt.copy()
        diol_positioned.translate([surf_center[0] + 2.0, surf_center[1], surf_top + 3.0] - 
                                diol_positioned.get_center_of_mass())
        
        # Combine structures
        combined.extend(diol_positioned)
        
        # Keep same constraints
        combined.set_constraint(FixAtoms(indices=list(range(len(self.surface_opt)//2))))
        
        # Optimize final system
        self.final_system, final_stats = optimize_structure(
            combined, self.calc_params,
            optimizer='BFGS', fmax=0.05, max_steps=120,
            name="Final_MLD_System"
        )
        
        self.structures['final'] = self.final_system
        self.energies['final'] = final_stats['final_energy']
        
        # Calculate total reaction energy
        self.reaction_energy = (self.energies['final'] - 
                              self.energies['surface'] - 
                              self.energies['tma'] - self.energies['diol'])
        
        if world.rank == 0:
            print(f"🧮 Total Reaction Energy: {self.reaction_energy:.6f} eV")
        
        self.stage_times['diol_addition'] = time.time() - stage_start
        self._advance_stage("Diol Added")
        self._update_overall_progress(95, "MLD Cycle Complete")
        
        if world.rank == 0:
            print(f"⏱️  Stage 5 completed in {self.stage_times['diol_addition']:.1f}s\n")
    
    def stage_6_analysis(self):
        """Stage 6: Analysis and reporting"""
        stage_start = time.time()
        
        if world.rank == 0:
            print("📊 Stage 6: Analysis and Reporting")
            print("-" * 40)
        
        self._update_overall_progress(98, "Analysis")
        
        # Perform bulk analysis if available
        try:
            from analysis.bulk_structure_analysis import BulkMLDAnalyzer
            
            if world.rank == 0:
                print("🔍 Running bulk structure analysis...")
            
            analyzer = BulkMLDAnalyzer('final_mld_system_optimized.xyz')
            analyzer.generate_report('mld_analysis_report.txt')
            
            if world.rank == 0:
                print("✅ Bulk analysis complete")
                
        except Exception as e:
            if world.rank == 0:
                print(f"⚠️  Bulk analysis failed: {e}")
        
        # Generate comprehensive report
        self._generate_final_report()
        
        self.stage_times['analysis'] = time.time() - stage_start
        self._advance_stage("Analysis Complete")
        self._update_overall_progress(100, "Complete!")
        
        if world.rank == 0:
            print(f"⏱️  Stage 6 completed in {self.stage_times['analysis']:.1f}s\n")
    
    def _generate_final_report(self):
        """Generate comprehensive final report"""
        if world.rank != 0:
            return
            
        total_time = time.time() - self.start_time
        
        report = f"""
{'='*80}
        MLD SIMULATION WORKFLOW REPORT
{'='*80}

📅 Simulation Details:
   Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
   Precision: {self.precision}
   Total Time: {total_time:.2f} s ({total_time/60:.2f} min)
   
🧪 Chemical System:
   Precursor 1: Trimethylaluminum (TMA)
   Precursor 2: 2-butyne-1,4-diol
   Substrate: Hydroxylated Si(100) surface
   
📊 Results Summary:
   TMA Energy: {self.energies.get('tma', 'N/A'):.6f} eV
   Diol Energy: {self.energies.get('diol', 'N/A'):.6f} eV
   Surface Energy: {self.energies.get('surface', 'N/A'):.6f} eV
   TMA Adsorption Energy: {getattr(self, 'adsorption_energy', 'N/A'):.6f} eV
   Total Reaction Energy: {getattr(self, 'reaction_energy', 'N/A'):.6f} eV
   
⏱️  Timing Breakdown:
"""
        
        for stage, time_val in self.stage_times.items():
            percentage = (time_val / total_time) * 100
            report += f"   {stage.replace('_', ' ').title()}: {time_val:.2f} s ({percentage:.1f}%)\n"
        
        report += f"""
📁 Output Files:
   Optimized Structures: *_optimized.xyz
   Trajectories: *_opt.traj
   GPAW Output: *_gpaw.out
   Analysis Report: mld_analysis_report.txt
   
🔬 Next Steps:
   1. Examine optimized structures with visualization software
   2. Run radiation damage simulations
   3. Study bulk properties and film growth
   4. Analyze reaction mechanisms and pathways
   
{'='*80}
"""
        
        with open('mld_workflow_report.txt', 'w') as f:
            f.write(report)
        
        print(report)
        print("📄 Full report saved: mld_workflow_report.txt")
    
    def _create_simple_tma(self):
        """Create simple TMA structure if file not found"""
        from ase import Atoms
        
        tma = Atoms('Al', positions=[[0, 0, 0]])
        # Add methyl groups (simplified)
        positions = [[0, 0, 2], [1.7, 0, -1], [-1.7, 0, -1]]
        for pos in positions:
            tma.extend(Atoms('CH3', positions=[pos]))
        
        return tma
    
    def _create_simple_diol(self):
        """Create simple diol structure if file not found"""
        from ase import Atoms
        
        # Create acetylene backbone with OH groups
        diol = Atoms('C4H6O2', positions=[
            [-0.6, 0, 0], [0.6, 0, 0],      # C≡C
            [-2.0, 0, 0], [2.0, 0, 0],      # CH2
            [-3.0, 0, 0], [3.0, 0, 0]       # OH
        ])
        
        return diol
    
    def _create_hydroxylated_surface(self):
        """Create hydroxylated Si(100) surface"""
        
        # Create Si surface
        if self.precision == 'fast':
            surf = surface('Si', (1, 0, 0), 3)  # 3 layers
            surf = surf.repeat((2, 2, 1))       # 2x2
        else:
            surf = surface('Si', (1, 0, 0), 4)  # 4 layers
            surf = surf.repeat((3, 3, 1))       # 3x3
        
        surf.center(vacuum=10.0, axis=2)
        
        # Add OH groups
        from ase import Atoms
        positions = surf.get_positions()
        top_z = positions[:, 2].max()
        
        oh_atoms = Atoms()
        for i, (pos, symbol) in enumerate(zip(positions, surf.get_chemical_symbols())):
            if symbol == 'Si' and abs(pos[2] - top_z) < 0.5:
                # Add OH group every other Si
                if i % 2 == 0:
                    o_pos = pos + [0, 0, 1.6]
                    h_pos = o_pos + [0, 0, 0.96]
                    oh_atoms.extend(Atoms('OH', positions=[o_pos, h_pos]))
        
        surf.extend(oh_atoms)
        return surf
    
    def run_complete_workflow(self):
        """Run the complete MLD workflow"""
        
        try:
            self.stage_1_preparation()
            self.stage_2_optimize_molecules()
            self.stage_3_optimize_surface()
            self.stage_4_tma_adsorption()
            self.stage_5_diol_addition()
            self.stage_6_analysis()
            
            if world.rank == 0:
                self.overall_progress.close()
                self.stage_progress.close()
                
                print(f"\n🎉 MLD Workflow Complete!")
                print(f"⏱️  Total Time: {time.time() - self.start_time:.2f} seconds")
                print(f"📁 Results saved in: {self.output_dir}")
                
        except KeyboardInterrupt:
            if world.rank == 0:
                print(f"\n⚠️  Workflow interrupted by user")
                self.overall_progress.close()
                self.stage_progress.close()
        except Exception as e:
            if world.rank == 0:
                print(f"\n❌ Workflow failed: {e}")
                self.overall_progress.close()
                self.stage_progress.close()
            raise

def main():
    """Main function with command line interface"""
    
    parser = argparse.ArgumentParser(description='MLD Simulation with Progress Tracking')
    parser.add_argument('--precision', choices=['fast', 'medium', 'production'], 
                       default='fast', help='Calculation precision level')
    parser.add_argument('--output', default='mld_results', 
                       help='Output directory')
    parser.add_argument('--molecule', choices=['tma', 'diol'], 
                       help='Optimize single molecule only')
    parser.add_argument('--test', action='store_true', 
                       help='Run quick test optimization')
    
    args = parser.parse_args()
    
    # Set GPAW environment
    gpaw_path = os.path.expanduser('~/.local/share/gpaw/gpaw-setups-24.11.0')
    os.environ['GPAW_SETUP_PATH'] = gpaw_path
    
    if args.test:
        # Run quick test
        from progress_optimization import quick_test_optimization
        quick_test_optimization()
        return
    
    if args.molecule:
        # Run single molecule optimization
        workflow = MLDWorkflowManager(args.precision, args.output)
        workflow.stage_1_preparation()
        
        if args.molecule == 'tma':
            workflow.stage_2_optimize_molecules()  # Just TMA part
        elif args.molecule == 'diol':
            workflow.stage_2_optimize_molecules()  # Just diol part
    else:
        # Run complete workflow
        workflow = MLDWorkflowManager(args.precision, args.output)
        workflow.run_complete_workflow()

if __name__ == "__main__":
    main()