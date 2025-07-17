#!/usr/bin/env python3
"""
Quick MLD execution script with user-friendly interface
Provides simple commands for common MLD simulation tasks
"""

import os
import sys
import argparse
from datetime import datetime
from pathlib import Path

def setup_environment():
    """Enhanced GPAW environment setup and dependency validation"""
    
    print("🔍 Validating MLD environment...")
    
    # Check conda environment
    conda_env = os.environ.get('CONDA_DEFAULT_ENV', 'unknown')
    if conda_env == 'mld_modeling':
        print("✅ In mld_modeling environment")
    else:
        print(f"⚠️  Not in mld_modeling environment (current: {conda_env})")
        print("   Recommend: conda activate mld_modeling")
    
    # Set GPAW path
    gpaw_path = os.path.expanduser('~/.local/share/gpaw/gpaw-setups-24.11.0')
    os.environ['GPAW_SETUP_PATH'] = gpaw_path
    
    # Check if GPAW datasets exist
    if not os.path.exists(gpaw_path):
        print("❌ GPAW datasets not found!")
        print("Solutions:")
        print("  1. Run: gpaw install-data ~/.local/share/gpaw")
        print("  2. Run: python validate_setup.py --fix-issues")
        print("  3. Run: ./setup_environment.sh")
        return False
    else:
        # Count dataset files for validation
        dataset_files = list(Path(gpaw_path).glob('*.py'))
        print(f"✅ GPAW datasets found ({len(dataset_files)} files)")
    
    # Check core packages with versions
    core_packages = {
        'ase': '3.22.0',
        'gpaw': '22.8.0',
        'numpy': '1.21.0',
        'scipy': '1.7.0',
        'tqdm': '4.0.0'
    }
    
    missing_packages = []
    for package, min_version in core_packages.items():
        try:
            module = __import__(package)
            version = getattr(module, '__version__', 'unknown')
            print(f"✅ {package}: {version}")
        except ImportError:
            missing_packages.append(package)
            print(f"❌ {package}: not installed")
    
    # Handle missing packages
    if missing_packages:
        print(f"\n❌ Missing packages: {', '.join(missing_packages)}")
        print("Solutions:")
        print("  1. Run: python validate_setup.py --fix-issues")
        print("  2. Run: ./setup_environment.sh")
        print(f"  3. Manual: pip install {' '.join(missing_packages)}")
        return False
    
    # Check structure files
    required_files = [
        'structures/tma_molecule.xyz',
        'structures/butyne_diol_molecule.xyz',
        'structures/simple_si_surface.xyz'
    ]
    
    missing_structures = []
    for file_path in required_files:
        if os.path.exists(file_path):
            size = os.path.getsize(file_path)
            print(f"✅ {file_path} ({size} bytes)")
        else:
            missing_structures.append(file_path)
            print(f"❌ {file_path}: missing")
    
    if missing_structures:
        print(f"\n❌ Missing structure files")
        print("Solution: Check if you're in the correct directory")
        return False
    
    # Check system resources
    try:
        import psutil
        memory_gb = psutil.virtual_memory().total / (1024**3)
        cpu_count = psutil.cpu_count()
        
        if memory_gb >= 4:
            print(f"✅ Memory: {memory_gb:.1f} GB")
        else:
            print(f"⚠️  Low memory: {memory_gb:.1f} GB (recommend 4+ GB)")
            
        if cpu_count >= 2:
            print(f"✅ CPUs: {cpu_count}")
        else:
            print(f"⚠️  Low CPU count: {cpu_count}")
            
    except ImportError:
        print("⚠️  psutil not available for system monitoring")
    
    print("✅ Environment validation completed")
    return True

def run_quick_test():
    """Run quick functionality test"""
    
    print("🧪 Running Quick MLD Test...")
    print("=" * 40)
    
    if not setup_environment():
        return
    
    from progress_optimization import quick_test_optimization
    
    try:
        atoms, stats = quick_test_optimization()
        
        print("\n✅ Quick test completed successfully!")
        print(f"📊 Statistics:")
        for key, value in stats.items():
            print(f"   {key}: {value}")
            
    except Exception as e:
        print(f"❌ Test failed: {e}")
        print("Check your GPAW installation and try again.")

def run_single_molecule(molecule_name, precision='fast'):
    """Optimize a single molecule"""
    
    print(f"🧬 Optimizing {molecule_name.upper()} Molecule...")
    print("=" * 40)
    
    if not setup_environment():
        return
    
    from run_mld_with_progress import MLDWorkflowManager
    
    try:
        # Create workflow manager
        output_dir = f"{molecule_name}_{precision}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        workflow = MLDWorkflowManager(precision, output_dir)
        
        # Run preparation and molecule optimization
        workflow.stage_1_preparation()
        
        if molecule_name.lower() == 'tma':
            workflow.tma.center(vacuum=8.0)
            calc_params = workflow.calc_params.copy()
            calc_params['kpts'] = (1, 1, 1)
            
            from progress_optimization import optimize_structure
            atoms, stats = optimize_structure(
                workflow.tma, calc_params,
                optimizer='BFGS', fmax=0.05, max_steps=100,
                name="TMA_Single"
            )
            
        elif molecule_name.lower() == 'diol':
            # Check if we need to use improved diol structure
            improved_diol_file = 'structures/diol_improved.xyz'
            if os.path.exists(improved_diol_file):
                print("🔧 Using improved diol structure...")
                from ase.io import read
                workflow.diol = read(improved_diol_file)
            
            workflow.diol.center(vacuum=8.0)
            calc_params = workflow.calc_params.copy()
            calc_params['kpts'] = (1, 1, 1)
            
            from progress_optimization import optimize_structure
            atoms, stats = optimize_structure(
                workflow.diol, calc_params,
                optimizer='BFGS', fmax=0.05, max_steps=100,
                name="Diol_Single"
            )
        
        print(f"\n✅ {molecule_name.upper()} optimization completed!")
        print(f"📁 Results in: {output_dir}")
        
    except Exception as e:
        print(f"❌ Optimization failed: {e}")

def run_surface_only(precision='fast'):
    """Optimize surface structure only"""
    
    print("🏔️  Optimizing Surface Structure...")
    print("=" * 40)
    
    if not setup_environment():
        return
    
    from run_mld_with_progress import MLDWorkflowManager
    
    try:
        output_dir = f"surface_{precision}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        workflow = MLDWorkflowManager(precision, output_dir)
        
        workflow.stage_1_preparation()
        workflow.stage_3_optimize_surface()
        
        print(f"\n✅ Surface optimization completed!")
        print(f"📁 Results in: {output_dir}")
        
    except Exception as e:
        print(f"❌ Surface optimization failed: {e}")

def run_full_cycle(precision='fast'):
    """Run complete MLD cycle"""
    
    print("🔬 Running Complete MLD Cycle...")
    print("=" * 40)
    
    if not setup_environment():
        return
    
    from run_mld_with_progress import MLDWorkflowManager
    
    try:
        output_dir = f"mld_full_{precision}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        workflow = MLDWorkflowManager(precision, output_dir)
        workflow.run_complete_workflow()
        
        print(f"\n✅ Complete MLD cycle finished!")
        print(f"📁 Results in: {output_dir}")
        
    except Exception as e:
        print(f"❌ MLD cycle failed: {e}")

def run_thermal_simulation(precision='fast'):
    """Run thermal MLD simulation with realistic temperature conditions"""
    
    print("🌡️  Running Thermal MLD Simulation...")
    print("=" * 40)
    
    if not setup_environment():
        return
    
    try:
        import subprocess
        import sys
        
        cmd = [
            sys.executable, 'thermal_mld_simulation.py',
            '--precision', precision,
            '--output', f'thermal_mld_{precision}_{datetime.now().strftime("%Y%m%d_%H%M%S")}'
        ]
        
        print("🔥 Starting thermal simulation with simplified conditions:")
        print("   - TMA: 25°C → 120°C (heats up in chamber)")
        print("   - Diol: 80°C → 120°C (vapor equilibrates in chamber)")
        print("   - Surface: 120°C (chamber temperature)")
        print("   - Final deposition: 120°C (everything at chamber temp)")
        print()
        
        result = subprocess.run(cmd, capture_output=False, text=True)
        
        if result.returncode == 0:
            print(f"\n✅ Thermal MLD simulation completed!")
        else:
            print(f"\n❌ Thermal simulation failed with return code {result.returncode}")
            
    except Exception as e:
        print(f"❌ Failed to run thermal simulation: {e}")

def run_radiation_test():
    """Run radiation damage simulation test"""
    
    print("☢️  Running Radiation Damage Simulation...")
    print("=" * 40)
    
    try:
        from calculations.radiation.radiation_damage import RadiationDamageSimulator
        
        # Use any available structure file
        structure_files = [
            'structures/tma_molecule.xyz',
            'structures/butyne_diol_molecule.xyz',
            'structures/simple_si_surface.xyz'
        ]
        
        structure_file = None
        for fname in structure_files:
            if os.path.exists(fname):
                structure_file = fname
                break
        
        if not structure_file:
            print("❌ No structure files found for radiation test")
            return
        
        print(f"📁 Using structure: {structure_file}")
        
        sim = RadiationDamageSimulator(structure_file)
        
        # UV test
        print("🌞 Testing UV photolysis (254 nm)...")
        uv_damage = sim.uv_photolysis_simulation(wavelength=254, fluence=1e15)
        print(f"   UV damage events: {len(uv_damage)}")
        
        # E-beam test
        print("⚡ Testing electron beam (10 keV)...")
        eb_damage = sim.electron_beam_simulation(energy_keV=10, dose=1e16)
        print(f"   E-beam damage events: {eb_damage['damage_events']}")
        
        sim.generate_damage_report('radiation_test_report.txt')
        
        print(f"\n✅ Radiation test completed!")
        print(f"📄 Report saved: radiation_test_report.txt")
        
    except Exception as e:
        print(f"❌ Radiation test failed: {e}")

def run_bulk_analysis():
    """Run bulk structure analysis"""
    
    print("📊 Running Bulk Structure Analysis...")
    print("=" * 40)
    
    try:
        from analysis.bulk_structure_analysis import BulkMLDAnalyzer
        
        # Find the most recent MLD result
        structure_files = []
        for root, dirs, files in os.walk('.'):
            for file in files:
                if file.endswith('_optimized.xyz') and 'mld' in file.lower():
                    structure_files.append(os.path.join(root, file))
        
        if not structure_files:
            print("❌ No optimized MLD structures found")
            print("Run a full MLD cycle first with: python quick_run.py --full-cycle")
            return
        
        # Use the most recent file
        structure_file = max(structure_files, key=os.path.getmtime)
        print(f"📁 Analyzing: {structure_file}")
        
        analyzer = BulkMLDAnalyzer(structure_file)
        analyzer.generate_report('bulk_analysis_report.txt')
        
        print(f"\n✅ Bulk analysis completed!")
        print(f"📄 Report saved: bulk_analysis_report.txt")
        
    except Exception as e:
        print(f"❌ Bulk analysis failed: {e}")

def run_diol_fix():
    """Run diol structure fix utility"""
    
    print("🔧 Running Diol Structure Fix...")
    print("=" * 40)
    
    try:
        # Import and run the fix utility
        import subprocess
        import sys
        
        result = subprocess.run([sys.executable, 'fix_diol_structure.py'], 
                              capture_output=True, text=True)
        
        print(result.stdout)
        if result.stderr:
            print("Errors:")
            print(result.stderr)
        
        if result.returncode == 0:
            print("\n✅ Diol structure fix completed!")
            print("Now try: python quick_run.py --molecule diol")
        else:
            print(f"\n❌ Diol fix failed with return code {result.returncode}")
            
    except Exception as e:
        print(f"❌ Failed to run diol fix: {e}")

def show_status():
    """Show current status and available results"""
    
    print("📋 MLD Project Status")
    print("=" * 40)
    
    # Check installation
    print("🔧 Installation Status:")
    
    try:
        import ase, gpaw, pymatgen
        print("   ✅ ASE, GPAW, PyMatGen installed")
    except ImportError as e:
        print(f"   ❌ Missing packages: {e}")
    
    gpaw_path = os.path.expanduser('~/.local/share/gpaw/gpaw-setups-24.11.0')
    if os.path.exists(gpaw_path):
        print("   ✅ GPAW datasets available")
    else:
        print("   ❌ GPAW datasets missing")
    
    # Check structure files
    print("\n📁 Structure Files:")
    structure_files = [
        'structures/tma_molecule.xyz',
        'structures/butyne_diol_molecule.xyz',
        'structures/simple_si_surface.xyz'
    ]
    
    for fname in structure_files:
        if os.path.exists(fname):
            print(f"   ✅ {fname}")
        else:
            print(f"   ❌ {fname}")
    
    # Check recent results
    print("\n📊 Recent Results:")
    result_dirs = [d for d in os.listdir('.') if os.path.isdir(d) and 
                   ('mld' in d.lower() or 'tma' in d.lower() or 'diol' in d.lower())]
    
    if result_dirs:
        for d in sorted(result_dirs)[-5:]:  # Show last 5
            print(f"   📂 {d}")
    else:
        print("   No recent calculations found")
    
    print(f"\n💡 Quick Commands:")
    print(f"   Test: python quick_run.py --test")
    print(f"   Single molecule: python quick_run.py --molecule tma")
    print(f"   Full cycle: python quick_run.py --full-cycle")
    print(f"   Radiation test: python quick_run.py --radiation")

def main():
    """Main command line interface"""
    
    parser = argparse.ArgumentParser(
        description='Quick MLD Simulation Interface',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python quick_run.py --test                    # Quick functionality test
  python quick_run.py --molecule tma            # Optimize TMA only
  python quick_run.py --molecule diol           # Optimize diol only
  python quick_run.py --surface                 # Optimize surface only
  python quick_run.py --full-cycle              # Complete MLD cycle
  python quick_run.py --thermal                 # Thermal MLD with realistic conditions
  python quick_run.py --thermal --precision medium  # Higher accuracy thermal
  python quick_run.py --radiation               # Radiation damage test
  python quick_run.py --analysis                # Bulk structure analysis
  python quick_run.py --status                  # Show project status
        """
    )
    
    # Main action arguments (mutually exclusive)
    action_group = parser.add_mutually_exclusive_group(required=True)
    action_group.add_argument('--test', action='store_true',
                             help='Run quick functionality test')
    action_group.add_argument('--molecule', choices=['tma', 'diol'],
                             help='Optimize single molecule')
    action_group.add_argument('--surface', action='store_true',
                             help='Optimize surface structure only')
    action_group.add_argument('--full-cycle', action='store_true',
                             help='Run complete MLD cycle')
    action_group.add_argument('--thermal', action='store_true',
                             help='Run thermal MLD simulation with realistic conditions')
    action_group.add_argument('--radiation', action='store_true',
                             help='Run radiation damage simulation')
    action_group.add_argument('--analysis', action='store_true',
                             help='Run bulk structure analysis')
    action_group.add_argument('--status', action='store_true',
                             help='Show project status')
    action_group.add_argument('--fix-diol', action='store_true',
                             help='Fix diol structure for optimization')
    
    # Optional arguments
    parser.add_argument('--precision', choices=['fast', 'medium', 'production'],
                       default='fast', help='Calculation precision (default: fast)')
    
    args = parser.parse_args()
    
    # Print header
    print("\n" + "="*60)
    print("🧪 MLD SIMULATION QUICK INTERFACE")
    print("="*60)
    print(f"🕒 Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"⚙️  Precision: {args.precision}")
    print("="*60)
    
    # Execute requested action
    try:
        if args.test:
            run_quick_test()
        elif args.molecule:
            run_single_molecule(args.molecule, args.precision)
        elif args.surface:
            run_surface_only(args.precision)
        elif args.full_cycle:
            run_full_cycle(args.precision)
        elif args.thermal:
            run_thermal_simulation(args.precision)
        elif args.radiation:
            run_radiation_test()
        elif args.analysis:
            run_bulk_analysis()
        elif args.status:
            show_status()
        elif args.fix_diol:
            run_diol_fix()
            
    except KeyboardInterrupt:
        print("\n⚠️  Operation interrupted by user")
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        print("Use --status to check your installation")
    
    print("\n" + "="*60)
    print("🏁 MLD Quick Interface Complete")
    print("="*60)

if __name__ == "__main__":
    main()