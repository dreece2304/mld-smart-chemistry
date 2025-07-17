#!/usr/bin/env python3
"""
Test current MLD setup capabilities
Works with whatever packages are available
"""

import os
import sys

def test_environment():
    """Test what's available in current environment"""
    
    print("🔍 TESTING CURRENT MLD SETUP")
    print("=" * 50)
    
    # Test basic Python
    print(f"✅ Python {sys.version.split()[0]}")
    print(f"✅ Working directory: {os.getcwd()}")
    
    # Test core packages
    packages_to_test = [
        ('ase', 'Atomic Simulation Environment'),
        ('gpaw', 'GPAW DFT calculator'),
        ('numpy', 'NumPy arrays'),
        ('matplotlib', 'Plotting'),
        ('scipy', 'Scientific computing'),
        ('tqdm', 'Progress bars'),
        ('mpi4py', 'MPI parallel computing'),
        ('psutil', 'System monitoring')
    ]
    
    available_packages = []
    missing_packages = []
    
    for package, description in packages_to_test:
        try:
            module = __import__(package)
            version = getattr(module, '__version__', 'unknown')
            print(f"✅ {package} {version} - {description}")
            available_packages.append(package)
        except ImportError:
            print(f"❌ {package} - {description}")
            missing_packages.append(package)
    
    # Test structure files
    print(f"\n📁 STRUCTURE FILES:")
    structure_files = [
        'structures/tma_molecule.xyz',
        'structures/butyne_diol_molecule.xyz',
        'structures/simple_si_surface.xyz',
        'structures/diol_improved.xyz'
    ]
    
    for file_path in structure_files:
        if os.path.exists(file_path):
            size = os.path.getsize(file_path)
            print(f"✅ {file_path} ({size} bytes)")
        else:
            print(f"❌ {file_path}")
    
    # Test key scripts
    print(f"\n🐍 SCRIPT FILES:")
    script_files = [
        'quick_run.py',
        'thermal_mld_simulation.py',
        'enhanced_monitor.py',
        'visualization_tools.py',
        'validate_setup.py'
    ]
    
    for script in script_files:
        if os.path.exists(script):
            print(f"✅ {script}")
        else:
            print(f"❌ {script}")
    
    # Summary
    print(f"\n📊 SUMMARY:")
    print(f"✅ Available packages: {len(available_packages)}/{len(packages_to_test)}")
    print(f"✅ Core functionality: {'ASE' in [p for p in available_packages]} and {'GPAW' in [p for p in available_packages]}")
    
    if 'ase' in available_packages and 'gpaw' in available_packages:
        print(f"🎉 READY FOR MLD SIMULATIONS!")
        print(f"\nNext steps:")
        print(f"1. Activate environment: conda activate mld_modeling")
        print(f"2. Test: python quick_run.py --test")
        print(f"3. Run: python thermal_mld_simulation.py --molecule tma --precision fast")
    else:
        print(f"⚠️  Need to activate mld_modeling environment first")
        print(f"Run: conda activate mld_modeling")
    
    return len(missing_packages) == 0

def test_visualization_basics():
    """Test basic visualization without complex dependencies"""
    
    print(f"\n🎨 TESTING BASIC VISUALIZATION:")
    
    try:
        import matplotlib.pyplot as plt
        import numpy as np
        
        # Create simple test plot
        x = np.linspace(0, 10, 100)
        y = np.sin(x)
        
        plt.figure(figsize=(8, 4))
        plt.plot(x, y, 'b-', linewidth=2, label='sin(x)')
        plt.title('MLD Visualization Test')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Save plot
        test_plot = 'test_visualization.png'
        plt.savefig(test_plot, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"✅ Created test plot: {test_plot}")
        return True
        
    except ImportError as e:
        print(f"❌ Matplotlib not available: {e}")
        return False
    except Exception as e:
        print(f"❌ Plotting error: {e}")
        return False

def test_structure_reading():
    """Test structure file reading if ASE available"""
    
    print(f"\n🧬 TESTING STRUCTURE READING:")
    
    try:
        from ase.io import read
        
        test_files = [
            'structures/tma_molecule.xyz',
            'structures/butyne_diol_molecule.xyz'
        ]
        
        for file_path in test_files:
            if os.path.exists(file_path):
                try:
                    atoms = read(file_path)
                    symbols = atoms.get_chemical_symbols()
                    print(f"✅ {file_path}: {len(atoms)} atoms ({', '.join(set(symbols))})")
                except Exception as e:
                    print(f"❌ {file_path}: Read error - {e}")
            else:
                print(f"⚠️  {file_path}: File not found")
        
        return True
        
    except ImportError:
        print(f"❌ ASE not available in current environment")
        return False

def main():
    """Main test function"""
    
    env_ok = test_environment()
    viz_ok = test_visualization_basics() 
    struct_ok = test_structure_reading()
    
    print(f"\n" + "=" * 50)
    print(f"🎯 OVERALL STATUS:")
    
    if env_ok and viz_ok and struct_ok:
        print(f"🎉 ALL TESTS PASSED - SYSTEM READY!")
    elif struct_ok:
        print(f"✅ CORE FUNCTIONALITY READY")
        print(f"⚠️  Some optional features missing")
    else:
        print(f"⚠️  NEED TO ACTIVATE mld_modeling ENVIRONMENT")
        print(f"Run: conda activate mld_modeling")
        print(f"Then: python test_current_setup.py")
    
    print(f"=" * 50)

if __name__ == "__main__":
    main()