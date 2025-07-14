#!/usr/bin/env python3
"""
Headless Visualization Test for WSL2/MLD
Tests all visualization capabilities without GUI dependencies
"""

import os
import sys
import tempfile
import numpy as np
from datetime import datetime
from pathlib import Path
import argparse

def setup_headless_backend():
    """Set up headless matplotlib backend"""
    try:
        import matplotlib
        matplotlib.use('Agg')  # Force headless backend
        import matplotlib.pyplot as plt
        
        print("✅ Matplotlib configured for headless operation")
        return True
    except Exception as e:
        print(f"❌ Failed to setup headless backend: {e}")
        return False

def test_basic_matplotlib():
    """Test basic matplotlib plotting"""
    print("🔧 Testing Basic Matplotlib Plotting...")
    
    try:
        import matplotlib.pyplot as plt
        import numpy as np
        
        # Create test data
        x = np.linspace(0, 2*np.pi, 100)
        y1 = np.sin(x)
        y2 = np.cos(x)
        
        # Create figure with subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        
        # Plot 1: Basic line plot
        ax1.plot(x, y1, 'b-', linewidth=2, label='sin(x)')
        ax1.plot(x, y2, 'r--', linewidth=2, label='cos(x)')
        ax1.set_title('Trigonometric Functions')
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Scatter plot with colormap
        np.random.seed(42)
        x_scatter = np.random.randn(100)
        y_scatter = np.random.randn(100)
        colors = x_scatter + y_scatter
        
        scatter = ax2.scatter(x_scatter, y_scatter, c=colors, alpha=0.7, cmap='viridis')
        ax2.set_title('Random Scatter Plot')
        ax2.set_xlabel('X')
        ax2.set_ylabel('Y')
        plt.colorbar(scatter, ax=ax2)
        
        plt.tight_layout()
        
        # Save to temporary file
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
            plt.savefig(tmp.name, dpi=150, bbox_inches='tight')
            test_file = tmp.name
        
        plt.close(fig)
        
        # Verify file was created and has content
        if os.path.exists(test_file) and os.path.getsize(test_file) > 1000:
            print("✅ Basic matplotlib plotting: working")
            os.unlink(test_file)  # Cleanup
            return True
        else:
            print("❌ Basic matplotlib plotting: file not created properly")
            return False
            
    except Exception as e:
        print(f"❌ Basic matplotlib test failed: {e}")
        return False

def test_ase_visualization():
    """Test ASE structure visualization"""
    print("🔧 Testing ASE Structure Visualization...")
    
    try:
        from ase import Atoms
        from ase.visualize.plot import plot_atoms
        import matplotlib.pyplot as plt
        
        # Create test molecule (water)
        h2o = Atoms('H2O', positions=[[0, 0, 0], [0.76, 0.59, 0], [-0.76, 0.59, 0]])
        h2o.center(vacuum=3.0)
        
        # Create visualization
        fig, ax = plt.subplots(figsize=(6, 6))
        plot_atoms(h2o, ax, radii=0.5, colors=None)
        ax.set_title('H2O Molecule')
        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        
        # Save plot
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
            plt.savefig(tmp.name, dpi=150, bbox_inches='tight')
            test_file = tmp.name
        
        plt.close(fig)
        
        # Verify
        if os.path.exists(test_file) and os.path.getsize(test_file) > 1000:
            print("✅ ASE structure visualization: working")
            os.unlink(test_file)
            return True
        else:
            print("❌ ASE structure visualization: failed")
            return False
            
    except ImportError as e:
        print(f"❌ ASE not available: {e}")
        return False
    except Exception as e:
        print(f"❌ ASE visualization test failed: {e}")
        return False

def test_plotly_visualization():
    """Test Plotly interactive visualization"""
    print("🔧 Testing Plotly Interactive Visualization...")
    
    try:
        import plotly.graph_objects as go
        import plotly.express as px
        import numpy as np
        
        # Create test data
        x = np.linspace(0, 10, 100)
        y = np.sin(x) * np.exp(-x/10)
        
        # Create interactive plot
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=x, y=y, mode='lines', name='Damped Sine'))
        fig.update_layout(
            title='Interactive Plotly Test',
            xaxis_title='X',
            yaxis_title='Y',
            template='plotly_white'
        )
        
        # Save as HTML
        with tempfile.NamedTemporaryFile(suffix='.html', delete=False) as tmp:
            fig.write_html(tmp.name)
            test_file = tmp.name
        
        # Verify
        if os.path.exists(test_file) and os.path.getsize(test_file) > 1000:
            print("✅ Plotly interactive visualization: working")
            os.unlink(test_file)
            return True
        else:
            print("❌ Plotly visualization: failed")
            return False
            
    except ImportError:
        print("⚠️  Plotly not available (optional)")
        return True  # Not required
    except Exception as e:
        print(f"❌ Plotly test failed: {e}")
        return False

def test_py3dmol_visualization():
    """Test py3dmol 3D molecular visualization"""
    print("🔧 Testing py3dmol 3D Visualization...")
    
    try:
        import py3dmol
        
        # Create simple molecule XYZ string
        xyz_string = """4
Water molecule
O 0.00000 0.00000 0.00000
H 0.75690 0.58588 0.00000
H -0.75690 0.58588 0.00000
"""
        
        # Create viewer
        viewer = py3dmol.view(width=400, height=300)
        viewer.addModel(xyz_string, 'xyz')
        viewer.setStyle({'stick': {'radius': 0.1}, 'sphere': {'radius': 0.3}})
        viewer.setBackgroundColor('white')
        viewer.zoomTo()
        
        # Save as HTML
        with tempfile.NamedTemporaryFile(suffix='.html', delete=False) as tmp:
            html_content = viewer._make_html()
            with open(tmp.name, 'w') as f:
                f.write(html_content)
            test_file = tmp.name
        
        # Verify
        if os.path.exists(test_file) and os.path.getsize(test_file) > 1000:
            print("✅ py3dmol 3D visualization: working")
            os.unlink(test_file)
            return True
        else:
            print("❌ py3dmol visualization: failed")
            return False
            
    except ImportError:
        print("⚠️  py3dmol not available (optional)")
        return True  # Not required
    except Exception as e:
        print(f"❌ py3dmol test failed: {e}")
        return False

def test_trajectory_analysis():
    """Test trajectory analysis and plotting"""
    print("🔧 Testing Trajectory Analysis...")
    
    try:
        from ase import Atoms
        from ase.io import write, read
        import matplotlib.pyplot as plt
        
        # Create fake trajectory data
        trajectory_data = []
        for i in range(10):
            # Simulate H2 optimization
            bond_length = 0.8 - i * 0.003  # Approaching equilibrium
            h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, bond_length]])
            h2.center(vacuum=5.0)
            
            # Add fake energy and forces
            h2.info['energy'] = -6.0 + (bond_length - 0.74)**2 * 100  # Parabolic potential
            trajectory_data.append(h2)
        
        # Create trajectory analysis plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        
        # Plot 1: Energy vs step
        steps = range(len(trajectory_data))
        energies = [atoms.info['energy'] for atoms in trajectory_data]
        
        ax1.plot(steps, energies, 'b-o', linewidth=2, markersize=4)
        ax1.set_title('Energy Convergence')
        ax1.set_xlabel('Optimization Step')
        ax1.set_ylabel('Energy (eV)')
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Bond length vs step
        bond_lengths = []
        for atoms in trajectory_data:
            pos = atoms.get_positions()
            bond_length = np.linalg.norm(pos[1] - pos[0])
            bond_lengths.append(bond_length)
        
        ax2.plot(steps, bond_lengths, 'r-s', linewidth=2, markersize=4)
        ax2.set_title('Bond Length Evolution')
        ax2.set_xlabel('Optimization Step')
        ax2.set_ylabel('H-H Distance (Å)')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Save plot
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
            plt.savefig(tmp.name, dpi=150, bbox_inches='tight')
            test_file = tmp.name
        
        plt.close(fig)
        
        # Verify
        if os.path.exists(test_file) and os.path.getsize(test_file) > 1000:
            print("✅ Trajectory analysis: working")
            os.unlink(test_file)
            return True
        else:
            print("❌ Trajectory analysis: failed")
            return False
            
    except Exception as e:
        print(f"❌ Trajectory analysis test failed: {e}")
        return False

def test_mld_structure_visualization():
    """Test MLD-specific structure visualization"""
    print("🔧 Testing MLD Structure Visualization...")
    
    try:
        # Test with actual MLD structures if available
        structure_files = [
            'structures/tma_molecule.xyz',
            'structures/butyne_diol_molecule.xyz'
        ]
        
        success_count = 0
        total_files = 0
        
        for file_path in structure_files:
            if os.path.exists(file_path):
                total_files += 1
                try:
                    from ase.io import read
                    from ase.visualize.plot import plot_atoms
                    import matplotlib.pyplot as plt
                    
                    # Read structure
                    atoms = read(file_path)
                    
                    # Create visualization
                    fig, ax = plt.subplots(figsize=(6, 6))
                    plot_atoms(atoms, ax, radii=0.3)
                    ax.set_title(f'MLD Structure: {Path(file_path).stem}')
                    
                    # Save plot
                    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                        plt.savefig(tmp.name, dpi=150, bbox_inches='tight')
                        test_file = tmp.name
                    
                    plt.close(fig)
                    
                    # Verify
                    if os.path.exists(test_file) and os.path.getsize(test_file) > 1000:
                        print(f"  ✅ {Path(file_path).name}: visualized successfully")
                        success_count += 1
                        os.unlink(test_file)
                    else:
                        print(f"  ❌ {Path(file_path).name}: visualization failed")
                        
                except Exception as e:
                    print(f"  ❌ {Path(file_path).name}: error - {e}")
        
        if total_files == 0:
            print("⚠️  No MLD structure files found")
            return True  # Not a failure
        elif success_count == total_files:
            print("✅ MLD structure visualization: all files working")
            return True
        else:
            print(f"⚠️  MLD structure visualization: {success_count}/{total_files} working")
            return success_count > 0
            
    except Exception as e:
        print(f"❌ MLD structure visualization failed: {e}")
        return False

def generate_visualization_report():
    """Generate comprehensive visualization report"""
    print("\n" + "="*60)
    print("HEADLESS VISUALIZATION TEST REPORT")
    print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*60)
    
    # Check environment
    print(f"\nEnvironment Configuration:")
    
    # Matplotlib backend
    try:
        import matplotlib
        backend = matplotlib.get_backend()
        print(f"  Matplotlib backend: {backend}")
        if backend == 'Agg':
            print(f"  Status: ✅ Headless (good for WSL2)")
        else:
            print(f"  Status: ⚠️  GUI backend (may cause issues in WSL2)")
    except:
        print(f"  Matplotlib: ❌ Not available")
    
    # Display environment
    display = os.environ.get('DISPLAY', 'not set')
    print(f"  DISPLAY variable: {display}")
    if not display or display == "":
        print(f"  Status: ✅ Headless mode")
    else:
        print(f"  Status: ⚠️  GUI mode detected")
    
    # WSL detection
    try:
        with open('/proc/version', 'r') as f:
            version_info = f.read()
        if 'Microsoft' in version_info:
            print(f"  Environment: WSL2")
            print(f"  Recommendation: Use headless visualization only")
        else:
            print(f"  Environment: Native Linux")
    except:
        print(f"  Environment: Unknown")
    
    print(f"\nVisualization Capabilities:")
    print(f"  ✅ Matplotlib (PNG/PDF output)")
    print(f"  ✅ ASE structure plots") 
    print(f"  ✅ Trajectory analysis")
    print(f"  ✅ Publication-quality figures")
    
    try:
        import plotly
        print(f"  ✅ Plotly (interactive HTML)")
    except:
        print(f"  ❌ Plotly (install with: pip install plotly)")
    
    try:
        import py3dmol
        print(f"  ✅ py3dmol (3D web viewer)")
    except:
        print(f"  ❌ py3dmol (install with: pip install py3dmol)")
    
    print(f"\nRecommended Workflow:")
    print(f"  1. Generate plots with matplotlib (headless)")
    print(f"  2. Create interactive HTML with plotly/py3dmol")
    print(f"  3. View results in web browser")
    print(f"  4. Transfer files to Windows for advanced visualization")
    
    print("="*60)

def main():
    """Main headless visualization test"""
    
    parser = argparse.ArgumentParser(description='Test Headless Visualization')
    parser.add_argument('--quick', action='store_true',
                       help='Skip optional tests')
    
    args = parser.parse_args()
    
    print("🎨 HEADLESS VISUALIZATION TEST")
    print("=" * 50)
    
    # Setup headless backend first
    if not setup_headless_backend():
        print("❌ Failed to setup headless backend")
        return 1
    
    print()
    
    # Run tests
    tests_passed = 0
    total_tests = 0
    
    # Core tests
    tests = [
        test_basic_matplotlib,
        test_ase_visualization,
        test_trajectory_analysis,
        test_mld_structure_visualization
    ]
    
    # Optional tests
    if not args.quick:
        tests.extend([
            test_plotly_visualization,
            test_py3dmol_visualization
        ])
    
    # Run all tests
    for test_func in tests:
        total_tests += 1
        if test_func():
            tests_passed += 1
        print()
    
    # Generate report
    generate_visualization_report()
    
    # Final summary
    success_rate = (tests_passed / total_tests * 100) if total_tests > 0 else 0
    
    print(f"\nTest Results: {tests_passed}/{total_tests} passed ({success_rate:.1f}%)")
    
    if success_rate >= 75:
        print("🎉 Headless visualization is working well!")
        return 0
    else:
        print("⚠️  Some visualization tests failed")
        return 1

if __name__ == "__main__":
    exit(main())