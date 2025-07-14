#!/usr/bin/env python3
"""
Comprehensive MLD Visualization Tools
Provides multiple visualization options for molecular structures and calculations
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse

def setup_visualization_environment():
    """Set up the best available visualization tools"""
    
    available_tools = {}
    
    # Check ASE GUI
    try:
        from ase.gui.gui import GUI
        available_tools['ase_gui'] = True
        print("✅ ASE GUI available")
    except ImportError:
        available_tools['ase_gui'] = False
        print("❌ ASE GUI not available")
    
    # Check OVITO
    try:
        import ovito
        available_tools['ovito'] = True
        print(f"✅ OVITO {ovito.__version__} available")
    except ImportError:
        available_tools['ovito'] = False
        print("❌ OVITO not available")
    
    # Check nglview (for Jupyter)
    try:
        import nglview
        available_tools['nglview'] = True
        print("✅ NGLView available (Jupyter)")
    except ImportError:
        available_tools['nglview'] = False
        print("❌ NGLView not available")
    
    # Check py3dmol
    try:
        import py3dmol
        available_tools['py3dmol'] = True
        print("✅ Py3DMol available")
    except ImportError:
        available_tools['py3dmol'] = False
        print("❌ Py3DMol not available")
    
    # Check plotly
    try:
        import plotly
        available_tools['plotly'] = True
        print("✅ Plotly available")
    except ImportError:
        available_tools['plotly'] = False
        print("❌ Plotly not available")
    
    return available_tools

def visualize_with_ase(structure_file, output_file=None):
    """Visualize structure using ASE matplotlib backend"""
    
    from ase.io import read
    from ase.visualize import view
    import matplotlib.pyplot as plt
    
    print(f"📊 Visualizing {structure_file} with ASE...")
    
    # Read structure
    atoms = read(structure_file)
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    # Plot 1: 3D structure
    from ase.visualize.plot import plot_atoms
    plot_atoms(atoms, ax1, radii=0.5, colors=None)
    ax1.set_title(f"Structure: {Path(structure_file).stem}")
    ax1.set_xlabel("X (Å)")
    ax1.set_ylabel("Y (Å)")
    
    # Plot 2: Bond lengths analysis
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    
    # Calculate bond lengths
    from scipy.spatial.distance import pdist, squareform
    distances = squareform(pdist(positions))
    
    # Find reasonable bonds (< 3 Å)
    bond_lengths = []
    bond_types = []
    
    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            dist = distances[i, j]
            if dist < 3.0:  # Reasonable bond cutoff
                bond_lengths.append(dist)
                bond_types.append(f"{symbols[i]}-{symbols[j]}")
    
    # Plot bond length distribution
    ax2.hist(bond_lengths, bins=20, alpha=0.7, edgecolor='black')
    ax2.set_title("Bond Length Distribution")
    ax2.set_xlabel("Distance (Å)")
    ax2.set_ylabel("Count")
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"💾 Saved plot: {output_file}")
    else:
        plt.show()
    
    return fig

def visualize_with_ovito(structure_file, output_file=None):
    """Visualize structure using OVITO"""
    
    try:
        from ovito.io import import_file
        from ovito.vis import Viewport
        import ovito
        
        print(f"🎨 Visualizing {structure_file} with OVITO...")
        
        # Import structure
        pipeline = import_file(structure_file)
        
        # Set up visualization
        pipeline.add_to_scene()
        
        # Render image
        if output_file:
            vp = Viewport()
            vp.type = Viewport.Type.Perspective
            vp.camera_pos = (10, 10, 10)
            vp.camera_dir = (-1, -1, -1)
            vp.render_image(size=(800, 600), filename=output_file)
            print(f"💾 Saved OVITO render: {output_file}")
        
        print("✅ OVITO visualization complete")
        return True
        
    except ImportError:
        print("❌ OVITO not available")
        return False
    except Exception as e:
        print(f"❌ OVITO error: {e}")
        return False

def create_interactive_plot(structure_file):
    """Create interactive plot with py3dmol"""
    
    try:
        import py3dmol
        from ase.io import read
        
        print(f"🌐 Creating interactive view of {structure_file}...")
        
        atoms = read(structure_file)
        
        # Convert to XYZ string
        xyz_string = f"{len(atoms)}\n"
        xyz_string += f"Structure from {structure_file}\n"
        
        for atom in atoms:
            pos = atom.position
            xyz_string += f"{atom.symbol} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n"
        
        # Create viewer
        viewer = py3dmol.view(width=800, height=600)
        viewer.addModel(xyz_string, 'xyz')
        viewer.setStyle({'stick': {'radius': 0.1}, 'sphere': {'radius': 0.3}})
        viewer.setBackgroundColor('white')
        viewer.zoomTo()
        
        # Save as HTML
        html_file = structure_file.replace('.xyz', '_interactive.html')
        viewer.write_html(html_file)
        print(f"💾 Saved interactive plot: {html_file}")
        
        return viewer
        
    except ImportError:
        print("❌ py3dmol not available")
        return None
    except Exception as e:
        print(f"❌ Interactive plot error: {e}")
        return None

def analyze_trajectory(traj_file, output_dir=None):
    """Analyze optimization trajectory"""
    
    from ase.io import read
    import matplotlib.pyplot as plt
    
    print(f"📈 Analyzing trajectory: {traj_file}")
    
    if output_dir is None:
        output_dir = Path(traj_file).parent
    
    # Read trajectory
    traj = read(traj_file, ':')
    
    # Extract data
    steps = range(len(traj))
    energies = [atoms.get_potential_energy() for atoms in traj if hasattr(atoms, 'get_potential_energy')]
    
    # Calculate forces if available
    forces = []
    max_forces = []
    
    for atoms in traj:
        try:
            force_array = atoms.get_forces()
            forces.append(force_array)
            max_forces.append(np.max(np.linalg.norm(force_array, axis=1)))
        except:
            max_forces.append(None)
    
    # Create analysis plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'Trajectory Analysis: {Path(traj_file).stem}')
    
    # Energy convergence
    if energies:
        axes[0,0].plot(steps[:len(energies)], energies, 'b-o', markersize=3)
        axes[0,0].set_title('Energy Convergence')
        axes[0,0].set_xlabel('Step')
        axes[0,0].set_ylabel('Energy (eV)')
        axes[0,0].grid(True, alpha=0.3)
    
    # Force convergence
    valid_forces = [f for f in max_forces if f is not None]
    if valid_forces:
        axes[0,1].semilogy(range(len(valid_forces)), valid_forces, 'r-o', markersize=3)
        axes[0,1].set_title('Force Convergence')
        axes[0,1].set_xlabel('Step')
        axes[0,1].set_ylabel('Max Force (eV/Å)')
        axes[0,1].grid(True, alpha=0.3)
        
        # Add convergence target
        axes[0,1].axhline(y=0.05, color='gray', linestyle='--', label='Target (0.05 eV/Å)')
        axes[0,1].legend()
    
    # Structural parameters
    if len(traj) > 1:
        # Track a bond length (if possible)
        try:
            positions = [atoms.get_positions() for atoms in traj]
            if len(positions[0]) >= 2:
                bond_lengths = []
                for pos in positions:
                    # Distance between first two atoms
                    dist = np.linalg.norm(pos[0] - pos[1])
                    bond_lengths.append(dist)
                
                axes[1,0].plot(steps[:len(bond_lengths)], bond_lengths, 'g-o', markersize=3)
                axes[1,0].set_title('Bond Length Evolution')
                axes[1,0].set_xlabel('Step')
                axes[1,0].set_ylabel('Distance (Å)')
                axes[1,0].grid(True, alpha=0.3)
        except:
            axes[1,0].text(0.5, 0.5, 'Bond analysis\nnot available', 
                          ha='center', va='center', transform=axes[1,0].transAxes)
    
    # Final structure
    if traj:
        final_atoms = traj[-1]
        try:
            from ase.visualize.plot import plot_atoms
            plot_atoms(final_atoms, axes[1,1], radii=0.3)
            axes[1,1].set_title('Final Structure')
        except:
            axes[1,1].text(0.5, 0.5, 'Structure plot\nnot available', 
                          ha='center', va='center', transform=axes[1,1].transAxes)
    
    plt.tight_layout()
    
    # Save analysis
    output_file = Path(output_dir) / f"{Path(traj_file).stem}_analysis.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"💾 Saved trajectory analysis: {output_file}")
    
    plt.show()
    
    return fig

def create_mld_visualization_report(structure_dir='.', output_file='mld_visualization_report.html'):
    """Create comprehensive HTML report with all visualizations"""
    
    structure_dir = Path(structure_dir)
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>MLD Visualization Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            .structure {{ border: 1px solid #ddd; margin: 20px 0; padding: 15px; }}
            .structure h3 {{ color: #2c5aa0; }}
            img {{ max-width: 100%; height: auto; }}
            .grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
        </style>
    </head>
    <body>
        <h1>🧪 MLD Simulation Visualization Report</h1>
        <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        
        <h2>📁 Structure Files</h2>
    """
    
    # Find all structure files
    xyz_files = list(structure_dir.glob('**/*.xyz'))
    traj_files = list(structure_dir.glob('**/*.traj'))
    
    # Visualize structures
    for xyz_file in xyz_files:
        html_content += f"""
        <div class="structure">
            <h3>🧬 {xyz_file.name}</h3>
            <p>File: {xyz_file}</p>
        """
        
        # Create visualizations
        plot_file = xyz_file.with_suffix('.png')
        try:
            visualize_with_ase(str(xyz_file), str(plot_file))
            html_content += f'<img src="{plot_file.name}" alt="Structure plot">'
        except Exception as e:
            html_content += f'<p>❌ Visualization failed: {e}</p>'
        
        html_content += "</div>"
    
    # Analyze trajectories
    for traj_file in traj_files:
        html_content += f"""
        <div class="structure">
            <h3>📈 {traj_file.name}</h3>
            <p>Trajectory: {traj_file}</p>
        """
        
        try:
            analyze_trajectory(str(traj_file), str(traj_file.parent))
            analysis_file = traj_file.with_name(f"{traj_file.stem}_analysis.png")
            if analysis_file.exists():
                html_content += f'<img src="{analysis_file.name}" alt="Trajectory analysis">'
        except Exception as e:
            html_content += f'<p>❌ Analysis failed: {e}</p>'
        
        html_content += "</div>"
    
    html_content += """
        </body>
    </html>
    """
    
    # Save report
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"📄 Saved visualization report: {output_file}")
    return output_file

def main():
    """Main visualization tool"""
    
    parser = argparse.ArgumentParser(description='MLD Visualization Tools')
    parser.add_argument('--structure', help='Structure file to visualize')
    parser.add_argument('--trajectory', help='Trajectory file to analyze')
    parser.add_argument('--report', action='store_true', help='Generate full report')
    parser.add_argument('--dir', default='.', help='Directory to scan')
    parser.add_argument('--output', help='Output file')
    parser.add_argument('--check', action='store_true', help='Check available tools')
    
    args = parser.parse_args()
    
    if args.check:
        print("🔍 Checking available visualization tools...")
        setup_visualization_environment()
        return
    
    if args.structure:
        print(f"🧬 Visualizing structure: {args.structure}")
        
        # Try multiple visualization methods
        visualize_with_ase(args.structure, args.output)
        visualize_with_ovito(args.structure, args.output)
        create_interactive_plot(args.structure)
        
    elif args.trajectory:
        print(f"📈 Analyzing trajectory: {args.trajectory}")
        analyze_trajectory(args.trajectory)
        
    elif args.report:
        print(f"📄 Generating visualization report for: {args.dir}")
        output_file = args.output or 'mld_visualization_report.html'
        create_mld_visualization_report(args.dir, output_file)
        
    else:
        print("🎨 MLD Visualization Tools")
        print("Usage examples:")
        print("  python visualization_tools.py --check")
        print("  python visualization_tools.py --structure tma_molecule.xyz")
        print("  python visualization_tools.py --trajectory optimization.traj")
        print("  python visualization_tools.py --report --dir results/")

if __name__ == "__main__":
    from datetime import datetime
    main()