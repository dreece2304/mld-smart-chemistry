#!/usr/bin/env python3
"""
Simple molecule visualization for optimized structures
Works in headless mode (perfect for your setup)
"""

import os
import matplotlib.pyplot as plt
from ase.io import read
from ase.visualize.plot import plot_atoms
import numpy as np

def visualize_molecule(filename, title=None, save_png=True):
    """Simple 2D visualization of molecule"""
    
    if not os.path.exists(filename):
        print(f"❌ File not found: {filename}")
        return False
    
    # Read molecule
    atoms = read(filename)
    
    # Set up matplotlib for headless
    plt.switch_backend('Agg')
    
    # Create visualization
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot atoms
    plot_atoms(atoms, ax, radii=0.5, colors=None)
    
    # Customize plot
    title = title or f"Molecule: {os.path.basename(filename)}"
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel('X (Å)', fontsize=12)
    ax.set_ylabel('Y (Å)', fontsize=12)
    
    # Add atom count and energy info
    info_text = f"Atoms: {len(atoms)}"
    try:
        energy = atoms.get_potential_energy()
        info_text += f"\nEnergy: {energy:.3f} eV"
    except:
        pass
    
    ax.text(0.02, 0.98, info_text, transform=ax.transAxes, 
            verticalalignment='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Save plot
    if save_png:
        output_file = filename.replace('.xyz', '_plot.png')
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"📊 Visualization saved: {output_file}")
    
    plt.close()
    return True

def compare_molecules(file1, file2, title="Molecule Comparison"):
    """Compare two molecular structures side by side"""
    
    # Read molecules
    try:
        mol1 = read(file1)
        mol2 = read(file2)
    except Exception as e:
        print(f"❌ Error reading files: {e}")
        return False
    
    # Set up matplotlib
    plt.switch_backend('Agg')
    
    # Create side-by-side plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot molecule 1
    plot_atoms(mol1, ax1, radii=0.5)
    ax1.set_title(f"{os.path.basename(file1)}\n({len(mol1)} atoms)")
    ax1.set_xlabel('X (Å)')
    ax1.set_ylabel('Y (Å)')
    
    # Plot molecule 2  
    plot_atoms(mol2, ax2, radii=0.5)
    ax2.set_title(f"{os.path.basename(file2)}\n({len(mol2)} atoms)")
    ax2.set_xlabel('X (Å)')
    ax2.set_ylabel('Y (Å)')
    
    plt.suptitle(title, fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    # Save comparison
    output_file = "molecule_comparison.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"📊 Comparison saved: {output_file}")
    
    plt.close()
    return True

def visualize_optimization_trajectory(traj_file, molecule_name="Molecule"):
    """Visualize optimization trajectory"""
    
    if not os.path.exists(traj_file):
        print(f"❌ Trajectory file not found: {traj_file}")
        return False
    
    from ase.io import read
    
    try:
        # Read trajectory
        trajectory = read(traj_file, ':')
        print(f"📈 Trajectory has {len(trajectory)} steps")
        
        # Extract energies and forces
        energies = []
        max_forces = []
        
        for atoms in trajectory:
            try:
                energies.append(atoms.get_potential_energy())
                forces = atoms.get_forces()
                max_forces.append(np.max(np.sqrt(np.sum(forces**2, axis=1))))
            except:
                pass
        
        # Plot convergence
        plt.switch_backend('Agg')
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
        
        # Energy plot
        if energies:
            ax1.plot(energies, 'b-o', linewidth=2, markersize=4)
            ax1.set_title(f'{molecule_name} Energy Convergence')
            ax1.set_xlabel('Optimization Step')
            ax1.set_ylabel('Energy (eV)')
            ax1.grid(True, alpha=0.3)
        
        # Force plot
        if max_forces:
            ax2.plot(max_forces, 'r-s', linewidth=2, markersize=4)
            ax2.axhline(y=0.05, color='green', linestyle='--', label='Convergence (0.05 eV/Å)')
            ax2.set_title(f'{molecule_name} Force Convergence')
            ax2.set_xlabel('Optimization Step')
            ax2.set_ylabel('Max Force (eV/Å)')
            ax2.legend()
            ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Save plot
        output_file = f"{molecule_name}_convergence.png"
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"📊 Convergence plot saved: {output_file}")
        
        plt.close()
        return True
        
    except Exception as e:
        print(f"❌ Error visualizing trajectory: {e}")
        return False

def main():
    """Demo visualization functions"""
    print("🎨 Molecule Visualization Tools")
    print("=" * 40)
    
    # Check for optimized molecules
    molecules = ['tma_optimized.xyz', 'trimethylaluminum_optimized.xyz', 'h2o_optimized.xyz', 'water_optimized.xyz']
    
    for mol_file in molecules:
        if os.path.exists(mol_file):
            print(f"\n📊 Visualizing {mol_file}...")
            visualize_molecule(mol_file)
    
    # Check for trajectory files
    trajectories = ['tma_dft.traj', 'trimethylaluminum_dft.traj', 'h2o_dft.traj', 'water_dft.traj']
    
    for traj_file in trajectories:
        if os.path.exists(traj_file):
            mol_name = traj_file.replace('_dft.traj', '').replace('_', ' ').title()
            print(f"\n📈 Analyzing {traj_file}...")
            visualize_optimization_trajectory(traj_file, mol_name)

if __name__ == "__main__":
    main()