#!/usr/bin/env python3
"""
Visualize MLD simulation results
Creates images and animations of the MLD growth process
"""

import os
import glob
import matplotlib.pyplot as plt
from ase.io import read
from ase.visualize.plot import plot_atoms
import numpy as np

def visualize_mld_cycle_progression():
    """Visualize the progression through MLD cycles"""
    
    print("🎬 Creating MLD cycle progression visualization...")
    
    # Find all MLD cycle files
    cycle_files = sorted(glob.glob('mld_cycle*_after_*.xyz'))
    
    if not cycle_files:
        print("❌ No MLD cycle files found!")
        print("Run step5_realistic_mld.py first")
        return
    
    print(f"Found {len(cycle_files)} MLD structures")
    
    # Set up matplotlib for headless operation
    plt.switch_backend('Agg')
    
    # Create grid layout for all cycles
    n_files = len(cycle_files)
    cols = min(4, n_files)
    rows = (n_files + cols - 1) // cols
    
    fig, axes = plt.subplots(rows, cols, figsize=(4*cols, 3*rows))
    if n_files == 1:
        axes = [axes]
    elif rows == 1:
        axes = axes
    else:
        axes = axes.flatten()
    
    for i, cycle_file in enumerate(cycle_files):
        if i >= len(axes):
            break
            
        try:
            atoms = read(cycle_file)
            
            # Plot structure
            ax = axes[i]
            plot_atoms(atoms, ax, radii=0.3)
            
            # Extract cycle info from filename
            filename = os.path.basename(cycle_file)
            if 'after_tma' in filename:
                phase = 'After TMA'
                color = 'lightblue'
            elif 'after_h2o' in filename:
                phase = 'After H2O'
                color = 'lightgreen'
            else:
                phase = 'Unknown'
                color = 'lightgray'
            
            # Get cycle number
            cycle_num = filename.split('cycle')[1].split('_')[0]
            
            ax.set_title(f'Cycle {cycle_num}\n{phase}\n({len(atoms)} atoms)', 
                        fontsize=10, bbox=dict(boxstyle='round', facecolor=color, alpha=0.7))
            ax.set_xlabel('X (Å)', fontsize=8)
            ax.set_ylabel('Y (Å)', fontsize=8)
            
        except Exception as e:
            print(f"Error plotting {cycle_file}: {e}")
    
    # Hide unused subplots
    for i in range(len(cycle_files), len(axes)):
        axes[i].set_visible(False)
    
    plt.suptitle('MLD Cycle Progression', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    # Save plot
    plt.savefig('mld_cycle_progression.png', dpi=200, bbox_inches='tight')
    print(f"📊 Saved: mld_cycle_progression.png")
    plt.close()

def visualize_growth_analysis():
    """Analyze and visualize MLD growth statistics"""
    
    print("📈 Creating MLD growth analysis...")
    
    cycle_files = sorted(glob.glob('mld_cycle*_after_h2o.xyz'))
    
    if len(cycle_files) < 2:
        print("❌ Need at least 2 completed cycles for growth analysis")
        return
    
    # Analyze growth
    cycles = []
    atom_counts = []
    
    for cycle_file in cycle_files:
        atoms = read(cycle_file)
        cycle_num = int(cycle_file.split('cycle')[1].split('_')[0])
        cycles.append(cycle_num)
        atom_counts.append(len(atoms))
    
    # Calculate growth per cycle
    if len(atom_counts) > 1:
        growth_per_cycle = [atom_counts[i] - atom_counts[i-1] for i in range(1, len(atom_counts))]
    else:
        growth_per_cycle = []
    
    plt.switch_backend('Agg')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot 1: Total atoms vs cycle
    ax1.plot(cycles, atom_counts, 'bo-', linewidth=2, markersize=8)
    ax1.set_xlabel('MLD Cycle')
    ax1.set_ylabel('Total Atoms')
    ax1.set_title('MLD Growth: Total Atoms')
    ax1.grid(True, alpha=0.3)
    
    # Add annotations
    for i, (cycle, count) in enumerate(zip(cycles, atom_counts)):
        ax1.annotate(f'{count}', (cycle, count), textcoords="offset points", 
                    xytext=(0,10), ha='center', fontsize=9)
    
    # Plot 2: Growth per cycle
    if growth_per_cycle:
        cycle_labels = [f'Cycle {cycles[i]} → {cycles[i+1]}' for i in range(len(growth_per_cycle))]
        ax2.bar(range(len(growth_per_cycle)), growth_per_cycle, 
               color=['lightblue' if i%2==0 else 'lightcoral' for i in range(len(growth_per_cycle))])
        ax2.set_xlabel('Cycle Transition')
        ax2.set_ylabel('Atoms Added')
        ax2.set_title('MLD Growth: Atoms Added Per Cycle')
        ax2.set_xticks(range(len(growth_per_cycle)))
        ax2.set_xticklabels(cycle_labels, rotation=45)
        
        # Add value labels on bars
        for i, v in enumerate(growth_per_cycle):
            ax2.text(i, v + 0.1, str(v), ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.savefig('mld_growth_analysis.png', dpi=200, bbox_inches='tight')
    print(f"📊 Saved: mld_growth_analysis.png")
    plt.close()

def compare_before_after():
    """Compare initial surface with final MLD structure"""
    
    print("🔄 Creating before/after comparison...")
    
    # Find initial and final structures
    initial_files = ['uv_ozone_si_surface.xyz', 'surface_si100_hydroxylated.xyz', 'simple_surface.xyz']
    final_file = 'mld_final_structure.xyz'
    
    initial_surface = None
    for init_file in initial_files:
        if os.path.exists(init_file):
            initial_surface = read(init_file)
            initial_name = init_file
            break
    
    if not initial_surface:
        print("❌ No initial surface found!")
        return
    
    if not os.path.exists(final_file):
        print("❌ No final MLD structure found!")
        return
    
    final_surface = read(final_file)
    
    plt.switch_backend('Agg')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot initial surface
    plot_atoms(initial_surface, ax1, radii=0.4)
    ax1.set_title(f'Initial Surface\n{os.path.basename(initial_name)}\n({len(initial_surface)} atoms)', 
                 fontsize=12)
    ax1.set_xlabel('X (Å)')
    ax1.set_ylabel('Y (Å)')
    
    # Plot final surface
    plot_atoms(final_surface, ax2, radii=0.4)
    ax2.set_title(f'After MLD\nmld_final_structure.xyz\n({len(final_surface)} atoms)', 
                 fontsize=12)
    ax2.set_xlabel('X (Å)')
    ax2.set_ylabel('Y (Å)')
    
    # Calculate growth
    atoms_added = len(final_surface) - len(initial_surface)
    
    plt.suptitle(f'MLD Before/After Comparison\nGrowth: +{atoms_added} atoms', 
                fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    plt.savefig('mld_before_after.png', dpi=200, bbox_inches='tight')
    print(f"📊 Saved: mld_before_after.png")
    plt.close()

def analyze_surface_composition():
    """Analyze the composition of MLD structures"""
    
    print("🧪 Analyzing surface composition...")
    
    final_file = 'mld_final_structure.xyz'
    if not os.path.exists(final_file):
        print("❌ No final structure found!")
        return
    
    atoms = read(final_file)
    symbols = atoms.get_chemical_symbols()
    
    # Count each element
    composition = {}
    for symbol in set(symbols):
        composition[symbol] = symbols.count(symbol)
    
    print(f"\n📊 Final structure composition:")
    total_atoms = len(atoms)
    for element, count in sorted(composition.items()):
        percentage = (count / total_atoms) * 100
        print(f"   {element}: {count} atoms ({percentage:.1f}%)")
    
    # Create composition pie chart
    plt.switch_backend('Agg')
    fig, ax = plt.subplots(figsize=(8, 8))
    
    elements = list(composition.keys())
    counts = list(composition.values())
    
    # Colors for common elements
    element_colors = {
        'Si': 'gray', 'O': 'red', 'H': 'white', 'Al': 'blue', 
        'C': 'black', 'N': 'lightblue'
    }
    colors = [element_colors.get(el, 'lightgray') for el in elements]
    
    wedges, texts, autotexts = ax.pie(counts, labels=elements, autopct='%1.1f%%', 
                                     colors=colors, startangle=90)
    
    ax.set_title('MLD Final Structure Composition', fontsize=16, fontweight='bold')
    
    plt.savefig('mld_composition.png', dpi=200, bbox_inches='tight')
    print(f"📊 Saved: mld_composition.png")
    plt.close()

def main():
    """Run all MLD visualization analyses"""
    
    print("🎨 MLD Results Visualization")
    print("=" * 50)
    
    # Check for MLD output files
    mld_files = glob.glob('mld_*.xyz')
    if not mld_files:
        print("❌ No MLD output files found!")
        print("Run step5_realistic_mld.py first")
        return 1
    
    print(f"Found {len(mld_files)} MLD output files")
    
    # Run visualizations
    visualize_mld_cycle_progression()
    visualize_growth_analysis()
    compare_before_after()
    analyze_surface_composition()
    
    print(f"\n✅ MLD visualization complete!")
    print(f"📁 Generated images:")
    print(f"   - mld_cycle_progression.png")
    print(f"   - mld_growth_analysis.png") 
    print(f"   - mld_before_after.png")
    print(f"   - mld_composition.png")
    
    print(f"\n💡 To view images:")
    print(f"   - Copy to local machine: scp user@server:~/path/*.png ./")
    print(f"   - Or use VSCode remote file browser")
    
    return 0

if __name__ == "__main__":
    exit(main())