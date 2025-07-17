#!/usr/bin/env python3
"""
Surface-Focused MLD Visualizer
Shows only the active surface layers and deposited molecules
Perfect for analyzing MLD chemistry without bulk Si clutter
"""

import os
import glob
import numpy as np
from ase.io import read
from collections import Counter
import sys

class SurfaceVisualizer:
    """Focus on surface layers and deposited chemistry"""
    
    def __init__(self, surface_depth=8.0):
        """
        Initialize surface visualizer
        
        Args:
            surface_depth: Depth in Å to consider as "surface" (default: 8 Å)
        """
        self.surface_depth = surface_depth
        
        # Terminal colors
        self.colors = {
            'red': '\033[91m', 'green': '\033[92m', 'yellow': '\033[93m',
            'blue': '\033[94m', 'magenta': '\033[95m', 'cyan': '\033[96m',
            'white': '\033[97m', 'bold': '\033[1m', 'end': '\033[0m'
        }
        
        # Surface-focused element styling
        self.surface_elements = {
            # Substrate (show but de-emphasize)
            'Si': {'symbol': '●', 'color': 'white', 'important': False},
            'O': {'symbol': '○', 'color': 'red', 'important': True},   # Surface oxide
            'H': {'symbol': '·', 'color': 'cyan', 'important': True},  # Surface OH
            
            # Deposited species (emphasize these!)
            'Al': {'symbol': '◆', 'color': 'blue', 'important': True},    # TMA aluminum
            'C': {'symbol': '◯', 'color': 'yellow', 'important': True},   # Methyl groups
            'N': {'symbol': '◉', 'color': 'magenta', 'important': True}   # If any amines
        }
    
    def colorize(self, text, color):
        """Add color to text"""
        return f"{self.colors[color]}{text}{self.colors['end']}"
    
    def extract_surface_region(self, atoms):
        """
        Extract only the surface region for visualization
        
        Args:
            atoms: Full structure
            
        Returns:
            (surface_atoms, surface_info)
        """
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        
        if len(positions) == 0:
            return atoms, {}
        
        # Find top of structure
        z_coords = positions[:, 2]
        z_max = z_coords.max()
        z_surface = z_max - self.surface_depth
        
        # Extract surface atoms
        surface_indices = []
        for i, pos in enumerate(positions):
            if pos[2] >= z_surface:
                surface_indices.append(i)
        
        surface_positions = positions[surface_indices]
        surface_symbols = [symbols[i] for i in surface_indices]
        
        # Create surface atoms object
        from ase import Atoms
        surface_atoms = Atoms(
            symbols=surface_symbols,
            positions=surface_positions
        )
        
        # Gather surface info
        surface_info = {
            'total_atoms': len(atoms),
            'surface_atoms': len(surface_atoms),
            'z_range': (z_surface, z_max),
            'surface_depth': self.surface_depth,
            'composition': Counter(surface_symbols)
        }
        
        return surface_atoms, surface_info
    
    def show_surface_layers(self, atoms, title="Surface Structure"):
        """
        Show surface structure in layers
        """
        surface_atoms, info = self.extract_surface_region(atoms)
        
        print(f"\n{self.colorize(title, 'bold')}")
        print(f"Surface region: {info['surface_atoms']}/{info['total_atoms']} atoms")
        print(f"Depth: {info['surface_depth']:.1f} Å (Z: {info['z_range'][0]:.1f} to {info['z_range'][1]:.1f})")
        
        if len(surface_atoms) == 0:
            print("No surface atoms found")
            return
        
        positions = surface_atoms.get_positions()
        symbols = surface_atoms.get_chemical_symbols()
        
        # Organize by Z layers (bins of 1 Å)
        z_coords = positions[:, 2]
        z_min, z_max = z_coords.min(), z_coords.max()
        
        layer_thickness = 1.0  # Å
        n_layers = max(1, int((z_max - z_min) / layer_thickness) + 1)
        
        print(f"\n{self.colorize('Surface Layer Analysis:', 'bold')}")
        print("─" * 60)
        
        for layer in range(n_layers):
            z_bottom = z_min + layer * layer_thickness
            z_top = z_bottom + layer_thickness
            
            # Find atoms in this layer
            layer_atoms = []
            for i, (pos, sym) in enumerate(zip(positions, symbols)):
                if z_bottom <= pos[2] < z_top:
                    layer_atoms.append((pos, sym))
            
            if not layer_atoms:
                continue
            
            # Analyze layer composition
            layer_symbols = [atom[1] for atom in layer_atoms]
            layer_comp = Counter(layer_symbols)
            
            print(f"\nLayer {layer+1} (Z: {z_bottom:.1f}-{z_top:.1f} Å) - {len(layer_atoms)} atoms")
            
            # Show composition with emphasis on important elements
            comp_str = ""
            for element, count in layer_comp.most_common():
                el_info = self.surface_elements.get(element, {'symbol': '?', 'color': 'white', 'important': False})
                symbol = el_info['symbol']
                color = el_info['color']
                
                if el_info['important']:
                    comp_str += f" {self.colorize(f'{symbol}×{count}', color)} {element}"
                else:
                    comp_str += f" {symbol}×{count} {element}"
            
            print(f"  Composition:{comp_str}")
            
            # Show 2D projection of this layer
            if len(layer_atoms) <= 50:  # Only for manageable sizes
                self.show_layer_2d_projection(layer_atoms, f"Layer {layer+1}")
    
    def show_layer_2d_projection(self, layer_atoms, layer_title):
        """Show 2D projection of a single layer"""
        if not layer_atoms:
            return
        
        positions = np.array([atom[0] for atom in layer_atoms])
        symbols = [atom[1] for atom in layer_atoms]
        
        # Use X-Y projection
        x_coords = positions[:, 0]
        y_coords = positions[:, 1]
        
        # Scale to terminal size
        width, height = 40, 15
        
        if len(x_coords) == 0:
            return
        
        x_min, x_max = x_coords.min(), x_coords.max()
        y_min, y_max = y_coords.min(), y_coords.max()
        
        x_range = x_max - x_min if x_max > x_min else 1
        y_range = y_max - y_min if y_max > y_min else 1
        
        # Create grid
        grid = [[' ' for _ in range(width)] for _ in range(height)]
        
        # Place atoms
        for i, (pos, symbol) in enumerate(layer_atoms):
            x, y = pos[0], pos[1]
            
            # Scale to grid
            grid_x = int(((x - x_min) / x_range) * (width - 1))
            grid_y = int(((y - y_min) / y_range) * (height - 1))
            
            if 0 <= grid_x < width and 0 <= grid_y < height:
                el_info = self.surface_elements.get(symbol, {'symbol': '?', 'color': 'white', 'important': False})
                atom_symbol = el_info['symbol']
                color = el_info['color']
                
                # Highlight important elements
                if el_info['important']:
                    grid[height - 1 - grid_y][grid_x] = self.colorize(atom_symbol, color)
                else:
                    grid[height - 1 - grid_y][grid_x] = atom_symbol
        
        print(f"    {layer_title} (top view):")
        for row in grid:
            print(f"    {''.join(row)}")
    
    def track_deposited_species(self, atoms_list, titles=None):
        """
        Track evolution of deposited species (Al, C) over cycles
        
        Args:
            atoms_list: List of Atoms objects
            titles: List of titles for each structure
        """
        print(f"\n{self.colorize('Deposited Species Evolution', 'bold')}")
        print("═" * 70)
        
        if titles is None:
            titles = [f"Structure {i+1}" for i in range(len(atoms_list))]
        
        print("Structure                   │ Al  │ C   │ Al:C │ Coverage")
        print("────────────────────────────┼─────┼─────┼──────┼─────────")
        
        for atoms, title in zip(atoms_list, titles):
            surface_atoms, info = self.extract_surface_region(atoms)
            composition = info['composition']
            
            al_count = composition.get('Al', 0)
            c_count = composition.get('C', 0)
            
            # Calculate Al:C ratio
            if al_count > 0 and c_count > 0:
                ratio = c_count / al_count
                ratio_str = f"{ratio:.1f}"
            else:
                ratio_str = "  -  "
            
            # Estimate coverage (very rough)
            si_count = composition.get('Si', 0)
            if si_count > 0:
                coverage = (al_count / si_count) * 100
                coverage_str = f"{coverage:.1f}%"
            else:
                coverage_str = "  -  "
            
            # Color code the species counts
            al_display = self.colorize(f"{al_count:3}", 'blue') if al_count > 0 else f"{al_count:3}"
            c_display = self.colorize(f"{c_count:3}", 'yellow') if c_count > 0 else f"{c_count:3}"
            
            print(f"{title[:27]:27} │ {al_display} │ {c_display} │ {ratio_str:4} │ {coverage_str:7}")
    
    def show_mld_surface_evolution(self):
        """Show evolution of surface during MLD cycles"""
        print(f"\n{self.colorize('MLD Surface Evolution Analysis', 'bold')}")
        print("═" * 80)
        
        # Find cycle files
        cycle_files = sorted(glob.glob('mld_cycle*_after_*.xyz'))
        
        if not cycle_files:
            print(f"{self.colorize('❌ No MLD cycle files found!', 'red')}")
            return
        
        structures = []
        titles = []
        
        for cycle_file in cycle_files:
            try:
                atoms = read(cycle_file)
                structures.append(atoms)
                
                # Create descriptive title
                filename = os.path.basename(cycle_file)
                cycle_num = filename.split('cycle')[1].split('_')[0]
                phase = 'TMA' if 'tma' in filename else 'H2O'
                titles.append(f"Cycle {cycle_num} after {phase}")
                
            except Exception as e:
                print(f"Warning: Could not read {cycle_file}")
        
        if not structures:
            print(f"{self.colorize('❌ No valid structures found!', 'red')}")
            return
        
        # Track deposited species
        self.track_deposited_species(structures, titles)
        
        # Show detailed analysis for key structures
        print(f"\n{self.colorize('Detailed Surface Analysis:', 'bold')}")
        
        # Show initial structure
        if len(structures) > 0:
            print(f"\n{self.colorize('═ Initial State ═', 'cyan')}")
            self.show_surface_layers(structures[0], titles[0])
        
        # Show final structure
        if len(structures) > 1:
            print(f"\n{self.colorize('═ Final State ═', 'green')}")
            self.show_surface_layers(structures[-1], titles[-1])
        
        # Show key intermediate if available
        if len(structures) > 2:
            mid_idx = len(structures) // 2
            print(f"\n{self.colorize('═ Intermediate State ═', 'yellow')}")
            self.show_surface_layers(structures[mid_idx], titles[mid_idx])
    
    def analyze_surface_chemistry(self):
        """Analyze chemical groups and bonding on surface"""
        print(f"\n{self.colorize('Surface Chemistry Analysis', 'bold')}")
        print("═" * 70)
        
        final_file = 'mld_final_structure.xyz'
        if not os.path.exists(final_file):
            print(f"{self.colorize('❌ No final structure found!', 'red')}")
            return
        
        atoms = read(final_file)
        surface_atoms, info = self.extract_surface_region(atoms)
        
        positions = surface_atoms.get_positions()
        symbols = surface_atoms.get_chemical_symbols()
        
        print(f"Surface atoms analyzed: {len(surface_atoms)}")
        print(f"Surface composition: {dict(info['composition'])}")
        
        # Identify chemical groups
        print(f"\n{self.colorize('Chemical Group Identification:', 'bold')}")
        
        # Find Al-C bonds (methyl groups)
        al_indices = [i for i, sym in enumerate(symbols) if sym == 'Al']
        c_indices = [i for i, sym in enumerate(symbols) if sym == 'C']
        
        methyl_groups = 0
        for al_idx in al_indices:
            al_pos = positions[al_idx]
            bonded_carbons = 0
            
            for c_idx in c_indices:
                c_pos = positions[c_idx]
                distance = np.linalg.norm(al_pos - c_pos)
                
                if 1.8 < distance < 2.1:  # Al-C bond length
                    bonded_carbons += 1
            
            methyl_groups += bonded_carbons
        
        # Find OH groups
        o_indices = [i for i, sym in enumerate(symbols) if sym == 'O']
        h_indices = [i for i, sym in enumerate(symbols) if sym == 'H']
        
        oh_groups = 0
        for o_idx in o_indices:
            o_pos = positions[o_idx]
            
            for h_idx in h_indices:
                h_pos = positions[h_idx]
                distance = np.linalg.norm(o_pos - h_pos)
                
                if 0.85 < distance < 1.15:  # O-H bond length
                    oh_groups += 1
                    break  # Each O can only bond to one H
        
        print(f"  {self.colorize('OH groups:', 'cyan')} {oh_groups}")
        print(f"  {self.colorize('Al-CH3 groups:', 'blue')} {methyl_groups}")
        print(f"  {self.colorize('Al atoms:', 'blue')} {len(al_indices)}")
        
        # Calculate average coordination
        if len(al_indices) > 0:
            avg_coordination = methyl_groups / len(al_indices)
            print(f"  Average Al coordination: {avg_coordination:.1f} CH3 groups per Al")
        
        # Surface density calculations
        if len(al_indices) > 0:
            # Rough surface area estimate
            x_range = positions[:, 0].max() - positions[:, 0].min()
            y_range = positions[:, 1].max() - positions[:, 1].min()
            surface_area = x_range * y_range * 1e-2  # Convert to nm²
            
            al_density = len(al_indices) / surface_area
            oh_density = oh_groups / surface_area
            
            print(f"\n{self.colorize('Surface Densities:', 'bold')}")
            print(f"  Al density: {al_density:.2f} Al/nm²")
            print(f"  OH density: {oh_density:.2f} OH/nm²")

def main():
    """Interactive surface visualization menu"""
    visualizer = SurfaceVisualizer()
    
    # Check for MLD files
    mld_files = glob.glob('mld_*.xyz')
    if not mld_files:
        print(f"{visualizer.colorize('❌ No MLD files found!', 'red')}")
        print("Run step5_realistic_mld.py first")
        return 1
    
    while True:
        print(f"\n{visualizer.colorize('═' * 60, 'cyan')}")
        print(f"{visualizer.colorize('SURFACE-FOCUSED MLD VISUALIZER', 'bold')}")
        print(f"{visualizer.colorize('═' * 60, 'cyan')}")
        
        print(f"\n{visualizer.colorize('Options:', 'bold')}")
        print(f"  {visualizer.colorize('1', 'cyan')} - Show MLD Surface Evolution")
        print(f"  {visualizer.colorize('2', 'cyan')} - Analyze Surface Chemistry")
        print(f"  {visualizer.colorize('3', 'cyan')} - Show Final Structure Layers")
        print(f"  {visualizer.colorize('4', 'cyan')} - Compare Before/After Surface")
        print(f"  {visualizer.colorize('5', 'cyan')} - Settings (Surface Depth)")
        print(f"  {visualizer.colorize('q', 'red')} - Quit")
        
        choice = input(f"\n{visualizer.colorize('Choice:', 'bold')} ").strip().lower()
        
        if choice == '1':
            visualizer.show_mld_surface_evolution()
        elif choice == '2':
            visualizer.analyze_surface_chemistry()
        elif choice == '3':
            final_file = 'mld_final_structure.xyz'
            if os.path.exists(final_file):
                atoms = read(final_file)
                visualizer.show_surface_layers(atoms, "Final MLD Structure")
            else:
                print(f"{visualizer.colorize('❌ No final structure found!', 'red')}")
        elif choice == '4':
            # Compare initial and final surfaces
            initial_files = ['uv_ozone_si_surface.xyz', 'surface_si100_hydroxylated.xyz']
            initial_atoms = None
            
            for init_file in initial_files:
                if os.path.exists(init_file):
                    initial_atoms = read(init_file)
                    break
            
            final_file = 'mld_final_structure.xyz'
            if initial_atoms and os.path.exists(final_file):
                final_atoms = read(final_file)
                
                print(f"\n{visualizer.colorize('═ BEFORE MLD ═', 'yellow')}")
                visualizer.show_surface_layers(initial_atoms, "Initial Surface")
                
                print(f"\n{visualizer.colorize('═ AFTER MLD ═', 'green')}")
                visualizer.show_surface_layers(final_atoms, "Final Surface")
                
                # Track changes
                structures = [initial_atoms, final_atoms]
                titles = ["Initial", "Final"]
                visualizer.track_deposited_species(structures, titles)
            else:
                print(f"{visualizer.colorize('❌ Missing initial or final structure!', 'red')}")
        elif choice == '5':
            current_depth = visualizer.surface_depth
            print(f"Current surface depth: {current_depth:.1f} Å")
            try:
                new_depth = float(input("Enter new surface depth (Å): "))
                visualizer.surface_depth = new_depth
                print(f"Surface depth set to {new_depth:.1f} Å")
            except ValueError:
                print("Invalid input, keeping current setting")
        elif choice == 'q':
            print(f"\n{visualizer.colorize('Thanks for using Surface Visualizer!', 'green')}")
            break
        else:
            print(f"{visualizer.colorize('Invalid choice!', 'red')}")
        
        input(f"\n{visualizer.colorize('Press Enter to continue...', 'yellow')}")
    
    return 0

if __name__ == "__main__":
    exit(main())