#!/usr/bin/env python3
"""
ASCII Structure Viewer for MLD
Creates detailed ASCII art representations of molecular structures
Perfect for terminal viewing of surface chemistry
"""

import numpy as np
from ase.io import read
import os

class ASCIIStructureViewer:
    """Create ASCII art molecular structures"""
    
    def __init__(self):
        # Enhanced element representation
        self.elements = {
            'Si': {'char': '●', 'color': '\033[90m', 'size': 'large'},      # Dark gray, large
            'O':  {'char': '○', 'color': '\033[91m', 'size': 'medium'},     # Red, medium  
            'H':  {'char': '·', 'color': '\033[97m', 'size': 'small'},      # White, small
            'Al': {'char': '◆', 'color': '\033[94m', 'size': 'large'},      # Blue, large
            'C':  {'char': '◯', 'color': '\033[93m', 'size': 'medium'},     # Yellow, medium
            'N':  {'char': '◉', 'color': '\033[95m', 'size': 'medium'}      # Magenta, medium
        }
        
        self.reset = '\033[0m'
        self.bold = '\033[1m'
        
        # Bond representations
        self.bonds = {
            'single': '─',
            'double': '═',
            'vertical': '│',
            'diagonal_up': '╱',
            'diagonal_down': '╲'
        }
    
    def get_element_display(self, element):
        """Get colored display character for element"""
        if element in self.elements:
            info = self.elements[element]
            return f"{info['color']}{info['char']}{self.reset}"
        else:
            return f"\033[37m?{self.reset}"  # White question mark for unknown
    
    def create_surface_cross_section(self, atoms, title="Surface Cross-Section"):
        """
        Create a cross-sectional view of the surface (side view)
        Shows Z vs X coordinates to see surface layers
        """
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        
        if len(positions) == 0:
            print("No atoms to display")
            return
        
        print(f"\n{self.bold}{title}{self.reset}")
        print("Side view (X vs Z coordinates)")
        print("─" * 60)
        
        # Use X and Z coordinates
        x_coords = positions[:, 0]
        z_coords = positions[:, 2]
        
        # Create ASCII grid
        width, height = 60, 20
        
        x_min, x_max = x_coords.min(), x_coords.max()
        z_min, z_max = z_coords.min(), z_coords.max()
        
        x_range = x_max - x_min if x_max > x_min else 1
        z_range = z_max - z_min if z_max > z_min else 1
        
        # Create grid
        grid = [[' ' for _ in range(width)] for _ in range(height)]
        
        # Place atoms
        for x, z, symbol in zip(x_coords, z_coords, symbols):
            grid_x = int(((x - x_min) / x_range) * (width - 1))
            grid_z = int(((z - z_min) / z_range) * (height - 1))
            
            if 0 <= grid_x < width and 0 <= grid_z < height:
                # Flip Z for proper display (higher Z at top)
                display_z = height - 1 - grid_z
                
                # Handle overlapping atoms
                if grid[display_z][grid_x] == ' ':
                    grid[display_z][grid_x] = self.get_element_display(symbol)
                else:
                    # Multiple atoms at same position - show as '+'
                    grid[display_z][grid_x] = f"\033[96m+{self.reset}"
        
        # Print with Z-axis labels
        for i, row in enumerate(grid):
            z_value = z_min + ((height - 1 - i) / (height - 1)) * z_range
            print(f"{z_value:6.1f} │{''.join(row)}")
        
        # X-axis
        x_axis = '─' * width
        print(f"{'':8}└{x_axis}")
        
        # X-axis labels
        x_label_line = '         '
        for i in range(0, width, 10):
            x_value = x_min + (i / (width - 1)) * x_range
            x_label_line += f"{x_value:5.1f}     "
        print(x_label_line[:width + 9])
        
        # Legend
        print(f"\n{self.bold}Legend:{self.reset}")
        unique_elements = sorted(set(symbols))
        legend_line = ""
        for element in unique_elements:
            count = symbols.count(element)
            legend_line += f"{self.get_element_display(element)}={element}({count}) "
        print(legend_line)
    
    def create_surface_top_view(self, atoms, title="Surface Top View", 
                              surface_only=True, depth=5.0):
        """
        Create top-down view of surface (X vs Y)
        
        Args:
            atoms: ASE Atoms object
            title: Display title
            surface_only: If True, show only top 'depth' Å
            depth: Depth in Å to consider as surface
        """
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        
        if len(positions) == 0:
            print("No atoms to display")
            return
        
        # Filter to surface atoms if requested
        if surface_only:
            z_coords = positions[:, 2]
            z_max = z_coords.max()
            z_surface = z_max - depth
            
            surface_mask = z_coords >= z_surface
            positions = positions[surface_mask]
            symbols = [symbols[i] for i in range(len(symbols)) if surface_mask[i]]
        
        print(f"\n{self.bold}{title}{self.reset}")
        if surface_only:
            print(f"Top view - surface layer only (top {depth:.1f} Å)")
        else:
            print("Top view (X vs Y coordinates)")
        print("─" * 60)
        
        if len(positions) == 0:
            print("No surface atoms found")
            return
        
        # Use X and Y coordinates
        x_coords = positions[:, 0]
        y_coords = positions[:, 1]
        
        # Create ASCII grid
        width, height = 50, 25
        
        x_min, x_max = x_coords.min(), x_coords.max()
        y_min, y_max = y_coords.min(), y_coords.max()
        
        x_range = x_max - x_min if x_max > x_min else 1
        y_range = y_max - y_min if y_max > y_min else 1
        
        # Create grid
        grid = [[' ' for _ in range(width)] for _ in range(height)]
        
        # Place atoms
        for x, y, symbol in zip(x_coords, y_coords, symbols):
            grid_x = int(((x - x_min) / x_range) * (width - 1))
            grid_y = int(((y - y_min) / y_range) * (height - 1))
            
            if 0 <= grid_x < width and 0 <= grid_y < height:
                # Flip Y for proper display
                display_y = height - 1 - grid_y
                
                if grid[display_y][grid_x] == ' ':
                    grid[display_y][grid_x] = self.get_element_display(symbol)
                else:
                    # Overlapping atoms
                    grid[display_y][grid_x] = f"\033[96m+{self.reset}"
        
        # Print grid
        for row in grid:
            print(''.join(row))
        
        # Legend and stats
        print(f"\n{self.bold}Surface Composition:{self.reset}")
        from collections import Counter
        composition = Counter(symbols)
        for element, count in composition.most_common():
            print(f"  {self.get_element_display(element)} {element}: {count} atoms")
    
    def show_mld_chemistry_focus(self, atoms, title="MLD Chemistry Focus"):
        """
        Focus on the interesting chemistry: Al, C, OH groups
        De-emphasize bulk Si
        """
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        
        print(f"\n{self.bold}{title}{self.reset}")
        print("Highlighting: Al (blue ◆), C (yellow ◯), OH groups (red ○, white ·)")
        print("─" * 70)
        
        # Find interesting atoms (not bulk Si)
        interesting_indices = []
        for i, symbol in enumerate(symbols):
            if symbol in ['Al', 'C', 'O', 'H']:
                interesting_indices.append(i)
            elif symbol == 'Si':
                # Include Si only if it's likely surface Si
                z = positions[i][2]
                z_max = positions[:, 2].max()
                if z > z_max - 3.0:  # Top 3 Å only
                    interesting_indices.append(i)
        
        if not interesting_indices:
            print("No interesting surface chemistry found")
            return
        
        # Extract interesting atoms
        interesting_pos = positions[interesting_indices]
        interesting_syms = [symbols[i] for i in interesting_indices]
        
        print(f"Showing {len(interesting_indices)} chemically relevant atoms")
        print(f"(Hiding {len(symbols) - len(interesting_indices)} bulk Si atoms)")
        
        # Create focused visualization
        width, height = 60, 20
        
        x_coords = interesting_pos[:, 0]
        y_coords = interesting_pos[:, 1]
        z_coords = interesting_pos[:, 2]
        
        # Use Z vs X for side view of chemistry
        x_min, x_max = x_coords.min(), x_coords.max()
        z_min, z_max = z_coords.min(), z_coords.max()
        
        x_range = x_max - x_min if x_max > x_min else 1
        z_range = z_max - z_min if z_max > z_min else 1
        
        grid = [[' ' for _ in range(width)] for _ in range(height)]
        
        for x, z, symbol in zip(x_coords, z_coords, interesting_syms):
            grid_x = int(((x - x_min) / x_range) * (width - 1))
            grid_z = int(((z - z_min) / z_range) * (height - 1))
            
            if 0 <= grid_x < width and 0 <= grid_z < height:
                display_z = height - 1 - grid_z
                
                if grid[display_z][grid_x] == ' ':
                    grid[display_z][grid_x] = self.get_element_display(symbol)
                else:
                    grid[display_z][grid_x] = f"\033[96m+{self.reset}"
        
        # Print with labels
        for i, row in enumerate(grid):
            z_value = z_min + ((height - 1 - i) / (height - 1)) * z_range
            print(f"{z_value:6.1f} │{''.join(row)}")
        
        print(f"{'':8}└{'─' * width}")
        
        # Chemistry analysis
        print(f"\n{self.bold}Chemical Group Analysis:{self.reset}")
        
        # Count specific groups
        al_count = interesting_syms.count('Al')
        c_count = interesting_syms.count('C')
        o_count = interesting_syms.count('O')
        h_count = interesting_syms.count('H')
        
        print(f"  Al atoms (TMA centers): {self.get_element_display('Al')} × {al_count}")
        print(f"  C atoms (methyl groups): {self.get_element_display('C')} × {c_count}")
        print(f"  O atoms (surface oxide): {self.get_element_display('O')} × {o_count}")
        print(f"  H atoms (OH, CH3): {self.get_element_display('H')} × {h_count}")
        
        if al_count > 0:
            c_per_al = c_count / al_count
            print(f"\n  Average C per Al: {c_per_al:.1f} (ideal TMA = 3.0)")
            
            if c_per_al < 2.5:
                print(f"  {self.bold}→ Significant methyl loss (hydrolysis/pyrolysis){self.reset}")
            elif c_per_al > 2.8:
                print(f"  {self.bold}→ Good methyl retention{self.reset}")
    
    def compare_mld_progression(self, file_list, titles=None):
        """
        Show progression of MLD cycles side by side
        """
        if not file_list:
            print("No files provided")
            return
        
        if titles is None:
            titles = [f"Step {i+1}" for i in range(len(file_list))]
        
        print(f"\n{self.bold}MLD Progression Comparison{self.reset}")
        print("═" * 80)
        
        structures = []
        for filename in file_list:
            if os.path.exists(filename):
                atoms = read(filename)
                structures.append(atoms)
            else:
                print(f"Warning: {filename} not found")
        
        # Show chemistry focus for each
        for atoms, title in zip(structures, titles):
            self.show_mld_chemistry_focus(atoms, f"{title} - Chemistry")
            print()

def main():
    """Interactive ASCII structure viewer"""
    viewer = ASCIIStructureViewer()
    
    print(f"{viewer.bold}ASCII Structure Viewer for MLD{viewer.reset}")
    print("═" * 50)
    
    while True:
        print(f"\n{viewer.bold}Visualization Options:{viewer.reset}")
        print("1. Surface cross-section (side view)")
        print("2. Surface top view") 
        print("3. MLD chemistry focus")
        print("4. Compare MLD progression")
        print("5. Analyze specific file")
        print("q. Quit")
        
        choice = input(f"\n{viewer.bold}Choice:{viewer.reset} ").strip()
        
        if choice == '1':
            # Surface cross-section
            final_file = 'mld_final_structure.xyz'
            if os.path.exists(final_file):
                atoms = read(final_file)
                viewer.create_surface_cross_section(atoms, "MLD Final Structure - Cross Section")
            else:
                print("No final structure found")
                
        elif choice == '2':
            # Surface top view
            final_file = 'mld_final_structure.xyz'
            if os.path.exists(final_file):
                atoms = read(final_file)
                viewer.create_surface_top_view(atoms, "MLD Final Structure - Top View")
            else:
                print("No final structure found")
                
        elif choice == '3':
            # Chemistry focus
            final_file = 'mld_final_structure.xyz'
            if os.path.exists(final_file):
                atoms = read(final_file)
                viewer.show_mld_chemistry_focus(atoms, "MLD Surface Chemistry")
            else:
                print("No final structure found")
                
        elif choice == '4':
            # Compare progression
            import glob
            cycle_files = sorted(glob.glob('mld_cycle*_after_h2o.xyz'))
            if cycle_files:
                titles = []
                for f in cycle_files:
                    cycle_num = f.split('cycle')[1].split('_')[0]
                    titles.append(f"Cycle {cycle_num}")
                
                viewer.compare_mld_progression(cycle_files[:4], titles[:4])  # Show first 4
            else:
                print("No cycle files found")
                
        elif choice == '5':
            # Analyze specific file
            filename = input("Enter filename: ").strip()
            if os.path.exists(filename):
                atoms = read(filename)
                print(f"\nAnalyzing {filename}...")
                viewer.create_surface_cross_section(atoms, f"{filename} - Cross Section")
                viewer.show_mld_chemistry_focus(atoms, f"{filename} - Chemistry")
            else:
                print(f"File {filename} not found")
                
        elif choice == 'q':
            print(f"\n{viewer.bold}Goodbye!{viewer.reset}")
            break
        else:
            print("Invalid choice")
        
        input(f"\n{viewer.bold}Press Enter to continue...{viewer.reset}")

if __name__ == "__main__":
    main()