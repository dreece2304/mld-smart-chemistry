#!/usr/bin/env python3
"""
Terminal-based MLD Visualization and Analysis Tool
Perfect for Ubuntu servers - no GUI required!

Features:
- ASCII molecular structure display
- Terminal plots using Unicode characters
- Real-time data analysis
- Interactive menus
- Color-coded output
"""

import os
import glob
import numpy as np
from ase.io import read
import matplotlib.pyplot as plt
from collections import Counter
import sys

class TerminalVisualizer:
    """Terminal-based visualization for MLD simulations"""
    
    def __init__(self):
        # Terminal colors
        self.colors = {
            'red': '\033[91m',
            'green': '\033[92m',
            'yellow': '\033[93m',
            'blue': '\033[94m',
            'magenta': '\033[95m',
            'cyan': '\033[96m',
            'white': '\033[97m',
            'bold': '\033[1m',
            'end': '\033[0m'
        }
        
        # Element symbols for visualization
        self.element_symbols = {
            'Si': '●', 'O': '○', 'H': '·', 'Al': '◆', 'C': '◯', 'N': '◉'
        }
        
        # Element colors
        self.element_colors = {
            'Si': 'cyan', 'O': 'red', 'H': 'white', 
            'Al': 'blue', 'C': 'yellow', 'N': 'magenta'
        }
    
    def colorize(self, text, color):
        """Add color to text"""
        return f"{self.colors[color]}{text}{self.colors['end']}"
    
    def print_header(self, title):
        """Print a fancy header"""
        width = max(60, len(title) + 10)
        border = "=" * width
        
        print(f"\n{self.colorize(border, 'cyan')}")
        print(f"{self.colorize(title.center(width), 'bold')}")
        print(f"{self.colorize(border, 'cyan')}\n")
    
    def plot_ascii_graph(self, data, labels=None, title="", width=60, height=15):
        """Create ASCII bar chart"""
        if not data:
            print("No data to plot")
            return
        
        max_val = max(data)
        if max_val == 0:
            print("All values are zero")
            return
        
        print(f"\n{self.colorize(title, 'bold')}")
        print("─" * width)
        
        # Scale data to fit height
        scaled_data = [int((val / max_val) * height) for val in data]
        
        # Print graph from top to bottom
        for row in range(height, 0, -1):
            line = ""
            for i, val in enumerate(scaled_data):
                if val >= row:
                    line += "█"
                else:
                    line += " "
                line += " "
            
            # Add y-axis values
            y_val = (row / height) * max_val
            print(f"{y_val:6.1f} │{line}")
        
        # X-axis
        x_axis = "─" * (len(data) * 2)
        print(f"{'':8}└{x_axis}")
        
        # Labels
        if labels:
            label_line = "        "
            for label in labels:
                label_line += f"{label[:1]} "
            print(label_line)
    
    def plot_ascii_line(self, data, title="", width=60, height=10):
        """Create ASCII line plot"""
        if len(data) < 2:
            print("Need at least 2 data points")
            return
        
        print(f"\n{self.colorize(title, 'bold')}")
        print("─" * width)
        
        max_val = max(data)
        min_val = min(data)
        range_val = max_val - min_val
        
        if range_val == 0:
            print("Data has no variation")
            return
        
        # Scale data to height
        scaled_data = [int(((val - min_val) / range_val) * (height - 1)) for val in data]
        
        # Create plot grid
        grid = [[' ' for _ in range(len(data))] for _ in range(height)]
        
        # Plot points and lines
        for i, y in enumerate(scaled_data):
            grid[height - 1 - y][i] = '●'
            
            # Connect with lines
            if i > 0:
                y_prev = scaled_data[i-1]
                y_curr = y
                
                if y_prev != y_curr:
                    # Draw connecting line
                    start_y = min(y_prev, y_curr)
                    end_y = max(y_prev, y_curr)
                    
                    for y_line in range(start_y, end_y + 1):
                        if grid[height - 1 - y_line][i] == ' ':
                            grid[height - 1 - y_line][i] = '│'
        
        # Print grid
        for row_idx, row in enumerate(grid):
            y_val = min_val + ((height - 1 - row_idx) / (height - 1)) * range_val
            line = ''.join(row)
            print(f"{y_val:6.1f} │{line}")
        
        # X-axis
        print(f"{'':8}└{'─' * len(data)}")
        
        # X labels
        x_labels = ""
        for i in range(len(data)):
            x_labels += str(i % 10)
        print(f"{'':9}{x_labels}")
    
    def show_molecule_structure(self, atoms, title="Molecule"):
        """Display 2D projection of molecular structure"""
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        
        print(f"\n{self.colorize(title, 'bold')}")
        print(f"Total atoms: {len(atoms)}")
        
        # Project to 2D (use X-Y plane)
        x_coords = positions[:, 0]
        y_coords = positions[:, 1]
        
        # Scale to terminal size
        width, height = 50, 20
        
        if len(x_coords) > 0:
            x_min, x_max = x_coords.min(), x_coords.max()
            y_min, y_max = y_coords.min(), y_coords.max()
            
            x_range = x_max - x_min
            y_range = y_max - y_min
            
            if x_range == 0:
                x_range = 1
            if y_range == 0:
                y_range = 1
            
            # Create display grid
            grid = [[' ' for _ in range(width)] for _ in range(height)]
            
            # Place atoms on grid
            for i, (x, y, symbol) in enumerate(zip(x_coords, y_coords, symbols)):
                # Scale to grid coordinates
                grid_x = int(((x - x_min) / x_range) * (width - 1))
                grid_y = int(((y - y_min) / y_range) * (height - 1))
                
                # Get symbol and color
                atom_symbol = self.element_symbols.get(symbol, '?')
                color = self.element_colors.get(symbol, 'white')
                
                # Place on grid (invert Y for proper display)
                if 0 <= grid_x < width and 0 <= grid_y < height:
                    grid[height - 1 - grid_y][grid_x] = self.colorize(atom_symbol, color)
            
            # Print grid
            for row in grid:
                print(''.join(row))
            
            # Legend
            print(f"\n{self.colorize('Legend:', 'bold')}")
            unique_elements = set(symbols)
            for element in sorted(unique_elements):
                symbol = self.element_symbols.get(element, '?')
                color = self.element_colors.get(element, 'white')
                count = symbols.count(element)
                print(f"  {self.colorize(symbol, color)} = {element} ({count} atoms)")
    
    def analyze_mld_cycles(self):
        """Analyze MLD cycle progression"""
        self.print_header("MLD Cycle Analysis")
        
        # Find cycle files
        cycle_files = sorted(glob.glob('mld_cycle*_after_*.xyz'))
        
        if not cycle_files:
            print(f"{self.colorize('❌ No MLD cycle files found!', 'red')}")
            print("Run step5_realistic_mld.py first")
            return
        
        print(f"{self.colorize('✅ Found', 'green')} {len(cycle_files)} cycle files\n")
        
        # Analyze each cycle
        cycle_data = {}
        
        for cycle_file in cycle_files:
            try:
                atoms = read(cycle_file)
                filename = os.path.basename(cycle_file)
                
                # Extract cycle number and phase
                cycle_num = int(filename.split('cycle')[1].split('_')[0])
                phase = 'TMA' if 'tma' in filename else 'H2O'
                
                if cycle_num not in cycle_data:
                    cycle_data[cycle_num] = {}
                
                cycle_data[cycle_num][phase] = {
                    'atoms': len(atoms),
                    'composition': Counter(atoms.get_chemical_symbols())
                }
                
            except Exception as e:
                print(f"Error reading {cycle_file}: {e}")
        
        # Display cycle progression
        print(f"{self.colorize('Cycle Progression:', 'bold')}")
        print("─" * 60)
        
        cycles = sorted(cycle_data.keys())
        atom_counts = []
        
        for cycle in cycles:
            print(f"\n{self.colorize(f'Cycle {cycle}:', 'yellow')}")
            
            for phase in ['TMA', 'H2O']:
                if phase in cycle_data[cycle]:
                    data = cycle_data[cycle][phase]
                    atoms = data['atoms']
                    comp = data['composition']
                    
                    print(f"  After {phase:3}: {atoms:4} atoms")
                    
                    # Show composition changes
                    comp_str = ", ".join([f"{el}:{count}" for el, count in comp.most_common()])
                    print(f"           ({comp_str})")
                    
                    if phase == 'H2O':  # End of cycle
                        atom_counts.append(atoms)
        
        # Plot growth
        if len(atom_counts) > 1:
            print(f"\n{self.colorize('Growth Trend:', 'bold')}")
            self.plot_ascii_line(atom_counts, "Total Atoms vs Cycle", height=8)
            
            # Calculate growth per cycle
            growth = [atom_counts[i] - atom_counts[i-1] for i in range(1, len(atom_counts))]
            if growth:
                print(f"\n{self.colorize('Growth per Cycle:', 'bold')}")
                self.plot_ascii_graph(growth, [f"C{i+2}" for i in range(len(growth))], 
                                    "Atoms Added per Cycle", height=6)
    
    def show_composition_analysis(self):
        """Show detailed composition analysis"""
        self.print_header("Composition Analysis")
        
        # Find final structure
        final_file = 'mld_final_structure.xyz'
        if not os.path.exists(final_file):
            print(f"{self.colorize('❌ No final structure found!', 'red')}")
            return
        
        atoms = read(final_file)
        symbols = atoms.get_chemical_symbols()
        composition = Counter(symbols)
        
        total_atoms = len(atoms)
        
        print(f"{self.colorize('Final Structure Composition:', 'bold')}")
        print(f"Total atoms: {total_atoms}\n")
        
        # Create composition table
        print("Element │ Count │ Percentage │ Visual")
        print("────────┼───────┼────────────┼──────────────────────")
        
        elements = sorted(composition.keys())
        counts = [composition[el] for el in elements]
        
        for element, count in zip(elements, counts):
            percentage = (count / total_atoms) * 100
            
            # Visual bar
            bar_length = int((count / max(counts)) * 20)
            color = self.element_colors.get(element, 'white')
            symbol = self.element_symbols.get(element, '?')
            bar = self.colorize(symbol * bar_length, color)
            
            print(f"{element:7} │ {count:5} │ {percentage:8.1f}% │ {bar}")
        
        # ASCII pie chart approximation
        print(f"\n{self.colorize('Composition Distribution:', 'bold')}")
        total_chars = 50
        
        for element, count in zip(elements, counts):
            chars = int((count / total_atoms) * total_chars)
            color = self.element_colors.get(element, 'white')
            symbol = self.element_symbols.get(element, '?')
            
            if chars > 0:
                bar = self.colorize(symbol * chars, color)
                print(f"{element}: {bar} ({count} atoms)")
    
    def compare_structures(self):
        """Compare initial and final structures"""
        self.print_header("Structure Comparison")
        
        # Find initial surface
        initial_files = ['uv_ozone_si_surface.xyz', 'surface_si100_hydroxylated.xyz', 'simple_surface.xyz']
        initial_surface = None
        initial_name = None
        
        for init_file in initial_files:
            if os.path.exists(init_file):
                initial_surface = read(init_file)
                initial_name = init_file
                break
        
        if not initial_surface:
            print(f"{self.colorize('❌ No initial surface found!', 'red')}")
            return
        
        # Find final surface
        final_file = 'mld_final_structure.xyz'
        if not os.path.exists(final_file):
            print(f"{self.colorize('❌ No final structure found!', 'red')}")
            return
        
        final_surface = read(final_file)
        
        # Show comparison
        print(f"{self.colorize('Before MLD:', 'yellow')}")
        self.show_molecule_structure(initial_surface, f"Initial: {os.path.basename(initial_name)}")
        
        print(f"\n{self.colorize('After MLD:', 'green')}")
        self.show_molecule_structure(final_surface, "Final: MLD Structure")
        
        # Growth statistics
        initial_atoms = len(initial_surface)
        final_atoms = len(final_surface)
        growth = final_atoms - initial_atoms
        growth_percent = (growth / initial_atoms) * 100
        
        print(f"\n{self.colorize('Growth Summary:', 'bold')}")
        print(f"Initial atoms:    {initial_atoms}")
        print(f"Final atoms:      {final_atoms}")
        print(f"Atoms added:      {self.colorize(str(growth), 'green')}")
        print(f"Growth:           {self.colorize(f'{growth_percent:.1f}%', 'green')}")
    
    def interactive_menu(self):
        """Interactive menu system"""
        while True:
            self.print_header("MLD Terminal Visualizer")
            
            print(f"{self.colorize('Available Options:', 'bold')}")
            print(f"  {self.colorize('1', 'cyan')} - Analyze MLD Cycles")
            print(f"  {self.colorize('2', 'cyan')} - Show Composition Analysis")
            print(f"  {self.colorize('3', 'cyan')} - Compare Before/After Structures")
            print(f"  {self.colorize('4', 'cyan')} - Show Final Structure")
            print(f"  {self.colorize('5', 'cyan')} - List All Files")
            print(f"  {self.colorize('q', 'red')} - Quit")
            
            choice = input(f"\n{self.colorize('Enter choice:', 'bold')} ").strip().lower()
            
            if choice == '1':
                self.analyze_mld_cycles()
            elif choice == '2':
                self.show_composition_analysis()
            elif choice == '3':
                self.compare_structures()
            elif choice == '4':
                final_file = 'mld_final_structure.xyz'
                if os.path.exists(final_file):
                    atoms = read(final_file)
                    self.show_molecule_structure(atoms, "Final MLD Structure")
                else:
                    print(f"{self.colorize('❌ No final structure found!', 'red')}")
            elif choice == '5':
                self.list_files()
            elif choice == 'q':
                print(f"\n{self.colorize('Thanks for using MLD Terminal Visualizer!', 'green')}")
                break
            else:
                print(f"{self.colorize('Invalid choice!', 'red')}")
            
            input(f"\n{self.colorize('Press Enter to continue...', 'yellow')}")
    
    def list_files(self):
        """List all relevant files"""
        self.print_header("Available Files")
        
        file_patterns = [
            ('MLD Cycles', 'mld_cycle*.xyz'),
            ('Final Structure', 'mld_final_structure.xyz'),
            ('Surfaces', '*surface*.xyz'),
            ('Molecules', '*optimized.xyz'),
            ('DFT Results', 'dft_*.xyz')
        ]
        
        for category, pattern in file_patterns:
            files = glob.glob(pattern)
            if files:
                print(f"\n{self.colorize(category + ':', 'bold')}")
                for f in sorted(files):
                    size = os.path.getsize(f) / 1024  # KB
                    print(f"  📄 {f} ({size:.1f} KB)")
            else:
                print(f"\n{self.colorize(category + ':', 'bold')} {self.colorize('No files found', 'yellow')}")

def main():
    """Main function"""
    visualizer = TerminalVisualizer()
    
    # Check if we have MLD results
    mld_files = glob.glob('mld_*.xyz')
    
    if not mld_files:
        visualizer.print_header("MLD Terminal Visualizer")
        print(f"{visualizer.colorize('❌ No MLD results found!', 'red')}")
        print("Please run one of these first:")
        print("  • python scripts/step5_realistic_mld.py")
        print("  • python scripts/step5_true_dft_mld.py")
        return 1
    
    # Start interactive menu
    try:
        visualizer.interactive_menu()
    except KeyboardInterrupt:
        print(f"\n\n{visualizer.colorize('Goodbye!', 'green')}")
    
    return 0

if __name__ == "__main__":
    exit(main())