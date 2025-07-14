#!/usr/bin/env python3
"""
Bulk structure analysis tools for MLD films
Analyzes crystallinity, density, mechanical properties, and layer thickness
"""

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.neighborlist import NeighborList
from ase.geometry import get_distances
import matplotlib.pyplot as plt
from pymatgen.core import Structure
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.analysis.bond_valence import BVAnalyzer
import spglib

class BulkMLDAnalyzer:
    """Comprehensive analysis of bulk MLD film properties"""
    
    def __init__(self, structure_file):
        self.atoms = read(structure_file)
        self.structure = Structure.from_file(structure_file)
        
    def analyze_crystallinity(self):
        """Determine crystal structure and symmetry"""
        analyzer = SpacegroupAnalyzer(self.structure)
        
        results = {
            'space_group': analyzer.get_space_group_symbol(),
            'space_group_number': analyzer.get_space_group_number(),
            'crystal_system': analyzer.get_crystal_system(),
            'point_group': analyzer.get_point_group_symbol(),
            'lattice_type': analyzer.get_lattice_type(),
            'is_ordered': analyzer.is_ordered()
        }
        
        return results
    
    def calculate_density(self):
        """Calculate film density and porosity"""
        volume = self.atoms.get_volume()
        mass = sum(self.atoms.get_masses())
        density = mass / volume  # g/cm³ (assuming amu/Å³)
        
        # Estimate theoretical density for comparison
        composition = self.get_composition()
        theo_density = self.estimate_theoretical_density(composition)
        porosity = 1 - (density / theo_density)
        
        return {
            'density': density,
            'theoretical_density': theo_density,
            'porosity': porosity,
            'volume': volume
        }
    
    def analyze_layer_structure(self):
        """Analyze layer thickness and interface properties"""
        positions = self.atoms.get_positions()
        z_coords = positions[:, 2]
        
        # Find layer boundaries (assumes growth in z-direction)
        layer_boundaries = self.find_layer_boundaries(z_coords)
        layer_thicknesses = np.diff(layer_boundaries)
        
        return {
            'total_thickness': z_coords.max() - z_coords.min(),
            'num_layers': len(layer_boundaries) - 1,
            'layer_thicknesses': layer_thicknesses,
            'average_layer_thickness': np.mean(layer_thicknesses),
            'thickness_uniformity': np.std(layer_thicknesses)
        }
    
    def calculate_mechanical_properties(self):
        """Estimate elastic constants and mechanical properties"""
        from ase.calculators.emt import EMT
        from ase.optimize import BFGS
        
        # Simple estimate using bond analysis
        bond_analysis = self.analyze_bonding()
        
        # Estimate bulk modulus from bond strengths
        avg_bond_strength = np.mean(bond_analysis['bond_strengths'])
        estimated_bulk_modulus = avg_bond_strength * 10  # Rough scaling
        
        return {
            'estimated_bulk_modulus': estimated_bulk_modulus,
            'bond_density': bond_analysis['bond_density'],
            'coordination_numbers': bond_analysis['coordination_numbers']
        }
    
    def analyze_bonding(self):
        """Analyze chemical bonding in the structure"""
        nl = NeighborList(
            [1.5] * len(self.atoms),  # Cutoff radii
            self_interaction=False,
            bothways=True
        )
        nl.update(self.atoms)
        
        bonds = []
        bond_strengths = []
        coordination_numbers = []
        
        for i in range(len(self.atoms)):
            indices, offsets = nl.get_neighbors(i)
            coordination_numbers.append(len(indices))
            
            for j, offset in zip(indices, offsets):
                if i < j:  # Avoid double counting
                    distance = self.atoms.get_distance(i, j, mic=True)
                    bond_strength = 1.0 / distance**2  # Simple estimate
                    bonds.append((i, j, distance))
                    bond_strengths.append(bond_strength)
        
        return {
            'bonds': bonds,
            'bond_strengths': bond_strengths,
            'bond_density': len(bonds) / self.atoms.get_volume(),
            'coordination_numbers': coordination_numbers,
            'average_coordination': np.mean(coordination_numbers)
        }
    
    def calculate_xrd_pattern(self):
        """Calculate theoretical XRD pattern"""
        calculator = XRDCalculator(wavelength="CuKa")
        pattern = calculator.get_pattern(self.structure)
        
        return {
            'two_theta': pattern.x,
            'intensities': pattern.y,
            'hkl_list': [hkl['hkl'] for hkl in pattern.hkls]
        }
    
    def find_layer_boundaries(self, z_coords, threshold=2.0):
        """Find layer boundaries based on density gaps"""
        hist, bin_edges = np.histogram(z_coords, bins=50)
        
        # Find gaps in density
        gaps = []
        for i in range(len(hist)-1):
            if hist[i] < threshold and hist[i+1] < threshold:
                gaps.append(bin_edges[i])
        
        return np.array([z_coords.min()] + gaps + [z_coords.max()])
    
    def get_composition(self):
        """Get elemental composition"""
        symbols = self.atoms.get_chemical_symbols()
        composition = {}
        for symbol in symbols:
            composition[symbol] = composition.get(symbol, 0) + 1
        return composition
    
    def estimate_theoretical_density(self, composition):
        """Estimate theoretical density based on composition"""
        # Atomic densities (rough estimates in g/cm³)
        atomic_densities = {
            'Al': 2.70, 'Si': 2.33, 'O': 1.43, 'C': 2.27, 'H': 0.09
        }
        
        total_mass = 0
        total_atoms = sum(composition.values())
        
        for element, count in composition.items():
            if element in atomic_densities:
                total_mass += count * atomic_densities[element]
        
        return total_mass / total_atoms
    
    def generate_report(self, output_file="bulk_analysis_report.txt"):
        """Generate comprehensive analysis report"""
        with open(output_file, 'w') as f:
            f.write("=== MLD Bulk Structure Analysis Report ===\n\n")
            
            # Crystallinity
            crystal_info = self.analyze_crystallinity()
            f.write("CRYSTALLINITY ANALYSIS:\n")
            for key, value in crystal_info.items():
                f.write(f"  {key}: {value}\n")
            f.write("\n")
            
            # Density
            density_info = self.calculate_density()
            f.write("DENSITY ANALYSIS:\n")
            for key, value in density_info.items():
                f.write(f"  {key}: {value:.3f}\n")
            f.write("\n")
            
            # Layer structure
            layer_info = self.analyze_layer_structure()
            f.write("LAYER STRUCTURE:\n")
            for key, value in layer_info.items():
                if isinstance(value, np.ndarray):
                    f.write(f"  {key}: {value}\n")
                else:
                    f.write(f"  {key}: {value:.3f}\n")
            f.write("\n")
            
            # Mechanical properties
            mech_props = self.calculate_mechanical_properties()
            f.write("MECHANICAL PROPERTIES:\n")
            for key, value in mech_props.items():
                if isinstance(value, (list, np.ndarray)):
                    f.write(f"  {key}: {np.mean(value):.3f} ± {np.std(value):.3f}\n")
                else:
                    f.write(f"  {key}: {value:.3f}\n")
        
        print(f"Analysis report saved to {output_file}")

def main():
    """Example usage"""
    # This would be called with actual structure files
    print("Bulk MLD structure analysis tools loaded.")
    print("Usage: analyzer = BulkMLDAnalyzer('structure.cif')")
    print("       report = analyzer.generate_report()")

if __name__ == "__main__":
    main()