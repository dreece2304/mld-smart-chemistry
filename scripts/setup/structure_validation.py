#!/usr/bin/env python3
"""
Structure file validation utilities for MLD simulations
Validates molecular geometries and catches common structural issues
"""

import os
import numpy as np
from ase.io import read, write
from ase import Atoms
from scipy.spatial.distance import pdist
import argparse

class StructureValidator:
    """Validates molecular structure files for MLD simulations"""
    
    def __init__(self, verbose=False):
        self.verbose = verbose
        self.issues = []
        self.warnings = []
        self.fixes_applied = []
    
    def print_status(self, message, status='info', indent=0):
        """Print status message"""
        colors = {
            'pass': '\033[92m✅',
            'fail': '\033[91m❌', 
            'warn': '\033[93m⚠️ ',
            'info': '\033[94mℹ️ ',
            'fix': '\033[95m🔧'
        }
        
        indent_str = "  " * indent
        color = colors.get(status, '')
        end_color = '\033[0m' if color else ''
        
        print(f"{indent_str}{color} {message}{end_color}")
    
    def validate_xyz_format(self, file_path):
        """Validate XYZ file format"""
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
            
            if len(lines) < 2:
                self.issues.append(f"{file_path}: Too few lines")
                return False
            
            # Check first line (number of atoms)
            try:
                n_atoms = int(lines[0].strip())
            except ValueError:
                self.issues.append(f"{file_path}: First line must be number of atoms")
                return False
            
            # Check if we have enough lines
            expected_lines = n_atoms + 2  # header + comment + atoms
            if len(lines) < expected_lines:
                self.issues.append(f"{file_path}: Not enough lines ({len(lines)} < {expected_lines})")
                return False
            
            # Validate atom lines
            for i, line in enumerate(lines[2:2+n_atoms], 2):
                parts = line.strip().split()
                if len(parts) < 4:
                    self.issues.append(f"{file_path}: Line {i+1} has too few columns")
                    return False
                
                # Check if coordinates are numbers
                try:
                    float(parts[1])  # x
                    float(parts[2])  # y  
                    float(parts[3])  # z
                except ValueError:
                    self.issues.append(f"{file_path}: Line {i+1} has invalid coordinates")
                    return False
            
            return True
            
        except Exception as e:
            self.issues.append(f"{file_path}: Could not read file - {e}")
            return False
    
    def validate_molecular_geometry(self, atoms, molecule_name):
        """Validate molecular geometry for physical reasonableness"""
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        
        # Calculate all pairwise distances
        distances = pdist(positions)
        
        # Check for atoms too close (< 0.5 Å)
        min_dist = np.min(distances)
        if min_dist < 0.5:
            self.issues.append(f"{molecule_name}: Atoms too close ({min_dist:.3f} Å)")
            return False
        
        # Check for isolated atoms (> 10 Å from nearest neighbor)
        n_atoms = len(atoms)
        isolated_atoms = []
        
        for i in range(n_atoms):
            min_neighbor_dist = float('inf')
            for j in range(n_atoms):
                if i != j:
                    dist = np.linalg.norm(positions[i] - positions[j])
                    if dist < min_neighbor_dist:
                        min_neighbor_dist = dist
            
            if min_neighbor_dist > 10.0:
                isolated_atoms.append(i)
        
        if isolated_atoms:
            self.warnings.append(f"{molecule_name}: Potentially isolated atoms: {isolated_atoms}")
        
        return True
    
    def validate_tma_structure(self, atoms):
        """Validate TMA (Al(CH3)3) specific structure"""
        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()
        
        # Check composition
        al_count = symbols.count('Al')
        c_count = symbols.count('C')
        h_count = symbols.count('H')
        
        if al_count != 1:
            self.issues.append(f"TMA: Expected 1 Al atom, found {al_count}")
            return False
        
        if c_count != 3:
            self.issues.append(f"TMA: Expected 3 C atoms, found {c_count}")
            return False
        
        if h_count != 9:
            self.issues.append(f"TMA: Expected 9 H atoms, found {h_count}")
            return False
        
        # Find Al atom
        al_idx = symbols.index('Al')
        al_pos = positions[al_idx]
        
        # Check Al-C distances (should be ~2.0 Å)
        c_indices = [i for i, s in enumerate(symbols) if s == 'C']
        al_c_distances = []
        
        for c_idx in c_indices:
            dist = np.linalg.norm(positions[c_idx] - al_pos)
            al_c_distances.append(dist)
        
        # Typical Al-C distance in TMA is ~1.97 Å
        for i, dist in enumerate(al_c_distances):
            if not (1.5 < dist < 2.5):
                self.warnings.append(f"TMA: Al-C{i+1} distance {dist:.3f} Å unusual (expect ~1.97 Å)")
        
        self.print_status(f"TMA Al-C distances: {[f'{d:.3f}' for d in al_c_distances]} Å", 'info', 1)
        
        return True
    
    def validate_diol_structure(self, atoms):
        """Validate 2-butyne-1,4-diol structure"""
        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()
        
        # Check composition: C4H6O2
        c_count = symbols.count('C')
        h_count = symbols.count('H')
        o_count = symbols.count('O')
        
        if c_count != 4:
            self.issues.append(f"Diol: Expected 4 C atoms, found {c_count}")
            return False
        
        if o_count != 2:
            self.issues.append(f"Diol: Expected 2 O atoms, found {o_count}")
            return False
        
        if h_count != 6:
            self.issues.append(f"Diol: Expected 6 H atoms, found {h_count}")
            return False
        
        # Check for reasonable C≡C triple bond distance
        c_indices = [i for i, s in enumerate(symbols) if s == 'C']
        c_positions = positions[c_indices]
        
        # Find shortest C-C distance (should be triple bond ~1.21 Å)
        c_c_distances = []
        for i in range(len(c_indices)):
            for j in range(i+1, len(c_indices)):
                dist = np.linalg.norm(c_positions[i] - c_positions[j])
                c_c_distances.append(dist)
        
        triple_bond_dist = min(c_c_distances)
        if not (1.10 < triple_bond_dist < 1.35):
            self.warnings.append(f"Diol: Shortest C-C distance {triple_bond_dist:.3f} Å (expect ~1.21 Å for C≡C)")
        
        self.print_status(f"Diol shortest C-C distance: {triple_bond_dist:.3f} Å", 'info', 1)
        
        return True
    
    def validate_surface_structure(self, atoms):
        """Validate Si surface structure"""
        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()
        
        # Check that we have silicon
        si_count = symbols.count('Si')
        if si_count == 0:
            self.issues.append("Surface: No Si atoms found")
            return False
        
        # Check for reasonable Si-Si distances
        si_indices = [i for i, s in enumerate(symbols) if s == 'Si']
        si_positions = positions[si_indices]
        
        if len(si_indices) > 1:
            si_si_distances = []
            for i in range(len(si_indices)):
                for j in range(i+1, len(si_indices)):
                    dist = np.linalg.norm(si_positions[i] - si_positions[j])
                    si_si_distances.append(dist)
            
            min_si_si = min(si_si_distances)
            max_si_si = max(si_si_distances)
            
            # Si-Si distance in bulk Si is ~2.35 Å
            if min_si_si < 1.8:
                self.warnings.append(f"Surface: Very short Si-Si distance {min_si_si:.3f} Å")
            
            self.print_status(f"Surface Si-Si distances: {min_si_si:.3f} - {max_si_si:.3f} Å", 'info', 1)
        
        # Check for hydroxyl groups if present
        h_count = symbols.count('H')
        o_count = symbols.count('O')
        
        if o_count > 0:
            self.print_status(f"Surface hydroxylation: {o_count} O atoms, {h_count} H atoms", 'info', 1)
        
        return True
    
    def fix_common_issues(self, atoms, molecule_name):
        """Fix common structural issues"""
        fixed = False
        
        # Add small random displacement to break symmetry
        positions = atoms.get_positions()
        if molecule_name.lower() in ['diol', 'butyne_diol']:
            # Add tiny random displacement
            np.random.seed(42)  # Reproducible
            noise = np.random.normal(0, 0.005, positions.shape)
            positions += noise
            atoms.set_positions(positions)
            fixed = True
            self.fixes_applied.append(f"{molecule_name}: Added symmetry-breaking displacement")
        
        # Center molecule with adequate vacuum
        original_cell = atoms.get_cell()
        atoms.center(vacuum=8.0)
        new_cell = atoms.get_cell()
        
        if not np.allclose(original_cell.diagonal(), new_cell.diagonal()):
            fixed = True
            self.fixes_applied.append(f"{molecule_name}: Centered with 8 Å vacuum")
        
        return fixed
    
    def validate_structure_file(self, file_path, fix_issues=False):
        """Validate a single structure file"""
        if not os.path.exists(file_path):
            self.issues.append(f"File not found: {file_path}")
            return False
        
        # Check file format
        if not self.validate_xyz_format(file_path):
            return False
        
        # Load structure
        try:
            atoms = read(file_path)
        except Exception as e:
            self.issues.append(f"{file_path}: Could not load structure - {e}")
            return False
        
        # Basic geometry validation
        molecule_name = os.path.basename(file_path).replace('.xyz', '')
        if not self.validate_molecular_geometry(atoms, molecule_name):
            return False
        
        # Molecule-specific validation
        success = True
        if 'tma' in molecule_name.lower():
            success = self.validate_tma_structure(atoms)
        elif 'diol' in molecule_name.lower() or 'butyne' in molecule_name.lower():
            success = self.validate_diol_structure(atoms)
        elif 'surface' in molecule_name.lower() or 'si' in molecule_name.lower():
            success = self.validate_surface_structure(atoms)
        
        # Apply fixes if requested
        if fix_issues and success:
            if self.fix_common_issues(atoms, molecule_name):
                # Save fixed structure
                fixed_file = file_path.replace('.xyz', '_fixed.xyz')
                write(fixed_file, atoms)
                self.print_status(f"Saved fixed structure: {fixed_file}", 'fix', 1)
        
        return success
    
    def validate_all_structures(self, fix_issues=False):
        """Validate all MLD structure files"""
        structure_files = [
            'structures/tma_molecule.xyz',
            'structures/butyne_diol_molecule.xyz',
            'structures/simple_si_surface.xyz',
            'structures/diol_improved.xyz'  # Optional
        ]
        
        self.print_status("Validating MLD Structure Files", 'info')
        
        validated_count = 0
        for file_path in structure_files:
            if os.path.exists(file_path):
                self.print_status(f"Validating {file_path}...", 'info', 1)
                if self.validate_structure_file(file_path, fix_issues):
                    self.print_status(f"{file_path}: Valid", 'pass', 1)
                    validated_count += 1
                else:
                    self.print_status(f"{file_path}: Issues found", 'fail', 1)
            else:
                if 'improved' not in file_path:  # Optional file
                    self.issues.append(f"Required file missing: {file_path}")
                    self.print_status(f"{file_path}: Missing", 'fail', 1)
                else:
                    self.print_status(f"{file_path}: Not found (optional)", 'warn', 1)
        
        return validated_count
    
    def generate_report(self):
        """Generate validation report"""
        print(f"\n{'='*60}")
        print("STRUCTURE VALIDATION REPORT")
        print(f"{'='*60}")
        
        if self.issues:
            print(f"\n❌ ISSUES ({len(self.issues)}):")
            for i, issue in enumerate(self.issues, 1):
                print(f"  {i}. {issue}")
        
        if self.warnings:
            print(f"\n⚠️  WARNINGS ({len(self.warnings)}):")
            for i, warning in enumerate(self.warnings, 1):
                print(f"  {i}. {warning}")
        
        if self.fixes_applied:
            print(f"\n🔧 FIXES APPLIED ({len(self.fixes_applied)}):")
            for i, fix in enumerate(self.fixes_applied, 1):
                print(f"  {i}. {fix}")
        
        # Overall status
        if not self.issues:
            print(f"\n✅ STRUCTURE VALIDATION PASSED")
        else:
            print(f"\n❌ STRUCTURE VALIDATION FAILED")
        
        return len(self.issues) == 0

def main():
    """Main structure validation function"""
    
    parser = argparse.ArgumentParser(
        description='Validate MLD structure files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python structure_validation.py                    # Validate all structures
  python structure_validation.py --fix-issues       # Validate and fix issues
  python structure_validation.py structures/tma_molecule.xyz  # Single file
        """
    )
    
    parser.add_argument('files', nargs='*', help='Specific files to validate')
    parser.add_argument('--fix-issues', action='store_true',
                       help='Attempt to fix common issues')
    parser.add_argument('--verbose', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    # Create validator
    validator = StructureValidator(verbose=args.verbose)
    
    print("🔍 MLD STRUCTURE VALIDATION")
    print("=" * 40)
    
    if args.files:
        # Validate specific files
        for file_path in args.files:
            validator.validate_structure_file(file_path, args.fix_issues)
    else:
        # Validate all structures
        validator.validate_all_structures(args.fix_issues)
    
    # Generate report
    success = validator.generate_report()
    
    return 0 if success else 1

if __name__ == "__main__":
    exit(main())