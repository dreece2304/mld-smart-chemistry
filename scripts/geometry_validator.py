#!/usr/bin/env python3
"""
Geometry Validator for DFT MLD Calculations
Prevents "atoms too close" errors and validates structures before expensive DFT

Essential for preventing the failure mode we encountered:
"Atoms are too close, e.g. 0.0 Å"
"""

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones
import os
from typing import Tuple, List, Dict, Optional

class GeometryValidator:
    """Validate and fix molecular geometries before DFT"""
    
    def __init__(self):
        """Initialize with literature-based validation criteria"""
        
        # Minimum distances based on covalent radii (Å)
        self.min_distances = {
            ('H', 'H'): 1.0,    # H-H minimum
            ('H', 'C'): 1.0,    # H-C minimum  
            ('H', 'O'): 0.8,    # H-O minimum
            ('H', 'Si'): 1.2,   # H-Si minimum
            ('H', 'Al'): 1.2,   # H-Al minimum
            ('C', 'C'): 1.3,    # C-C minimum
            ('C', 'O'): 1.2,    # C-O minimum
            ('C', 'Si'): 1.7,   # C-Si minimum
            ('C', 'Al'): 1.8,   # C-Al minimum
            ('O', 'O'): 1.2,    # O-O minimum
            ('O', 'Si'): 1.4,   # O-Si minimum
            ('O', 'Al'): 1.6,   # O-Al minimum
            ('Si', 'Si'): 2.0,  # Si-Si minimum
            ('Si', 'Al'): 2.2,  # Si-Al minimum
            ('Al', 'Al'): 2.4,  # Al-Al minimum
        }
        
        # Maximum reasonable distances (Å) - beyond this, atoms aren't interacting
        self.max_distances = {
            ('H', 'H'): 4.0,
            ('H', 'C'): 3.0,
            ('H', 'O'): 2.5,
            ('H', 'Si'): 3.0,
            ('H', 'Al'): 3.0,
            ('C', 'C'): 4.0,
            ('C', 'O'): 3.5,
            ('C', 'Si'): 4.0,
            ('C', 'Al'): 4.0,
            ('O', 'O'): 4.0,
            ('O', 'Si'): 3.5,
            ('O', 'Al'): 3.5,
            ('Si', 'Si'): 5.0,
            ('Si', 'Al'): 5.0,
            ('Al', 'Al'): 5.0,
        }
        
        # Expected bond lengths for validation (Å)
        self.expected_bonds = {
            ('H', 'C'): 1.09,   # C-H in methyl
            ('H', 'O'): 0.96,   # O-H hydroxyl
            ('C', 'Al'): 1.96,  # Al-C in TMA
            ('O', 'Si'): 1.65,  # Si-O surface
            ('O', 'Al'): 1.80,  # Al-O bond
            ('Si', 'Si'): 2.35, # Si-Si bulk
        }
        
        print("🔍 Geometry Validator for DFT MLD")
        print("   Prevents 'atoms too close' DFT failures")
        print("   Validates bond lengths and geometries")
    
    def get_min_distance(self, sym1: str, sym2: str) -> float:
        """Get minimum allowed distance between two elements"""
        key = tuple(sorted([sym1, sym2]))
        return self.min_distances.get(key, 1.0)  # Default 1.0 Å
    
    def get_max_distance(self, sym1: str, sym2: str) -> float:
        """Get maximum reasonable distance between two elements"""
        key = tuple(sorted([sym1, sym2]))
        return self.max_distances.get(key, 5.0)  # Default 5.0 Å
    
    def check_minimum_distances(self, atoms: Atoms) -> Tuple[bool, List[str]]:
        """
        Check if any atoms are too close together
        
        Returns:
            (is_valid, list_of_issues)
        """
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        issues = []
        
        n_atoms = len(atoms)
        min_distance_found = float('inf')
        closest_pair = None
        
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                sym1, sym2 = symbols[i], symbols[j]
                pos1, pos2 = positions[i], positions[j]
                
                distance = np.linalg.norm(pos1 - pos2)
                min_allowed = self.get_min_distance(sym1, sym2)
                
                if distance < min_distance_found:
                    min_distance_found = distance
                    closest_pair = (i, j, sym1, sym2, distance)
                
                if distance < min_allowed:
                    issues.append(
                        f"Atoms {i}({sym1}) and {j}({sym2}) too close: "
                        f"{distance:.3f} Å < {min_allowed:.3f} Å minimum"
                    )
        
        print(f"\n🔍 Distance Check Results:")
        print(f"   Closest pair: {closest_pair[2]}-{closest_pair[3]} = {closest_pair[4]:.3f} Å")
        
        if issues:
            print(f"   ❌ {len(issues)} distance violations found")
            for issue in issues[:5]:  # Show first 5
                print(f"      {issue}")
            if len(issues) > 5:
                print(f"      ... and {len(issues) - 5} more")
        else:
            print(f"   ✅ All distances OK (minimum: {min_distance_found:.3f} Å)")
        
        return len(issues) == 0, issues
    
    def check_bond_lengths(self, atoms: Atoms) -> Tuple[bool, List[str]]:
        """
        Check if bond lengths are reasonable
        
        Returns:
            (is_valid, list_of_warnings)
        """
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        warnings = []
        
        n_atoms = len(atoms)
        bond_count = 0
        
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                sym1, sym2 = symbols[i], symbols[j]
                pos1, pos2 = positions[i], positions[j]
                
                distance = np.linalg.norm(pos1 - pos2)
                max_bond = self.get_max_distance(sym1, sym2)
                
                # Check if this is a reasonable bond
                if distance < max_bond:
                    bond_count += 1
                    
                    # Check against expected bond length if available
                    key = tuple(sorted([sym1, sym2]))
                    if key in self.expected_bonds:
                        expected = self.expected_bonds[key]
                        deviation = abs(distance - expected)
                        
                        if deviation > 0.3:  # More than 0.3 Å deviation
                            warnings.append(
                                f"Bond {sym1}-{sym2}: {distance:.3f} Å "
                                f"(expected ~{expected:.3f} Å, deviation: {deviation:.3f} Å)"
                            )
        
        print(f"\n🔗 Bond Length Analysis:")
        print(f"   Bonds found: {bond_count}")
        
        if warnings:
            print(f"   ⚠️  {len(warnings)} unusual bond lengths")
            for warning in warnings[:5]:
                print(f"      {warning}")
        else:
            print(f"   ✅ All bond lengths reasonable")
        
        return len(warnings) == 0, warnings
    
    def check_overlapping_atoms(self, atoms: Atoms, threshold: float = 0.1) -> Tuple[bool, List[str]]:
        """
        Check for exactly overlapping atoms (critical DFT failure mode)
        
        Args:
            threshold: Distance below which atoms are considered overlapping
            
        Returns:
            (is_valid, list_of_overlaps)
        """
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        overlaps = []
        
        n_atoms = len(atoms)
        
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                distance = np.linalg.norm(positions[i] - positions[j])
                
                if distance < threshold:
                    overlaps.append(
                        f"Atoms {i}({symbols[i]}) and {j}({symbols[j]}) "
                        f"overlapping: {distance:.4f} Å"
                    )
        
        print(f"\n💥 Overlap Check (< {threshold} Å):")
        if overlaps:
            print(f"   ❌ {len(overlaps)} overlapping atoms found!")
            for overlap in overlaps:
                print(f"      {overlap}")
        else:
            print(f"   ✅ No overlapping atoms")
        
        return len(overlaps) == 0, overlaps
    
    def fix_close_contacts(self, atoms: Atoms, min_sep: float = 1.0) -> Atoms:
        """
        Attempt to fix atoms that are too close by small movements
        
        Args:
            atoms: Structure to fix
            min_sep: Minimum separation to enforce
            
        Returns:
            Fixed structure
        """
        print(f"\n🔧 Attempting to fix close contacts...")
        
        atoms_fixed = atoms.copy()
        positions = atoms_fixed.get_positions()
        symbols = atoms_fixed.get_chemical_symbols()
        
        n_fixes = 0
        max_iterations = 10
        
        for iteration in range(max_iterations):
            has_close_contacts = False
            
            for i in range(len(atoms_fixed)):
                for j in range(i + 1, len(atoms_fixed)):
                    distance = np.linalg.norm(positions[i] - positions[j])
                    min_allowed = max(min_sep, self.get_min_distance(symbols[i], symbols[j]))
                    
                    if distance < min_allowed:
                        has_close_contacts = True
                        
                        # Move atoms apart along their connecting vector
                        vector = positions[j] - positions[i]
                        if np.linalg.norm(vector) > 1e-6:
                            vector = vector / np.linalg.norm(vector)
                        else:
                            # Random direction if atoms are exactly overlapping
                            vector = np.random.randn(3)
                            vector = vector / np.linalg.norm(vector)
                        
                        # Move both atoms away from each other
                        move_distance = (min_allowed - distance) / 2 + 0.1
                        positions[i] -= vector * move_distance
                        positions[j] += vector * move_distance
                        n_fixes += 1
            
            if not has_close_contacts:
                break
        
        atoms_fixed.set_positions(positions)
        
        print(f"   Fixes applied: {n_fixes}")
        print(f"   Iterations: {iteration + 1}")
        
        # Validate fix
        is_valid, _ = self.check_minimum_distances(atoms_fixed)
        if is_valid:
            print(f"   ✅ Close contacts resolved")
        else:
            print(f"   ⚠️  Some close contacts remain")
        
        return atoms_fixed
    
    def pre_optimize_structure(self, atoms: Atoms, steps: int = 100) -> Atoms:
        """
        Pre-optimize structure with simple force field to fix bad geometries
        
        Args:
            atoms: Structure to optimize
            steps: Maximum optimization steps
            
        Returns:
            Pre-optimized structure
        """
        print(f"\n⚙️  Pre-optimizing with simple force field...")
        
        atoms_opt = atoms.copy()
        
        try:
            # Try Lennard-Jones for simple pre-optimization
            atoms_opt.calc = LennardJones()
            
            # Simple optimization
            from ase.optimize import BFGS
            opt = BFGS(atoms_opt, logfile=None)
            
            initial_energy = atoms_opt.get_potential_energy()
            opt.run(fmax=0.1, steps=steps)
            final_energy = atoms_opt.get_potential_energy()
            
            print(f"   Initial energy: {initial_energy:.3f} eV")
            print(f"   Final energy: {final_energy:.3f} eV")
            print(f"   Energy change: {final_energy - initial_energy:.3f} eV")
            print(f"   ✅ Pre-optimization complete")
            
        except Exception as e:
            print(f"   ⚠️  Pre-optimization failed: {e}")
            print(f"   Returning original structure")
            return atoms
        
        return atoms_opt
    
    def validate_structure(self, atoms: Atoms, name: str = "Structure") -> Dict:
        """
        Comprehensive structure validation
        
        Returns:
            Dictionary with validation results
        """
        print(f"\n{'='*60}")
        print(f"🔍 VALIDATING: {name}")
        print(f"{'='*60}")
        print(f"Atoms: {len(atoms)}")
        print(f"Elements: {set(atoms.get_chemical_symbols())}")
        
        results = {
            'name': name,
            'n_atoms': len(atoms),
            'valid': True,
            'issues': [],
            'warnings': []
        }
        
        # 1. Check for overlapping atoms (critical)
        overlaps_ok, overlaps = self.check_overlapping_atoms(atoms)
        if not overlaps_ok:
            results['valid'] = False
            results['issues'].extend(overlaps)
        
        # 2. Check minimum distances
        distances_ok, distance_issues = self.check_minimum_distances(atoms)
        if not distances_ok:
            results['valid'] = False
            results['issues'].extend(distance_issues)
        
        # 3. Check bond lengths (warnings only)
        bonds_ok, bond_warnings = self.check_bond_lengths(atoms)
        results['warnings'].extend(bond_warnings)
        
        # Summary
        print(f"\n📊 Validation Summary:")
        if results['valid']:
            print(f"   ✅ STRUCTURE IS VALID for DFT")
        else:
            print(f"   ❌ STRUCTURE HAS ISSUES")
            print(f"   Critical issues: {len(results['issues'])}")
        
        if results['warnings']:
            print(f"   ⚠️  Warnings: {len(results['warnings'])}")
        
        return results
    
    def fix_structure(self, atoms: Atoms, name: str = "Structure") -> Tuple[Atoms, bool]:
        """
        Attempt to fix structural issues
        
        Returns:
            (fixed_atoms, success)
        """
        print(f"\n🔧 FIXING STRUCTURE: {name}")
        
        # Start with original
        fixed_atoms = atoms.copy()
        
        # Step 1: Fix close contacts
        validation = self.validate_structure(fixed_atoms, f"{name} (original)")
        
        if not validation['valid']:
            print(f"\n🔧 Applying fixes...")
            
            # Fix overlapping/close atoms
            fixed_atoms = self.fix_close_contacts(fixed_atoms)
            
            # Pre-optimize with force field
            fixed_atoms = self.pre_optimize_structure(fixed_atoms)
            
            # Re-validate
            final_validation = self.validate_structure(fixed_atoms, f"{name} (fixed)")
            
            if final_validation['valid']:
                print(f"\n✅ Structure successfully fixed!")
                return fixed_atoms, True
            else:
                print(f"\n⚠️  Structure partially fixed but issues remain")
                return fixed_atoms, False
        else:
            print(f"\n✅ Structure was already valid")
            return fixed_atoms, True

def main():
    """Validate geometries for DFT calculations"""
    
    validator = GeometryValidator()
    
    print("\n" + "="*60)
    print("GEOMETRY VALIDATION FOR DFT MLD")
    print("="*60)
    
    # Look for structure files to validate
    structure_files = [
        'small_si100_2x2_4L_OH.xyz',
        'small_si100_2x2_6L_OH.xyz', 
        'trimethylaluminum_optimized.xyz',
        'water_optimized.xyz',
        'mld_final_structure.xyz'
    ]
    
    validated_files = []
    
    for filename in structure_files:
        if os.path.exists(filename):
            print(f"\n{'='*40}")
            print(f"Validating: {filename}")
            print("="*40)
            
            try:
                atoms = read(filename)
                validation = validator.validate_structure(atoms, filename)
                
                if validation['valid']:
                    print(f"✅ {filename} is ready for DFT")
                    validated_files.append((filename, 'VALID'))
                else:
                    print(f"❌ {filename} has issues")
                    
                    # Attempt to fix
                    fixed_atoms, success = validator.fix_structure(atoms, filename)
                    
                    if success:
                        fixed_filename = filename.replace('.xyz', '_fixed.xyz')
                        write(fixed_filename, fixed_atoms)
                        print(f"💾 Saved fixed structure: {fixed_filename}")
                        validated_files.append((fixed_filename, 'FIXED'))
                    else:
                        validated_files.append((filename, 'ISSUES'))
                        
            except Exception as e:
                print(f"❌ Error reading {filename}: {e}")
                validated_files.append((filename, 'ERROR'))
    
    # Summary
    print(f"\n{'='*60}")
    print("VALIDATION SUMMARY")
    print("="*60)
    
    if validated_files:
        for filename, status in validated_files:
            if status == 'VALID':
                print(f"✅ {filename} - Ready for DFT")
            elif status == 'FIXED':
                print(f"🔧 {filename} - Fixed and ready")
            elif status == 'ISSUES':
                print(f"⚠️  {filename} - Has remaining issues")
            else:
                print(f"❌ {filename} - Could not process")
    else:
        print("No structure files found to validate")
        print("\nTo validate structures:")
        print("1. Run create_small_surfaces.py first")
        print("2. Then run this validator")
    
    print(f"\n💡 Next steps:")
    print(f"  • Use validated structures in DFT calculations")
    print(f"  • Avoid structures marked with issues")
    print(f"  • Start with smallest valid structures first")
    
    return 0

if __name__ == "__main__":
    exit(main())