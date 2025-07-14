#!/usr/bin/env python3
"""
Step 3: Smart H2O optimization using the intelligent hierarchy
Much simpler than the complex diol molecule!
"""

from smart_chemistry import optimize_molecule

def main():
    print("💧 Step 3: Smart H2O Optimization")
    print("=" * 50)
    
    # Try smart optimization for water
    atoms, log = optimize_molecule(
        molecule_name="water", 
        structure_file="structures/h2o_molecule.xyz"
    )
    
    if atoms:
        print("\n🎉 H2O optimization completed successfully!")
        print("\nOptimization pathway:")
        for step in log[-5:]:  # Show last 5 steps
            print(f"   {step}")
        
        # Show final structure info
        print(f"\n📊 Final H2O structure:")
        print(f"   Atoms: {len(atoms)}")
        print(f"   Energy: {atoms.get_potential_energy():.3f} eV")
        
    else:
        print("\n❌ H2O optimization failed!")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())