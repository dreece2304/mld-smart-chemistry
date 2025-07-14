#!/usr/bin/env python3
"""
Step 2: Smart TMA optimization using the intelligent hierarchy
Database → Classical → DFT
"""

from smart_chemistry import optimize_molecule

def main():
    print("🧠 Step 2: Smart TMA Optimization")
    print("=" * 50)
    
    # Try smart optimization (will attempt database → classical → DFT)
    atoms, log = optimize_molecule(
        molecule_name="trimethylaluminum", 
        structure_file="structures/tma_molecule.xyz"
    )
    
    if atoms:
        print("\n🎉 TMA optimization completed successfully!")
        print("\nOptimization pathway:")
        for step in log[-5:]:  # Show last 5 steps
            print(f"   {step}")
    else:
        print("\n❌ TMA optimization failed!")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())