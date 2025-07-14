#!/usr/bin/env python3
"""
Step 5: Simple TMA + H2O gas phase reaction
Al(CH3)3 + H2O → Al(CH3)2OH + CH4
"""

from ase.io import read, write
from ase import Atoms
from smart_chemistry import SmartOptimizer
import numpy as np

def create_reaction_system():
    """Combine TMA and H2O for reaction simulation"""
    
    print("⚛️  Step 5: TMA + H2O Gas Phase Reaction")
    print("=" * 50)
    print("🧪 Reaction: Al(CH₃)₃ + H₂O → Al(CH₃)₂OH + CH₄")
    
    # Load optimized molecules
    try:
        tma = read('trimethylaluminum_optimized.xyz')
        h2o = read('water_optimized.xyz')
        print(f"✅ Loaded TMA ({len(tma)} atoms) and H2O ({len(h2o)} atoms)")
    except FileNotFoundError as e:
        print(f"❌ Missing optimized molecules: {e}")
        return None
    
    # Position molecules for reaction
    # Put H2O about 3 Å from TMA (close enough to react)
    print("📐 Positioning molecules for reaction...")
    
    # Center TMA at origin
    tma_pos = tma.get_positions()
    tma_center = tma_pos.mean(axis=0)
    tma.translate(-tma_center)
    
    # Position H2O near Al atom (find Al in TMA)
    al_idx = None
    for i, symbol in enumerate(tma.get_chemical_symbols()):
        if symbol == 'Al':
            al_idx = i
            break
    
    if al_idx is None:
        print("❌ No Al atom found in TMA!")
        return None
    
    al_pos = tma.get_positions()[al_idx]
    
    # Place H2O 3 Å away from Al
    h2o_offset = al_pos + [3.0, 0.0, 0.0]  # 3 Å in X direction
    h2o.translate(h2o_offset - h2o.get_positions().mean(axis=0))
    
    # Combine into reaction system
    reaction_positions = np.vstack([tma.get_positions(), h2o.get_positions()])
    reaction_symbols = list(tma.get_chemical_symbols()) + list(h2o.get_chemical_symbols())
    
    reaction_system = Atoms(
        symbols=reaction_symbols,
        positions=reaction_positions
    )
    
    # Add large vacuum for gas phase
    reaction_system.center(vacuum=12.0)
    
    print(f"✅ Created reaction system:")
    print(f"   Total atoms: {len(reaction_system)}")
    print(f"   Al-O distance: {np.linalg.norm(al_pos - h2o_offset):.2f} Å")
    
    # Save initial system
    write('tma_h2o_initial.xyz', reaction_system)
    print(f"💾 Saved initial system: tma_h2o_initial.xyz")
    
    return reaction_system

def optimize_reaction_system():
    """Optimize the TMA + H2O system to find reaction products"""
    
    reaction_system = create_reaction_system()
    if reaction_system is None:
        return False
    
    print(f"\n🔧 Optimizing reaction system...")
    
    # Use our smart optimizer for the reaction
    optimizer = SmartOptimizer("tma_h2o_reaction", verbose=True)
    
    # Skip database lookup (this is a custom system)
    # Go straight to DFT optimization
    try:
        # Use looser convergence for reaction system (more complex)
        atoms_opt, method = optimizer.dft_optimization(reaction_system, max_steps=50)
        
        # Save optimized reaction system
        write('tma_h2o_optimized.xyz', atoms_opt)
        
        print(f"\n✅ Reaction optimization complete!")
        print(f"💾 Saved: tma_h2o_optimized.xyz")
        
        # Basic analysis
        initial_energy = reaction_system.get_potential_energy() if hasattr(reaction_system, 'calc') else None
        final_energy = atoms_opt.get_potential_energy()
        
        print(f"\n📊 Reaction Analysis:")
        print(f"   Final energy: {final_energy:.3f} eV")
        if initial_energy:
            reaction_energy = final_energy - initial_energy
            print(f"   Reaction energy: {reaction_energy:.3f} eV")
        
        return True
        
    except Exception as e:
        print(f"❌ Reaction optimization failed: {e}")
        return False

def main():
    """Run TMA + H2O reaction simulation"""
    
    success = optimize_reaction_system()
    
    if success:
        print(f"\n🎉 TMA + H2O reaction simulation complete!")
        print(f"📊 Visualize with: python visualize_molecules.py")
        print(f"\nNext: Step 6 - TMA + Surface reaction")
    else:
        print(f"\n❌ Reaction simulation failed!")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())