#!/usr/bin/env python3
"""
Step 2: Simple TMA optimization with accurate progress bar
"""

from ase.io import read, write
from ase.optimize import BFGS
from gpaw import GPAW, PW
from tqdm import tqdm
import time

class ProgressBFGS(BFGS):
    """BFGS with progress bar"""
    
    def __init__(self, atoms, trajectory=None, max_steps=100):
        super().__init__(atoms, trajectory=trajectory)
        self.max_steps = max_steps
        self.pbar = tqdm(total=max_steps, desc="🔧 TMA Optimization", 
                        unit="step", colour="blue")
        self.step_count = 0
    
    def step(self, forces=None):
        result = super().step(forces)
        self.step_count += 1
        
        # Update progress bar
        fmax = max(abs(forces.flatten())) if forces is not None else 0
        self.pbar.set_postfix({
            'Step': self.step_count,
            'Max_Force': f'{fmax:.3f}',
            'Target': '0.050'
        })
        self.pbar.update(1)
        
        return result
    
    def converged(self, forces=None):
        converged = super().converged(forces)
        if converged:
            self.pbar.set_description("✅ TMA Converged")
            self.pbar.close()
        return converged

def optimize_tma():
    print("🔧 Step 2: Optimizing TMA molecule")
    
    # Show threading info
    import os
    from gpaw.mpi import world
    print(f"   🖥️  Available CPU cores: {os.cpu_count()}")
    print(f"   🔄 GPAW MPI processes: {world.size}")
    
    # Read TMA
    tma = read('structures/tma_molecule.xyz')
    tma.center(vacuum=8.0)
    
    print(f"   🧬 TMA atoms: {len(tma)}")
    
    # Simple calculator with parallel settings
    calc = GPAW(
        mode=PW(300),
        xc='PBE',
        txt='tma_opt.log',
        parallel={'domain': min(world.size, 2)}  # Use domain decomposition
    )
    
    tma.calc = calc
    
    # Initial SCF can take time - let's do initial energy calculation with feedback
    print("🔄 Running initial SCF calculation (this may take 2-5 minutes)...")
    initial_energy = tma.get_potential_energy()
    print(f"✅ Initial SCF complete! Energy: {initial_energy:.3f} eV")
    print("🚀 Starting optimization steps...")
    
    # Optimize with progress bar
    opt = ProgressBFGS(tma, trajectory='tma_opt.traj', max_steps=100)
    
    start_time = time.time()
    opt.run(fmax=0.05)  # Converge forces to 0.05 eV/Å
    elapsed = time.time() - start_time
    
    # Get final energy
    energy = tma.get_potential_energy()
    
    # Save result
    write('tma_optimized.xyz', tma)
    
    print(f"✅ TMA optimization complete")
    print(f"   Energy: {energy:.3f} eV")
    print(f"   Time: {elapsed:.1f} seconds")
    print(f"   Saved: tma_optimized.xyz")

if __name__ == "__main__":
    optimize_tma()