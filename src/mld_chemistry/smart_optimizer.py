#!/usr/bin/env python3
"""
Smart Chemistry Optimization Framework
Hierarchy: Database → Classical → DFT
Use for ALL molecular optimization scripts
"""

import os
import time
import requests
from ase.io import read, write
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from gpaw import GPAW, PW
from tqdm import tqdm

class SmartOptimizer:
    """Smart molecular optimization with fallback hierarchy"""
    
    def __init__(self, molecule_name, verbose=True):
        self.molecule_name = molecule_name
        self.verbose = verbose
        self.optimization_log = []
    
    def log(self, message, level="INFO"):
        """Log optimization steps"""
        if self.verbose:
            icons = {"INFO": "ℹ️", "SUCCESS": "✅", "WARNING": "⚠️", "ERROR": "❌"}
            print(f"{icons.get(level, 'ℹ️')} {message}")
        self.optimization_log.append(f"{level}: {message}")
    
    def try_database_lookups(self):
        """Method 1: Try multiple databases in order of quality"""
        
        # Database priority: Best quality first
        databases = [
            ("NIST", self._try_nist_lookup),
            ("PubChem_3D", self._try_pubchem_3d_lookup), 
            ("PubChem", self._try_pubchem_lookup),
            ("ChEBI", self._try_chebi_lookup)
        ]
        
        for db_name, lookup_func in databases:
            self.log(f"Trying {db_name} lookup for {self.molecule_name}...")
            try:
                atoms, method = lookup_func()
                if atoms is not None:
                    self.log(f"Found {self.molecule_name} in {db_name}! ({len(atoms)} atoms)", "SUCCESS")
                    return atoms, f"{db_name.lower()}_db"
            except Exception as e:
                self.log(f"{db_name} lookup failed: {e}", "WARNING")
        
        return None, None
    
    def _try_nist_lookup(self):
        """Try NIST WebBook (highest quality thermochemical data)"""
        # NIST has limited API, but very high quality
        # For now, placeholder - could add NIST API calls
        return None, None
    
    def _try_pubchem_3d_lookup(self):
        """Try PubChem 3D conformer (better than 2D)"""
        try:
            # PubChem 3D conformer API
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{self.molecule_name}/conformers/SDF"
            response = requests.get(url, timeout=10)
            
            if response.status_code == 200:
                with open(f"{self.molecule_name}_pubchem3d.sdf", "w") as f:
                    f.write(response.text)
                
                atoms = read(f"{self.molecule_name}_pubchem3d.sdf")
                atoms.center(vacuum=8.0)
                os.remove(f"{self.molecule_name}_pubchem3d.sdf")
                
                return atoms, "pubchem_3d"
                
        except Exception:
            pass
        return None, None
    
    def _try_pubchem_lookup(self):
        """Try standard PubChem lookup"""
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{self.molecule_name}/SDF"
            response = requests.get(url, timeout=10)
            
            if response.status_code == 200:
                with open(f"{self.molecule_name}_pubchem.sdf", "w") as f:
                    f.write(response.text)
                
                atoms = read(f"{self.molecule_name}_pubchem.sdf")
                atoms.center(vacuum=8.0)
                os.remove(f"{self.molecule_name}_pubchem.sdf")
                
                return atoms, "pubchem"
                
        except Exception:
            pass
        return None, None
    
    def _try_chebi_lookup(self):
        """Try ChEBI database (high quality, smaller coverage)"""
        # ChEBI has very high quality but limited API
        # Placeholder for future implementation
        return None, None
    
    def try_classical_optimization(self, atoms):
        """Method 2: Classical force field optimization"""
        self.log("Running classical optimization (EMT force field)...")
        
        try:
            # Use EMT (fast classical calculator)
            atoms_copy = atoms.copy()
            atoms_copy.calc = EMT()
            
            # Quick classical optimization
            opt = BFGS(atoms_copy, trajectory=f"{self.molecule_name}_classical.traj")
            
            start_time = time.time()
            opt.run(fmax=0.1, steps=50)  # Loose convergence, few steps
            elapsed = time.time() - start_time
            
            energy = atoms_copy.get_potential_energy()
            
            self.log(f"Classical optimization complete! E={energy:.3f} eV, t={elapsed:.1f}s", "SUCCESS")
            return atoms_copy, "classical"
            
        except Exception as e:
            self.log(f"Classical optimization failed: {e}", "WARNING")
            return atoms, "initial"
    
    def dft_optimization(self, atoms, max_steps=100):
        """Method 3: Full DFT optimization (fallback)"""
        self.log("Running DFT optimization (GPAW)...")
        
        # DFT calculator
        calc = GPAW(
            mode=PW(300),
            xc='PBE',
            txt=f'{self.molecule_name}_dft.log'
        )
        
        atoms.calc = calc
        
        # DFT optimization with progress bar
        class ProgressBFGS(BFGS):
            def __init__(self, atoms, trajectory=None, max_steps=100, name="DFT"):
                super().__init__(atoms, trajectory=trajectory)
                self.pbar = tqdm(total=max_steps, desc=f"🔬 {name} Optimization", 
                               unit="step", colour="green")
                self.step_count = 0
            
            def step(self, forces=None):
                result = super().step(forces)
                self.step_count += 1
                
                # Get forces directly from atoms if not provided
                try:
                    if forces is None:
                        forces = self.atoms.get_forces()
                    
                    fmax = (forces**2).sum(axis=1).max()**0.5  # Maximum force magnitude
                    # Show convergence progress as percentage  
                    convergence_pct = max(0, min(100, (0.5 - fmax) / 0.5 * 100))
                    self.pbar.set_postfix({
                        'Step': self.step_count, 
                        'Max_Force': f'{fmax:.3f}',
                        'Target': '0.050',
                        'Conv%': f'{convergence_pct:.0f}%'
                    })
                except:
                    # Fallback if force calculation fails
                    self.pbar.set_postfix({
                        'Step': self.step_count,
                        'Status': 'calculating...'
                    })
                
                self.pbar.update(1)
                return result
            
            def converged(self, forces=None):
                converged = super().converged(forces)
                if converged:
                    if forces is not None:
                        final_fmax = max(abs(forces.flatten()))
                        self.pbar.set_description(f"✅ Converged (fmax={final_fmax:.3f})")
                    else:
                        self.pbar.set_description("✅ DFT Converged")
                    self.pbar.close()
                return converged
        
        # Initial SCF
        self.log("Running initial SCF calculation...")
        initial_energy = atoms.get_potential_energy()
        self.log(f"Initial energy: {initial_energy:.3f} eV")
        
        # Optimization
        opt = ProgressBFGS(atoms, trajectory=f"{self.molecule_name}_dft.traj", 
                          max_steps=max_steps, name=self.molecule_name.upper())
        
        start_time = time.time()
        opt.run(fmax=0.05)
        elapsed = time.time() - start_time
        
        final_energy = atoms.get_potential_energy()
        
        self.log(f"DFT optimization complete! E={final_energy:.3f} eV, t={elapsed:.1f}s", "SUCCESS")
        return atoms, "dft"
    
    def smart_optimize(self, fallback_structure_file=None):
        """Smart optimization with full hierarchy"""
        
        self.log(f"🧠 Starting smart optimization for {self.molecule_name}")
        self.log("Strategy: Database → Classical → DFT")
        
        # Method 1: Try multiple databases
        atoms, method = self.try_database_lookups()
        
        # Method 2: If no database hit, start from file + classical
        if atoms is None and fallback_structure_file:
            try:
                self.log(f"Loading from file: {fallback_structure_file}")
                atoms = read(fallback_structure_file)
                atoms.center(vacuum=8.0)
                atoms, method = self.try_classical_optimization(atoms)
            except Exception as e:
                self.log(f"File loading failed: {e}", "ERROR")
                return None, None
        
        # Method 3: DFT refinement (always do this for accuracy)
        if atoms is not None:
            # Only do limited DFT if we had a good starting point
            max_dft_steps = 20 if method in ["pubchem", "classical"] else 100
            atoms_final, final_method = self.dft_optimization(atoms, max_dft_steps)
            
            # Save final result
            output_file = f"{self.molecule_name}_optimized.xyz"
            write(output_file, atoms_final)
            
            self.log(f"Final result saved to: {output_file}", "SUCCESS")
            self.log(f"Optimization path: {method} → {final_method}")
            
            return atoms_final, self.optimization_log
        
        self.log("All optimization methods failed!", "ERROR")
        return None, self.optimization_log

# Convenience functions for common molecules
def optimize_molecule(molecule_name, structure_file=None):
    """One-line smart optimization"""
    optimizer = SmartOptimizer(molecule_name)
    return optimizer.smart_optimize(structure_file)

def quick_test():
    """Test the framework"""
    print("🧪 Testing Smart Chemistry Framework")
    
    # Test with water (should find in PubChem)
    atoms, log = optimize_molecule("water")
    if atoms:
        print("✅ Water optimization successful!")
    else:
        print("❌ Water optimization failed!")

if __name__ == "__main__":
    quick_test()