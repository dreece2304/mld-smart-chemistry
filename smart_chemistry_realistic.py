#!/usr/bin/env python3
"""
Realistic Smart Chemistry Module for MLD Simulations
Improved version with proper force fields and error handling
"""

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS, LBFGS
from ase.constraints import FixAtoms
from ase.calculators.calculator import Calculator
from typing import Tuple, Optional, Dict, List
import requests
import os
import json
from datetime import datetime

# Try to import various calculators
try:
    from gpaw import GPAW, PW
    GPAW_AVAILABLE = True
except ImportError:
    GPAW_AVAILABLE = False
    print("Warning: GPAW not available, DFT calculations disabled")

try:
    from ase.calculators.lj import LennardJones
    LJ_AVAILABLE = True
except ImportError:
    LJ_AVAILABLE = False

try:
    from ase.calculators.uff import UFF
    UFF_AVAILABLE = True
except ImportError:
    UFF_AVAILABLE = False
    print("Info: UFF not available, using simpler force fields")

class RealisticOptimizer:
    """
    Realistic molecular optimizer with proper force fields and convergence criteria
    """
    
    def __init__(self, molecule_name: str, work_dir: str = ".", verbose: bool = True):
        """
        Initialize optimizer
        
        Args:
            molecule_name: Name of molecule for logging
            work_dir: Working directory for files
            verbose: Print detailed progress
        """
        self.molecule_name = molecule_name
        self.work_dir = work_dir
        self.verbose = verbose
        self.log_file = os.path.join(work_dir, f"{molecule_name}_optimization.log")
        
        # Optimization parameters (based on best practices)
        self.mm_fmax = 0.1      # eV/Å for molecular mechanics
        self.dft_fmax = 0.02    # eV/Å for DFT (tighter)
        self.mm_max_steps = 500
        self.dft_max_steps = 200
        
        # Initialize log
        self.log_entries = []
        self.log(f"Initialized RealisticOptimizer for {molecule_name}")
    
    def log(self, message: str, level: str = "INFO"):
        """Log messages with timestamp"""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_entry = f"[{timestamp}] [{level}] {message}"
        self.log_entries.append(log_entry)
        
        if self.verbose:
            # Color coding for terminal
            colors = {
                "INFO": "\033[0m",      # Normal
                "SUCCESS": "\033[92m",   # Green
                "WARNING": "\033[93m",   # Yellow
                "ERROR": "\033[91m"      # Red
            }
            color = colors.get(level, "\033[0m")
            print(f"{color}{log_entry}\033[0m")
    
    def save_log(self):
        """Save log to file"""
        with open(self.log_file, 'w') as f:
            f.write("\n".join(self.log_entries))
    
    def get_molecular_ff_calculator(self) -> Optional[Calculator]:
        """
        Get appropriate molecular mechanics calculator
        
        Returns:
            Calculator object or None
        """
        if UFF_AVAILABLE:
            self.log("Using UFF (Universal Force Field) calculator")
            return UFF()
        elif LJ_AVAILABLE:
            self.log("Using Lennard-Jones calculator (less accurate)", "WARNING")
            return LennardJones()
        else:
            self.log("No molecular force field available!", "ERROR")
            return None
    
    def get_dft_calculator(self, atoms: Atoms) -> Optional[Calculator]:
        """
        Get DFT calculator with appropriate settings
        
        Args:
            atoms: Structure to calculate
            
        Returns:
            GPAW calculator or None
        """
        if not GPAW_AVAILABLE:
            self.log("GPAW not available for DFT calculations", "ERROR")
            return None
        
        # Calculate appropriate cell size
        positions = atoms.get_positions()
        if positions.size > 0:
            pmin = positions.min(axis=0)
            pmax = positions.max(axis=0)
            size = pmax - pmin
            
            # Ensure adequate vacuum (at least 8 Å, but scale with molecule size)
            vacuum = max(8.0, 0.5 * size.max())
            atoms.center(vacuum=vacuum)
            
            self.log(f"Molecule size: {size} Å, using vacuum: {vacuum:.1f} Å")
        
        # DFT parameters
        try:
            calc = GPAW(
                mode=PW(400),           # Higher cutoff for better accuracy
                xc='PBE',               # Standard functional
                txt=f'{self.molecule_name}_gpaw.txt',
                symmetry='off',         # For molecules
                nbands='nao',           # Enough bands
                convergence={
                    'energy': 1e-5,     # eV
                    'density': 1e-4,
                    'eigenstates': 1e-4,
                },
                mixer={'beta': 0.25}    # Slower mixing for stability
            )
            self.log("Created GPAW calculator with PBE/PW(400)")
            return calc
        except Exception as e:
            self.log(f"Failed to create GPAW calculator: {e}", "ERROR")
            return None
    
    def optimize_with_calculator(self, atoms: Atoms, calculator: Calculator, 
                               fmax: float, max_steps: int, 
                               label: str = "optimization") -> Tuple[bool, Atoms]:
        """
        Optimize structure with given calculator
        
        Args:
            atoms: Initial structure
            calculator: ASE calculator
            fmax: Force convergence criterion
            max_steps: Maximum optimization steps
            label: Label for logging
            
        Returns:
            (converged, optimized_atoms)
        """
        atoms = atoms.copy()
        atoms.calc = calculator
        
        # Choose optimizer
        if label == "DFT":
            opt = LBFGS(atoms, logfile=f"{self.molecule_name}_{label}.log")
        else:
            opt = BFGS(atoms, logfile=f"{self.molecule_name}_{label}.log")
        
        self.log(f"Starting {label} optimization (fmax={fmax} eV/Å)")
        
        # Custom convergence monitoring
        step = 0
        converged = False
        
        def check_convergence():
            nonlocal step, converged
            forces = atoms.get_forces()
            fmax_current = np.sqrt((forces**2).sum(axis=1).max())
            
            if step % 10 == 0:
                energy = atoms.get_potential_energy()
                self.log(f"  Step {step}: E = {energy:.4f} eV, Fmax = {fmax_current:.4f} eV/Å")
            
            step += 1
            
            if fmax_current < fmax:
                converged = True
                return True
            
            return step >= max_steps
        
        # Run optimization
        try:
            opt.run(fmax=fmax, steps=max_steps, callback=check_convergence)
            
            if converged:
                self.log(f"{label} optimization converged in {step} steps", "SUCCESS")
            else:
                self.log(f"{label} optimization stopped at {step} steps (not converged)", "WARNING")
                
        except Exception as e:
            self.log(f"{label} optimization failed: {e}", "ERROR")
            converged = False
        
        return converged, atoms
    
    def database_lookup(self) -> Optional[Atoms]:
        """
        Look up molecule structure in databases
        
        Returns:
            Atoms object or None
        """
        self.log(f"Searching databases for {self.molecule_name}")
        
        # Common molecules shortcuts
        common_molecules = {
            'water': 'H2O',
            'h2o': 'H2O',
            'tma': 'trimethylaluminum',
            'trimethylaluminum': 'Al(CH3)3',
            'ammonia': 'NH3',
            'methane': 'CH4'
        }
        
        search_name = common_molecules.get(self.molecule_name.lower(), self.molecule_name)
        
        # Try PubChem
        try:
            # Search by name
            search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{search_name}/cids/JSON"
            response = requests.get(search_url, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
                    cid = data['IdentifierList']['CID'][0]
                    self.log(f"Found PubChem CID: {cid}")
                    
                    # Get 3D structure
                    sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{cid}/record/SDF?record_type=3d"
                    sdf_response = requests.get(sdf_url, timeout=10)
                    
                    if sdf_response.status_code == 200:
                        # Save SDF temporarily
                        sdf_file = f"{self.molecule_name}_pubchem.sdf"
                        with open(sdf_file, 'w') as f:
                            f.write(sdf_response.text)
                        
                        # Read with ASE
                        atoms = read(sdf_file, format='sdf')
                        os.remove(sdf_file)  # Clean up
                        
                        self.log(f"Retrieved 3D structure from PubChem ({len(atoms)} atoms)", "SUCCESS")
                        return atoms
                        
        except Exception as e:
            self.log(f"Database lookup failed: {e}", "WARNING")
        
        return None
    
    def build_from_scratch(self) -> Optional[Atoms]:
        """
        Build simple molecules from scratch
        
        Returns:
            Atoms object or None
        """
        self.log("Building molecule from scratch")
        
        # Simple molecule library
        molecules = {
            'h2o': {
                'symbols': ['O', 'H', 'H'],
                'positions': [
                    [0.0, 0.0, 0.0],
                    [0.757, 0.586, 0.0],
                    [-0.757, 0.586, 0.0]
                ]
            },
            'nh3': {
                'symbols': ['N', 'H', 'H', 'H'],
                'positions': [
                    [0.0, 0.0, 0.0],
                    [0.942, 0.0, -0.331],
                    [-0.471, 0.816, -0.331],
                    [-0.471, -0.816, -0.331]
                ]
            },
            'ch4': {
                'symbols': ['C', 'H', 'H', 'H', 'H'],
                'positions': [
                    [0.0, 0.0, 0.0],
                    [0.629, 0.629, 0.629],
                    [-0.629, -0.629, 0.629],
                    [-0.629, 0.629, -0.629],
                    [0.629, -0.629, -0.629]
                ]
            }
        }
        
        mol_key = self.molecule_name.lower()
        if mol_key in molecules:
            mol_data = molecules[mol_key]
            atoms = Atoms(symbols=mol_data['symbols'], 
                         positions=mol_data['positions'])
            self.log(f"Built {mol_key} from template", "SUCCESS")
            return atoms
        
        self.log(f"No template available for {self.molecule_name}", "WARNING")
        return None
    
    def optimize(self, initial_structure: Optional[Atoms] = None) -> Tuple[Atoms, str]:
        """
        Main optimization workflow
        
        Args:
            initial_structure: Optional starting structure
            
        Returns:
            (optimized_atoms, method_used)
        """
        self.log(f"Starting optimization workflow for {self.molecule_name}")
        
        # 1. Get initial structure
        if initial_structure is not None:
            atoms = initial_structure.copy()
            self.log("Using provided initial structure")
        else:
            # Try database lookup
            atoms = self.database_lookup()
            if atoms is None:
                # Try building from scratch
                atoms = self.build_from_scratch()
                if atoms is None:
                    self.log("Failed to obtain initial structure", "ERROR")
                    self.save_log()
                    raise ValueError(f"Cannot find or build {self.molecule_name}")
        
        # Save initial structure
        write(f"{self.molecule_name}_initial.xyz", atoms)
        
        # 2. Pre-optimize with molecular mechanics (if not from database)
        mm_calc = self.get_molecular_ff_calculator()
        if mm_calc is not None and initial_structure is not None:
            self.log("Pre-optimizing with molecular mechanics")
            converged, atoms = self.optimize_with_calculator(
                atoms, mm_calc, self.mm_fmax, self.mm_max_steps, "MM"
            )
            write(f"{self.molecule_name}_mm_optimized.xyz", atoms)
        
        # 3. DFT optimization
        dft_calc = self.get_dft_calculator(atoms)
        if dft_calc is not None:
            self.log("Performing DFT optimization")
            converged, atoms = self.optimize_with_calculator(
                atoms, dft_calc, self.dft_fmax, self.dft_max_steps, "DFT"
            )
            
            if converged:
                write(f"{self.molecule_name}_dft_optimized.xyz", atoms)
                method = "DFT/PBE"
            else:
                self.log("DFT did not fully converge, results may be approximate", "WARNING")
                method = "DFT/PBE (partial)"
        else:
            self.log("DFT not available, using MM-optimized structure", "WARNING")
            method = "MM"
        
        # 4. Final analysis
        if hasattr(atoms, 'calc') and atoms.calc is not None:
            energy = atoms.get_potential_energy()
            forces = atoms.get_forces()
            fmax = np.sqrt((forces**2).sum(axis=1).max())
            
            self.log(f"\nFinal Results:", "SUCCESS")
            self.log(f"  Method: {method}")
            self.log(f"  Energy: {energy:.4f} eV")
            self.log(f"  Max force: {fmax:.4f} eV/Å")
            self.log(f"  Atoms: {len(atoms)}")
        
        # Save final structure and log
        write(f"{self.molecule_name}_optimized.xyz", atoms)
        self.save_log()
        
        return atoms, method

# Convenience function for backward compatibility
def SmartOptimizer(molecule_name: str, verbose: bool = True):
    """Create RealisticOptimizer instance (backward compatible)"""
    return RealisticOptimizer(molecule_name, verbose=verbose)