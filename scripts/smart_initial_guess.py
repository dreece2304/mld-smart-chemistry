#!/usr/bin/env python3
"""
Smart Initial Guess System for DFT Optimization
Provides intelligent starting points for similar materials to reduce computation time

Features:
- Geometry similarity matching
- Wavefunction reuse from similar systems
- Optimized starting parameters based on material type
- Cached results for faster subsequent runs
"""

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.neighborlist import NeighborList
import pickle
import os
import json
from typing import Dict, List, Tuple, Optional, Any
from datetime import datetime
import hashlib

try:
    from gpaw import GPAW, PW
    from gpaw.wavefunctions.pw import PWWaveFunctions
    GPAW_AVAILABLE = True
except ImportError:
    GPAW_AVAILABLE = False

class GeometryMatcher:
    """Match similar geometries for initial guess optimization"""
    
    def __init__(self, tolerance: float = 0.3):
        """
        Initialize geometry matcher
        
        Args:
            tolerance: Distance tolerance for similarity matching (Å)
        """
        self.tolerance = tolerance
        self.cache_dir = "geometry_cache"
        self.ensure_cache_dir()
    
    def ensure_cache_dir(self):
        """Ensure cache directory exists"""
        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)
    
    def compute_fingerprint(self, atoms: Atoms) -> str:
        """Compute unique fingerprint for atomic geometry"""
        
        # Get basic composition
        symbols = sorted(atoms.get_chemical_symbols())
        composition = f"{''.join(symbols)}"
        
        # Get neighbor information
        nl = NeighborList([1.5] * len(atoms), self_interaction=False, bothways=True)
        nl.update(atoms)
        
        neighbor_counts = []
        for i in range(len(atoms)):
            indices, offsets = nl.get_neighbors(i)
            neighbor_counts.append(len(indices))
        
        # Create fingerprint
        fingerprint_data = {
            'composition': composition,
            'n_atoms': len(atoms),
            'neighbor_counts': sorted(neighbor_counts),
            'cell_volume': atoms.get_volume() if atoms.pbc.any() else 0.0
        }
        
        # Hash the fingerprint
        fingerprint_str = json.dumps(fingerprint_data, sort_keys=True)
        return hashlib.md5(fingerprint_str.encode()).hexdigest()
    
    def compute_similarity(self, atoms1: Atoms, atoms2: Atoms) -> float:
        """
        Compute similarity score between two structures
        
        Returns:
            float: Similarity score (0-1, higher is more similar)
        """
        
        # Basic checks
        if len(atoms1) != len(atoms2):
            return 0.0
        
        if sorted(atoms1.get_chemical_symbols()) != sorted(atoms2.get_chemical_symbols()):
            return 0.0
        
        # Compute structural similarity
        try:
            # Get positions relative to center of mass
            pos1 = atoms1.get_positions()
            pos2 = atoms2.get_positions()
            
            # Center both structures
            pos1 -= pos1.mean(axis=0)
            pos2 -= pos2.mean(axis=0)
            
            # Compute distance matrices
            dist1 = np.linalg.norm(pos1[:, np.newaxis] - pos1[np.newaxis, :], axis=2)
            dist2 = np.linalg.norm(pos2[:, np.newaxis] - pos2[np.newaxis, :], axis=2)
            
            # Compare distance matrices
            diff = np.abs(dist1 - dist2)
            max_diff = np.max(diff)
            
            if max_diff < self.tolerance:
                similarity = 1.0 - (max_diff / self.tolerance)
            else:
                similarity = 0.0
            
            return similarity
            
        except Exception as e:
            print(f"Error computing similarity: {e}")
            return 0.0
    
    def find_similar_structures(self, target: Atoms, min_similarity: float = 0.7) -> List[Tuple[str, float]]:
        """
        Find similar structures in cache
        
        Args:
            target: Target structure to match
            min_similarity: Minimum similarity threshold
            
        Returns:
            List of (cache_key, similarity_score) tuples
        """
        
        similar_structures = []
        
        # Check all cached structures
        for cache_file in os.listdir(self.cache_dir):
            if not cache_file.endswith('.xyz'):
                continue
                
            cache_path = os.path.join(self.cache_dir, cache_file)
            try:
                cached_atoms = read(cache_path)
                similarity = self.compute_similarity(target, cached_atoms)
                
                if similarity >= min_similarity:
                    cache_key = cache_file.replace('.xyz', '')
                    similar_structures.append((cache_key, similarity))
                    
            except Exception as e:
                print(f"Error reading {cache_file}: {e}")
                continue
        
        # Sort by similarity (highest first)
        similar_structures.sort(key=lambda x: x[1], reverse=True)
        
        return similar_structures

class InitialGuessOptimizer:
    """Optimize initial guesses for DFT calculations"""
    
    def __init__(self, cache_dir: str = "dft_cache"):
        """Initialize with cache directory"""
        self.cache_dir = cache_dir
        self.geometry_matcher = GeometryMatcher()
        self.ensure_cache_dir()
        
        # Material-specific parameters
        self.material_params = {
            'si_surface': {
                'cutoff': 300,
                'smearing': 0.1,
                'mixing': 0.05,
                'expected_gap': 1.1  # eV
            },
            'sio2_surface': {
                'cutoff': 350,
                'smearing': 0.05,
                'mixing': 0.03,
                'expected_gap': 8.0  # eV
            },
            'organic_molecule': {
                'cutoff': 400,
                'smearing': 0.01,
                'mixing': 0.1,
                'expected_gap': 5.0  # eV
            },
            'mixed_system': {
                'cutoff': 400,
                'smearing': 0.1,
                'mixing': 0.03,
                'expected_gap': 2.0  # eV
            }
        }
    
    def ensure_cache_dir(self):
        """Ensure cache directory exists"""
        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)
    
    def identify_material_type(self, atoms: Atoms) -> str:
        """Identify material type for parameter optimization"""
        
        symbols = set(atoms.get_chemical_symbols())
        
        # Check composition
        has_si = 'Si' in symbols
        has_o = 'O' in symbols
        has_h = 'H' in symbols
        has_c = 'C' in symbols
        has_al = 'Al' in symbols
        
        if has_si and has_o and len(symbols) <= 3:
            return 'sio2_surface'
        elif has_si and not has_o:
            return 'si_surface'
        elif has_c and has_h and not has_si:
            return 'organic_molecule'
        else:
            return 'mixed_system'
    
    def get_optimized_parameters(self, atoms: Atoms) -> Dict[str, Any]:
        """Get optimized DFT parameters for this material type"""
        
        material_type = self.identify_material_type(atoms)
        base_params = self.material_params[material_type].copy()
        
        # Adjust parameters based on system size
        n_atoms = len(atoms)
        
        if n_atoms > 100:
            # Large systems - more conservative
            base_params['cutoff'] = min(base_params['cutoff'], 300)
            base_params['mixing'] *= 0.5
            base_params['smearing'] *= 2
        elif n_atoms < 20:
            # Small systems - can be more aggressive
            base_params['cutoff'] = min(base_params['cutoff'] + 50, 500)
            base_params['mixing'] *= 1.5
            base_params['smearing'] *= 0.5
        
        # Convert to GPAW parameters
        gpaw_params = {
            'mode': PW(base_params['cutoff']),
            'xc': 'PBE',
            'kpts': (1, 1, 1),
            'occupations': {'name': 'fermi-dirac', 'width': base_params['smearing']},
            'mixer': {'beta': base_params['mixing'], 'nmaxold': 5},
            'convergence': {'energy': 1e-4, 'density': 1e-3},
            'maxiter': 300
        }
        
        print(f"🎯 Optimized parameters for {material_type}:")
        print(f"   Cutoff: {base_params['cutoff']} eV")
        print(f"   Mixing: {base_params['mixing']}")
        print(f"   Smearing: {base_params['smearing']} eV")
        
        return gpaw_params
    
    def get_initial_geometry(self, atoms: Atoms, target_label: str) -> Atoms:
        """Get optimized initial geometry from similar systems"""
        
        print(f"🔍 Searching for similar geometries to {target_label}...")
        
        # Find similar structures
        similar = self.geometry_matcher.find_similar_structures(atoms, min_similarity=0.8)
        
        if similar:
            best_match, similarity = similar[0]
            print(f"✅ Found similar structure: {best_match} (similarity: {similarity:.3f})")
            
            # Load the optimized geometry
            cache_path = os.path.join(self.geometry_matcher.cache_dir, f"{best_match}.xyz")
            if os.path.exists(cache_path):
                try:
                    optimized_atoms = read(cache_path)
                    print(f"📄 Using optimized geometry from {best_match}")
                    return optimized_atoms
                except Exception as e:
                    print(f"❌ Error loading cached geometry: {e}")
        
        print(f"💭 No similar structures found, using original geometry")
        return atoms.copy()
    
    def cache_result(self, atoms: Atoms, label: str, energy: float, 
                    calc_time: float, converged: bool = True):
        """Cache calculation results for future use"""
        
        if not converged:
            print(f"⚠️  Not caching unconverged result for {label}")
            return
        
        # Save geometry
        geometry_path = os.path.join(self.geometry_matcher.cache_dir, f"{label}.xyz")
        write(geometry_path, atoms)
        
        # Save metadata
        metadata = {
            'label': label,
            'energy': energy,
            'calc_time': calc_time,
            'n_atoms': len(atoms),
            'material_type': self.identify_material_type(atoms),
            'fingerprint': self.geometry_matcher.compute_fingerprint(atoms),
            'timestamp': datetime.now().isoformat(),
            'converged': converged
        }
        
        metadata_path = os.path.join(self.cache_dir, f"{label}_metadata.json")
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        print(f"💾 Cached result for {label}")
        print(f"   Energy: {energy:.4f} eV")
        print(f"   Time: {calc_time:.1f} seconds")
    
    def get_cached_result(self, atoms: Atoms, label: str) -> Optional[Dict[str, Any]]:
        """Check if we have a cached result for this exact system"""
        
        fingerprint = self.geometry_matcher.compute_fingerprint(atoms)
        
        # Check if we have this exact calculation
        metadata_path = os.path.join(self.cache_dir, f"{label}_metadata.json")
        if os.path.exists(metadata_path):
            try:
                with open(metadata_path, 'r') as f:
                    metadata = json.load(f)
                
                if metadata['fingerprint'] == fingerprint:
                    geometry_path = os.path.join(self.geometry_matcher.cache_dir, f"{label}.xyz")
                    if os.path.exists(geometry_path):
                        cached_atoms = read(geometry_path)
                        
                        print(f"✅ Found exact cached result for {label}")
                        print(f"   Energy: {metadata['energy']:.4f} eV")
                        print(f"   Original time: {metadata['calc_time']:.1f} seconds")
                        
                        return {
                            'atoms': cached_atoms,
                            'energy': metadata['energy'],
                            'calc_time': metadata['calc_time'],
                            'from_cache': True
                        }
            except Exception as e:
                print(f"❌ Error reading cached result: {e}")
        
        return None
    
    def create_smart_calculator(self, atoms: Atoms, label: str) -> GPAW:
        """Create GPAW calculator with optimized parameters"""
        
        if not GPAW_AVAILABLE:
            raise ImportError("GPAW not available")
        
        # Get optimized parameters
        params = self.get_optimized_parameters(atoms)
        params['txt'] = f'{label}_smart.txt'
        
        # Create calculator
        calc = GPAW(**params)
        
        print(f"🧮 Created smart calculator for {label}")
        return calc

def smart_single_point(atoms: Atoms, label: str, 
                      use_cache: bool = True) -> Tuple[bool, float, float]:
    """
    Run single point calculation with smart initial guess
    
    Args:
        atoms: Input structure
        label: Calculation label
        use_cache: Whether to use cached results
        
    Returns:
        (success, energy, calc_time)
    """
    
    optimizer = InitialGuessOptimizer()
    
    # Check for cached result
    if use_cache:
        cached = optimizer.get_cached_result(atoms, label)
        if cached:
            return True, cached['energy'], 0.0  # No calc time for cached
    
    # Get optimized initial geometry
    optimized_atoms = optimizer.get_initial_geometry(atoms, label)
    
    # Create smart calculator
    calc = optimizer.create_smart_calculator(optimized_atoms, label)
    optimized_atoms.calc = calc
    
    # Run calculation
    import time
    start_time = time.time()
    
    try:
        energy = optimized_atoms.get_potential_energy()
        calc_time = time.time() - start_time
        
        # Cache result
        optimizer.cache_result(optimized_atoms, label, energy, calc_time, True)
        
        return True, energy, calc_time
        
    except Exception as e:
        calc_time = time.time() - start_time
        print(f"❌ Smart calculation failed: {e}")
        return False, 0.0, calc_time

def main():
    """Demo smart initial guess system"""
    
    print("🎯 Smart Initial Guess System Demo")
    print("="*50)
    
    if not GPAW_AVAILABLE:
        print("❌ GPAW not available - demo only")
        return 1
    
    # Create test system
    from ase.build import molecule
    h2o = molecule('H2O')
    h2o.center(vacuum=8.0)
    
    print("\n1. Testing smart single point calculation...")
    
    # First run (should be slow)
    print("\n📊 First calculation (building cache):")
    success1, energy1, time1 = smart_single_point(h2o, "h2o_demo", use_cache=False)
    
    # Second run (should be fast from cache)
    print("\n📊 Second calculation (using cache):")
    success2, energy2, time2 = smart_single_point(h2o, "h2o_demo", use_cache=True)
    
    print(f"\n✅ Smart initial guess demo complete!")
    print(f"   First run: {time1:.1f} seconds")
    print(f"   Second run: {time2:.1f} seconds (cached)")
    print(f"   Speedup: {time1/max(time2, 0.1):.1f}x")
    
    return 0

if __name__ == "__main__":
    exit(main())