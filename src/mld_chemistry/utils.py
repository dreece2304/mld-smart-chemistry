"""
Utility functions for MLD Chemistry Framework
"""

import os
import sys
import platform
import subprocess
from typing import Dict, Any, Optional
from pathlib import Path

def get_system_info() -> Dict[str, Any]:
    """Get comprehensive system information"""
    
    info = {
        'platform': platform.system(),
        'platform_version': platform.platform(),
        'python_version': platform.python_version(),
        'architecture': platform.architecture()[0],
        'processor': platform.processor(),
        'cpu_count': os.cpu_count(),
        'is_wsl': 'Microsoft' in platform.uname().release,
        'conda_env': os.environ.get('CONDA_DEFAULT_ENV', 'none'),
        'gpaw_setup_path': os.environ.get('GPAW_SETUP_PATH', ''),
        'mpl_backend': os.environ.get('MPLBACKEND', 'default')
    }
    
    # Memory information (if available)
    try:
        import psutil
        memory = psutil.virtual_memory()
        info['memory_total_gb'] = round(memory.total / (1024**3), 1)
        info['memory_available_gb'] = round(memory.available / (1024**3), 1)
    except ImportError:
        info['memory_total_gb'] = 'unknown'
        info['memory_available_gb'] = 'unknown'
    
    return info

def setup_environment(
    force_headless: bool = True,
    set_omp_threads: Optional[int] = None,
    verbose: bool = True
) -> bool:
    """Setup optimal environment for MLD calculations"""
    
    if verbose:
        print("🔧 Setting up MLD environment...")
    
    # Force headless backends for WSL compatibility
    if force_headless:
        os.environ['MPLBACKEND'] = 'Agg'
        os.environ['QT_QPA_PLATFORM'] = 'offscreen'
        if verbose:
            print("   ✅ Set headless visualization backends")
    
    # Set OpenMP threads for parallel efficiency
    if set_omp_threads is not None:
        os.environ['OMP_NUM_THREADS'] = str(set_omp_threads)
        if verbose:
            print(f"   ✅ Set OMP_NUM_THREADS={set_omp_threads}")
    
    # WSL-specific optimizations
    if 'Microsoft' in platform.uname().release:
        os.environ['LIBGL_ALWAYS_INDIRECT'] = '1'
        os.environ['MESA_GL_VERSION_OVERRIDE'] = '3.3'
        os.environ['OMPI_MCA_btl_vader_single_copy_mechanism'] = 'none'
        if verbose:
            print("   ✅ Applied WSL2 optimizations")
    
    # Ensure GPAW setup path is set
    if not os.environ.get('GPAW_SETUP_PATH'):
        # Try common locations
        possible_paths = [
            '/home/dreece23/miniconda3/envs/mld_modeling/share/gpaw',
            os.path.expanduser('~/.local/share/gpaw/gpaw-setups-24.11.0'),
            os.path.expanduser('~/gpaw-setups')
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                os.environ['GPAW_SETUP_PATH'] = path
                if verbose:
                    print(f"   ✅ Found GPAW datasets: {path}")
                break
        else:
            if verbose:
                print("   ⚠️  GPAW_SETUP_PATH not found - may need manual setup")
    
    return True

def validate_installation(verbose: bool = True) -> Dict[str, bool]:
    """Validate MLD installation components"""
    
    if verbose:
        print("🔍 Validating MLD installation...")
    
    results = {}
    
    # Test core imports
    core_packages = ['numpy', 'scipy', 'matplotlib', 'ase', 'gpaw', 'tqdm']
    
    for package in core_packages:
        try:
            __import__(package)
            results[package] = True
            if verbose:
                print(f"   ✅ {package}")
        except ImportError:
            results[package] = False
            if verbose:
                print(f"   ❌ {package}")
    
    # Test GPAW setup
    try:
        import gpaw
        gpaw_path = os.environ.get('GPAW_SETUP_PATH', '')
        if gpaw_path and os.path.exists(gpaw_path):
            setup_files = list(Path(gpaw_path).glob('*.gz'))
            results['gpaw_datasets'] = len(setup_files) > 0
            if verbose:
                print(f"   ✅ GPAW datasets: {len(setup_files)} files")
        else:
            results['gpaw_datasets'] = False
            if verbose:
                print(f"   ❌ GPAW datasets not found")
    except ImportError:
        results['gpaw_datasets'] = False
    
    # Test MPI
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        results['mpi'] = True
        if verbose:
            print(f"   ✅ MPI: {comm.Get_size()} process(es)")
    except ImportError:
        results['mpi'] = False
        if verbose:
            print(f"   ❌ MPI not available")
    
    # Overall success
    results['overall'] = all([
        results.get('numpy', False),
        results.get('scipy', False), 
        results.get('ase', False),
        results.get('gpaw', False),
        results.get('gpaw_datasets', False)
    ])
    
    if verbose:
        success_rate = sum(results.values()) / len(results) * 100
        print(f"   📊 Overall: {success_rate:.0f}% components working")
    
    return results

def create_results_directory(base_name: str = "mld_results") -> Path:
    """Create timestamped results directory"""
    
    from datetime import datetime
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_dir = Path(f"{base_name}_{timestamp}")
    results_dir.mkdir(exist_ok=True)
    
    return results_dir

def get_molecule_info(atoms) -> Dict[str, Any]:
    """Get information about molecular structure"""
    
    if atoms is None:
        return {}
    
    symbols = atoms.get_chemical_symbols()
    unique_elements = list(set(symbols))
    element_counts = {elem: symbols.count(elem) for elem in unique_elements}
    
    info = {
        'total_atoms': len(atoms),
        'elements': unique_elements,
        'element_counts': element_counts,
        'molecular_formula': ''.join([f"{elem}{count}" if count > 1 else elem 
                                    for elem, count in sorted(element_counts.items())])
    }
    
    # Try to get energy if available
    try:
        info['energy_ev'] = atoms.get_potential_energy()
    except:
        info['energy_ev'] = None
    
    return info

def check_disk_space(path: str = ".", min_gb: float = 1.0) -> bool:
    """Check if sufficient disk space is available"""
    
    try:
        import shutil
        total, used, free = shutil.disk_usage(path)
        free_gb = free / (1024**3)
        return free_gb >= min_gb
    except:
        return True  # Assume OK if can't check

def print_system_summary():
    """Print comprehensive system summary"""
    
    info = get_system_info()
    
    print("🖥️  System Information")
    print("=" * 40)
    print(f"Platform: {info['platform']} ({info['architecture']})")
    print(f"Python: {info['python_version']}")
    print(f"CPU cores: {info['cpu_count']}")
    print(f"Memory: {info['memory_total_gb']} GB total, {info['memory_available_gb']} GB available")
    print(f"Conda env: {info['conda_env']}")
    
    if info['is_wsl']:
        print("Environment: WSL2 (optimized for headless)")
    else:
        print("Environment: Native Linux")
    
    print("=" * 40)