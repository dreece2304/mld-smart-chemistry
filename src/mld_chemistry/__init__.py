"""
MLD Smart Chemistry Framework

A comprehensive computational chemistry framework for Molecular Layer Deposition (MLD) 
simulations using intelligent optimization strategies.
"""

__version__ = "1.0.0"
__author__ = "MLD Chemistry Team"
__email__ = "contact@mld-chemistry.org"

# Import main functions for easy access
from .smart_optimizer import SmartOptimizer, optimize_molecule
from .visualization import (
    visualize_molecule, 
    compare_molecules, 
    visualize_optimization_trajectory,
    create_html_report
)
from .utils import (
    setup_environment,
    validate_installation,
    get_system_info
)

# Define what gets imported with "from mld_chemistry import *"
__all__ = [
    'SmartOptimizer',
    'optimize_molecule',
    'visualize_molecule',
    'compare_molecules', 
    'visualize_optimization_trajectory',
    'create_html_report',
    'setup_environment',
    'validate_installation',
    'get_system_info'
]

# Package metadata
PACKAGE_INFO = {
    'name': 'mld_chemistry',
    'version': __version__,
    'description': 'Smart computational chemistry framework for MLD simulations',
    'supported_methods': ['DFT', 'Classical Force Fields', 'Database Lookup'],
    'supported_codes': ['GPAW', 'ASE', 'EMT'],
    'supported_databases': ['PubChem', 'PubChem 3D', 'NIST', 'ChEBI']
}

def get_package_info():
    """Get package information"""
    return PACKAGE_INFO.copy()

def print_welcome():
    """Print welcome message with package information"""
    print(f"🧪 MLD Smart Chemistry Framework v{__version__}")
    print("=" * 50)
    print("Intelligent computational chemistry for MLD simulations")
    print(f"Supported methods: {', '.join(PACKAGE_INFO['supported_methods'])}")
    print(f"Supported databases: {', '.join(PACKAGE_INFO['supported_databases'])}")
    print("=" * 50)