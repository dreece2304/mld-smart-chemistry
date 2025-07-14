#!/usr/bin/env python3
"""
Setup script for MLD Smart Chemistry Framework
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README for long description
readme_file = Path(__file__).parent / "README_NEW.md"
if readme_file.exists():
    with open(readme_file, "r", encoding="utf-8") as f:
        long_description = f.read()
else:
    long_description = "MLD Smart Chemistry Framework for computational chemistry simulations"

# Read requirements
requirements_file = Path(__file__).parent / "requirements.txt"
if requirements_file.exists():
    with open(requirements_file, "r") as f:
        requirements = [line.strip() for line in f if line.strip() and not line.startswith("#")]
else:
    requirements = [
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "matplotlib>=3.5.0",
        "ase>=3.22.0",
        "gpaw>=22.8.0",
        "mpi4py>=3.1.0",
        "tqdm>=4.0.0",
        "psutil>=5.8.0",
        "requests>=2.25.0"
    ]

setup(
    name="mld-smart-chemistry",
    version="1.0.0",
    author="MLD Chemistry Team",
    author_email="contact@mld-chemistry.org",
    description="Smart computational chemistry framework for MLD simulations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/mld-smart-chemistry",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/mld-smart-chemistry/issues",
        "Documentation": "https://github.com/yourusername/mld-smart-chemistry/docs",
        "Source": "https://github.com/yourusername/mld-smart-chemistry",
    },
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires=">=3.8,<3.12",
    install_requires=requirements,
    extras_require={
        "visualization": ["plotly>=5.0.0", "py3dmol>=2.0.0", "nglview>=3.0.0"],
        "dev": ["pytest>=6.0", "black", "flake8", "mypy"],
        "docs": ["sphinx", "sphinx-rtd-theme", "myst-parser"],
    },
    entry_points={
        "console_scripts": [
            "mld-optimize=mld_chemistry.smart_optimizer:main",
            "mld-visualize=mld_chemistry.visualization:main",
            "mld-validate=install.validate_installation:main",
        ],
    },
    include_package_data=True,
    package_data={
        "mld_chemistry": ["data/*.xyz", "templates/*.html"],
    },
    zip_safe=False,
    keywords=[
        "computational chemistry",
        "density functional theory", 
        "molecular layer deposition",
        "GPAW",
        "ASE",
        "DFT",
        "materials science"
    ],
)