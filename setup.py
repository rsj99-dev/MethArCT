#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MethArCT - Methanogenic Archaea Metabolic Pathway Analysis Tool
Setup script for Python package distribution
"""

from setuptools import setup, find_packages
import os
import sys

# Read README file
def read_readme():
    try:
        with open("README.md", "r", encoding="utf-8") as fh:
            return fh.read()
    except FileNotFoundError:
        return "MethArCT - Methanogenic Archaea Metabolic Pathway Analysis Tool"

# Read requirements file
def read_requirements():
    try:
        with open("requirements.txt", "r", encoding="utf-8") as fh:
            return [line.strip() for line in fh if line.strip() and not line.startswith("#")]
    except FileNotFoundError:
        return [
            "biopython>=1.79",
            "pandas>=1.3.0",
            "numpy>=1.21.0",
            "scikit-learn>=1.0.0",
            "pyyaml>=5.4",
            "tqdm>=4.60.0",
            "requests>=2.25.0"
        ]

# Get version information
def get_version():
    """Read version from metharct/__init__.py dynamically."""
    version_file = os.path.join("metharct", "__init__.py")
    if os.path.exists(version_file):
        with open(version_file, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("__version__"):
                    return line.split("=")[1].strip().strip('"').strip("'")
    return "0.7.3"


def collect_database_files():
    """Collect reference databases for installation under share/metharct.

    The databases live in the repository root ``data/databases`` directory.
    Installing them with ``data_files`` makes ``pip install`` usable outside a
    source checkout; ``Config.resolve_database_base_dir`` looks them up in
    ``<prefix>/share/metharct/databases``.
    """
    data_files = []
    source_root = os.path.join("data", "databases")
    if not os.path.isdir(source_root):
        return data_files

    for root, _dirs, filenames in os.walk(source_root):
        if not filenames:
            continue
        rel_dir = os.path.relpath(root, source_root)
        if rel_dir == ".":
            target = os.path.join("share", "metharct", "databases")
        else:
            target = os.path.join("share", "metharct", "databases", rel_dir)
        data_files.append(
            (target, [os.path.join(root, name) for name in filenames])
        )

    return data_files


setup(
    name="metharct",
    version=get_version(),
    author="rsj99",
    author_email="rsj1999@njtech.edu.cn",
    description="Methanogenic Archaea Metabolic Pathway and Environmental Preference Analysis Tool",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/MethArCT/metharct",
    project_urls={
        "Bug Reports": "https://github.com/rsj99-dev/MethArCT/issues",
        "Source": "https://github.com/rsj99-dev/MethArCT",
        "Documentation": "https://github.com/rsj99-dev/MethArCT",
    },
    packages=find_packages(include=['metharct', 'metharct.*']),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=read_requirements(),
    extras_require={
        "dev": [
            "pytest>=7.0",
            "pytest-cov>=4.0",
            "black>=23.0",
            "flake8>=6.0",
            "mypy>=1.0",
        ],
        "docs": [
            "sphinx>=7.0",
            "sphinx-rtd-theme>=1.0",
            "myst-parser>=2.0",
        ],
        "full": [
            "diamond-annotator>=1.0",
            "checkm2>=1.0",
        ]
    },
    entry_points={
        "console_scripts": [
            "metharct=metharct.cli.main:main",
        ],
    },
    include_package_data=True,
    package_data={
        "metharct": [
            "data/*.yaml",
            "data/*.yml",
            # Supported if the databases are later moved inside the package
            "data/databases/*",
            "data/databases/**/*",
            "core/susha/models/*.pkl",
            "core/mlhuainanzi/models/*.pkl",
            "core/ph_predictor/models/*.joblib",
            "core/ph_predictor/hmm/*.joblib",
        ],
    },
    data_files=[
        ("share/metharct", ["metharct_config.yaml"]),
    ] + collect_database_files(),
    keywords=[
        "bioinformatics",
        "metagenomics",
        "methanogenic archaea",
        "metabolic pathways",
        "diamond",
        "genome analysis",
        "salinity prediction",
        "temperature prediction",
        "huainanzi",
        "pH prediction",
        "checkm2",
        "susha",
    ],
    zip_safe=False,
)
