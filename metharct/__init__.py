#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MethArCT - Methanogenic Archaea Metabolic Pathway Analysis Tool

A comprehensive bioinformatics tool for analyzing metabolic pathways
in methanogenic archaea using Diamond BLAST and other analysis tools.

Author: MethArCT Development Team
Version: 0.5.5
License: MIT
"""

__version__ = "0.5.5"
__author__ = "MethArCT Development Team"
__email__ = "rsj1999@njtech.edu.cn"
__license__ = "MIT"
__description__ = "Methanogenic Archaea Metabolic Pathway Analysis Tool"

# Import main components
from .core.diamond_analyzer import DiamondAnalyzer
from .core.tome_analyzer import TomeAnalyzer
from .core.checkm2_analyzer import CheckM2Analyzer
from .core.antibiotic_analyzer import AntibioticAnalyzer
from .core.pathway_predictor import PathwayPredictor
from .core.ph_analyzer import PHAnalyzer
from .utils.config import Config
from .utils.logger import setup_logger

# Define public API
__all__ = [
    "DiamondAnalyzer",
    "TomeAnalyzer", 
    "CheckM2Analyzer",
    "AntibioticAnalyzer",
    "PathwayPredictor",
    "PHAnalyzer",
    "Config",
    "setup_logger",
    "__version__",
    "__author__",
    "__email__",
    "__license__",
    "__description__",
]

# Package metadata
PACKAGE_NAME = "metharct"
PACKAGE_VERSION = __version__
PACKAGE_AUTHOR = __author__
PACKAGE_EMAIL = __email__
PACKAGE_LICENSE = __license__
PACKAGE_DESCRIPTION = __description__

# Version info
VERSION_INFO = {
    "version": __version__,
    "author": __author__,
    "email": __email__,
    "license": __license__,
    "description": __description__,
}

def get_version():
    """Get package version."""
    return __version__

def get_info():
    """Get package information."""
    return VERSION_INFO.copy()

def check_installation():
    """Check if MethArCT is properly installed."""
    try:
        from .utils.config import Config
        from .core.diamond_analyzer import DiamondAnalyzer
        return True
    except ImportError:
        return False

# Package initialization message
def _print_welcome():
    """Print welcome message when package is imported."""
    print(f"Welcome to MethArCT v0.5.5")
    print("Methanogenic Archaea Metabolic Pathway Analysis Tool")
    print("For help and documentation, visit: https://github.com/rsj99-dev/MethArCT")

# Optional: Print welcome message (can be disabled)
# _print_welcome()
