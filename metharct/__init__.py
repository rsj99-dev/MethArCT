#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MethArCT - Methanogenic Archaea Metabolic Pathway Analysis Tool

A comprehensive bioinformatics tool for analyzing metabolic pathways
in methanogenic archaea using Diamond BLAST and other analysis tools.

Author: MethArCT Development Team
Version: 0.6.3
License: GPL-3.0
"""

__version__ = "0.6.3"
__author__ = "MethArCT Development Team"
__email__ = "rsj1999@njtech.edu.cn"
__license__ = "GPL-3.0"
__description__ = "Methanogenic Archaea Metabolic Pathway Analysis Tool"

def __getattr__(name):
    """Lazy import of heavy modules to avoid import errors when dependencies are missing."""
    _lazy_imports = {
        'DiamondAnalyzer': '.core.diamond_analyzer',
        'TomeAnalyzer': '.core.tome_analyzer',
        'CheckM2Analyzer': '.core.checkm2_analyzer',
        'AntibioticAnalyzer': '.core.antibiotic_analyzer',
        'PathwayPredictor': '.core.pathway_predictor',
        'PHAnalyzer': '.core.ph_analyzer',
        'Config': '.utils.config',
        'setup_logger': '.utils.logger',
    }
    if name in _lazy_imports:
        import importlib
        module = importlib.import_module(_lazy_imports[name], __package__)
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

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

def get_version():
    """Get package version."""
    return __version__

def get_info():
    """Get package information."""
    return {
        "version": __version__,
        "author": __author__,
        "email": __email__,
        "license": __license__,
        "description": __description__,
    }

def check_installation():
    """Check if MethArCT is properly installed."""
    try:
        from .utils.config import Config
        from .core.diamond_analyzer import DiamondAnalyzer
        return True
    except ImportError:
        return False
