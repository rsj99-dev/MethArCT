#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MethArCT Core Analysis Modules

This module contains the core analysis components for MethArCT:
- DiamondAnalyzer: Metabolic pathway prediction and cultivability analysis
- CultivationAnalyzer: Metabolic pathway analysis for cultivability assessment
- SuShaAnalyzer: Salinity adaptation prediction
- HuainanziAnalyzer: Growth temperature range prediction
- PHAnalyzer: Growth pH preference prediction
- CheckM2Analyzer: Genome quality assessment and cultivability evaluation
- PathwayPredictor: Comprehensive analysis integrating all tools
"""

from .diamond_analyzer import DiamondAnalyzer
from .cultivation_analyzer import CultivationAnalyzer
from .susha_analyzer import SuShaAnalyzer
from .huainanzi_analyzer import HuainanziAnalyzer
from .ph_analyzer import PHAnalyzer
from .checkm2_analyzer import CheckM2Analyzer
from .antibiotic_analyzer import AntibioticAnalyzer
from .pathway_predictor import PathwayPredictor

__all__ = [
    'DiamondAnalyzer',
    'CultivationAnalyzer',
    'SuShaAnalyzer',
    'HuainanziAnalyzer',
    'PHAnalyzer',
    'CheckM2Analyzer',
    'AntibioticAnalyzer',
    'PathwayPredictor'
]

# Version information - use metharct.__version__ for the canonical version
__author__ = 'rsj1999'