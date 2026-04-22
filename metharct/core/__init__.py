#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MethArCT Core Analysis Modules

This module contains the core analysis components for MethArCT:
- DiamondAnalyzer: Metabolic pathway prediction and salt tolerance analysis
- TomeAnalyzer: Optimal growth temperature prediction
- CheckM2Analyzer: Genome quality assessment and cultivability evaluation
- PathwayPredictor: Comprehensive analysis integrating all tools
"""

from .diamond_analyzer import DiamondAnalyzer
from .tome_analyzer import TomeAnalyzer
from .checkm2_analyzer import CheckM2Analyzer
from .pathway_predictor import PathwayPredictor

__all__ = [
    'DiamondAnalyzer',
    'TomeAnalyzer', 
    'CheckM2Analyzer',
    'PathwayPredictor'
]

# Version information
__version__ = '1.0.0'
__author__ = 'MethArCT Development Team'