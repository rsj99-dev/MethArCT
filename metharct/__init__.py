#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MethArCT - Methanogenic Archaeal Culturomics Toolkit

A comprehensive toolkit for analyzing methanogenic archaea culturomics data.
Includes metabolic pathway prediction, salt tolerance assessment, 
optimal growth temperature prediction, and cultivability evaluation.

Author: Chen Jiuyu, Ren Shijie, Zhou Jun
Email: rsj1999@njtech.edu.cn
Version: 1.0.0
"""

__version__ = "1.0.0"
__author__ = "Ren Shijie"
__email__ = "rsj1999@njtech.edu.cn"
__description__ = "Methanogenic Archaeal Culturomics Toolkit"

# Import main classes
from .core.diamond_analyzer import DiamondAnalyzer
from .core.tome_analyzer import TomeAnalyzer
from .core.checkm2_analyzer import CheckM2Analyzer
from .core.pathway_predictor import PathwayPredictor
from .utils.config import Config
from .utils.logger import setup_logger

# Main analysis class
from .analyzer import MethArCTAnalyzer

__all__ = [
    'DiamondAnalyzer',
    'TomeAnalyzer', 
    'CheckM2Analyzer',
    'PathwayPredictor',
    'MethArCTAnalyzer',
    'Config',
    'setup_logger'
]