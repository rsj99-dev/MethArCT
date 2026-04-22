#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MethArCT Command Line Interface

Provides command-line tools for running MethArCT analyses.
"""

from .main import main, create_parser
from .commands import (
    diamond_command,
    tome_command,
    checkm2_command,
    comprehensive_command
)

__all__ = [
    'main',
    'create_parser',
    'diamond_command',
    'tome_command', 
    'checkm2_command',
    'comprehensive_command'
]