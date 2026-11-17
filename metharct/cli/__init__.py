#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MethArCT Command Line Interface

Provides command-line tools for running MethArCT analyses.

Public names are resolved lazily so that ``python -m metharct.cli.main`` does
not import the module twice (which raised a ``runpy`` warning) and so that
``import metharct.cli`` stays cheap.
"""

from typing import Any

_LAZY_IMPORTS = {
    'main': '.main',
    'create_parser': '.main',
    'diamond_command': '.commands',
    'checkm2_command': '.commands',
    'comprehensive_command': '.commands',
    'susha_command': '.commands',
    'mlhuainanzi_command': '.commands',
    'ph_command': '.commands',
    'antibiotic_command': '.commands',
    'cultivation_command': '.commands',
    'batch_command': '.commands',
    'AI_COMMANDS': '.ai_commands',
    'add_ai_subparsers': '.ai_commands',
    'run_ai_command': '.ai_commands',
}

__all__ = list(_LAZY_IMPORTS)


def __getattr__(name: str) -> Any:
    """Resolve CLI entry points on first access (PEP 562)."""
    module_path = _LAZY_IMPORTS.get(name)
    if module_path is None:
        raise AttributeError("module {0!r} has no attribute {1!r}".format(__name__, name))
    import importlib

    module = importlib.import_module(module_path, __package__)
    return getattr(module, name)
