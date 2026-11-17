# -*- coding: utf-8 -*-
"""Regression tests: the legacy CLI and Python API must stay intact.

These tests do not exercise the AI layer; they freeze the pre-existing
contract so the additive upgrade cannot silently change it.
"""

import argparse
import os

import pytest


def test_package_version_is_exposed():
    import metharct

    assert isinstance(metharct.__version__, str)
    assert metharct.__version__.count(".") >= 1


def test_legacy_public_names_resolve():
    import metharct

    for name in (
        "MethArCTAnalyzer",
        "DiamondAnalyzer",
        "CultivationAnalyzer",
        "HuainanziAnalyzer",
        "CheckM2Analyzer",
        "AntibioticAnalyzer",
        "PathwayPredictor",
        "PHAnalyzer",
        "Config",
        "setup_logger",
    ):
        assert getattr(metharct, name) is not None, name


def test_core_package_exports_analysis_classes():
    from metharct import core

    for name in (
        "DiamondAnalyzer",
        "CultivationAnalyzer",
        "SuShaAnalyzer",
        "HuainanziAnalyzer",
        "PHAnalyzer",
        "CheckM2Analyzer",
        "AntibioticAnalyzer",
        "PathwayPredictor",
    ):
        assert hasattr(core, name), name


def test_cli_parser_exposes_legacy_subcommands():
    from metharct.cli import create_parser

    parser = create_parser()
    assert isinstance(parser, argparse.ArgumentParser)

    actions = [action for action in parser._actions if isinstance(action, argparse._SubParsersAction)]
    assert actions, "the CLI must define subcommands"
    choices = actions[0].choices
    for command in (
        "comprehensive",
        "diamond",
        "checkm2",
        "cultivation",
        "susha",
        "mlhuainanzi",
        "ph",
        "antibiotic",
        "batch",
    ):
        assert command in choices, command


def test_legacy_comprehensive_arguments_are_unchanged():
    from metharct.cli import create_parser

    parser = create_parser()
    args = parser.parse_args(
        [
            "comprehensive",
            "input.faa",
            "-o",
            "results",
            "--skip-checkm2",
            "--skip-susha",
            "--threads",
            "8",
            "--evalue",
            "1e-10",
        ]
    )
    assert args.command == "comprehensive"
    assert args.input == "input.faa"
    assert args.output == "results"
    assert args.skip_checkm2 is True
    assert args.skip_susha is True
    assert args.threads == 8
    assert args.evalue == 1e-10


def test_legacy_diamond_defaults_are_unchanged():
    from metharct.cli import create_parser

    parser = create_parser()
    args = parser.parse_args(["diamond", "input.faa", "-o", "out"])
    assert args.threads == 4
    assert args.evalue == 1e-5
    assert args.max_target_seqs == 10


def test_default_config_still_loads_reference_paths():
    from metharct.utils.config import Config

    config = Config()
    base_dir = config.get("databases.base_dir")
    assert os.path.isdir(base_dir), base_dir
    assert config.get("tools.diamond.evalue") == 1e-5
    assert config.get("tools.diamond.threads") == 4


def test_analyzer_signatures_are_unchanged():
    import inspect

    from metharct.core.antibiotic_analyzer import AntibioticAnalyzer
    from metharct.core.checkm2_analyzer import CheckM2Analyzer
    from metharct.core.cultivation_analyzer import CultivationAnalyzer
    from metharct.core.diamond_analyzer import DiamondAnalyzer
    from metharct.core.huainanzi_analyzer import HuainanziAnalyzer
    from metharct.core.ph_analyzer import PHAnalyzer
    from metharct.core.susha_analyzer import SuShaAnalyzer

    expected = {
        DiamondAnalyzer.analyze_sequence: ["self", "input_file", "output_prefix", "databases"],
        CultivationAnalyzer.analyze_genome: ["self", "genome_file"],
        SuShaAnalyzer.predict_salinity: ["self", "input_file", "output_prefix"],
        HuainanziAnalyzer.predict_temperature: ["self", "input_file", "output_prefix"],
        PHAnalyzer.predict_ph: ["self", "input_file", "output_prefix"],
        AntibioticAnalyzer.predict_antibiotics: ["self", "input_file", "output_prefix"],
        CheckM2Analyzer.analyze_genome_quality: ["self", "input_path", "output_prefix", "input_type"],
    }
    for method, parameters in expected.items():
        assert list(inspect.signature(method).parameters) == parameters, method.__qualname__


def test_ai_layer_is_optional_for_legacy_use():
    """Importing the legacy API must not import the AI layer."""
    import subprocess
    import sys

    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    code = "import sys; import metharct; print('metharct.ai' in sys.modules)"
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=repo_root,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=120,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.decode().strip() == "False"


@pytest.mark.parametrize("module_name", ["metharct.ai", "metharct.ai.constants", "metharct.ai.errors"])
def test_ai_modules_import_without_optional_dependencies(module_name):
    __import__(module_name)
