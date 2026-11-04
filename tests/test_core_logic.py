# -*- coding: utf-8 -*-
"""Regression tests for issues fixed in MethArCT v0.7.3.

These tests are intentionally independent of the DIAMOND executable so they
can run on any development machine.
"""

import os
import types

import pandas as pd
import pytest

from metharct.core.cultivation_analyzer import (
    CultivationAnalyzer,
    decode_encoded_filename,
)
from metharct.core.diamond_analyzer import DiamondAnalyzer
from metharct.utils.config import Config


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_CONFIG = os.path.join(REPO_ROOT, 'metharct_config.yaml')


def _hits(subject_ids, bitscore=200.0, evalue=1e-60, pident=90.0):
    """Build a DataFrame shaped like the DIAMOND outfmt 6 output."""
    rows = []
    for index, subject_id in enumerate(subject_ids):
        rows.append({
            'qseqid': f'q{index}',
            'sseqid': subject_id,
            'pident': pident,
            'length': 100,
            'mismatch': 0,
            'gapopen': 0,
            'qstart': 1,
            'qend': 100,
            'sstart': 1,
            'send': 100,
            'evalue': evalue,
            'bitscore': bitscore,
        })
    return pd.DataFrame(rows, columns=[
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
    ])


# ---------------------------------------------------------------- config

def test_repo_config_uses_unified_schema():
    """metharct_config.yaml must not turn cultivation into a dict for paths."""
    config = Config(REPO_CONFIG)

    paths = config.get_all_database_paths()
    assert 'CULTIVATION' in paths
    assert isinstance(paths['CULTIVATION'], str)

    for pathway_type in ('amino_acid', 'vitamin', 'nucleotide', 'atp'):
        directory = config.get_cultivation_pathway_dir(pathway_type)
        assert os.path.isdir(directory), f"{pathway_type}: {directory}"


def test_legacy_cultivation_dict_is_supported(tmp_path):
    legacy_config = tmp_path / 'legacy.yaml'
    legacy_config.write_text(
        "databases:\n"
        "  base_dir: '%s'\n"
        "  cultivation:\n"
        "    base_dir: cultivation\n"
        "    amino_acids: path.aa\n"
        "    vitamins: path.vc\n"
        "    nucleotides: path_hesuan\n"
        "    atp: path_atp\n"
        % os.path.join(REPO_ROOT, 'data', 'databases'),
        encoding='utf-8',
    )

    config = Config(str(legacy_config))
    paths = config.get_all_database_paths()
    assert paths['CULTIVATION'].endswith('cultivation')
    assert config.get_cultivation_pathway_dir('amino_acid').endswith('path.aa')


# ------------------------------------------------------------ cultivation

def test_decode_encoded_filename():
    encoded = '#U5f02#U4eae#U6c28#U9178#U751f#U7269#U5408#U6210-#U82cf#U6c28#U9178'
    assert decode_encoded_filename(encoded) == '异亮氨酸生物合成-苏氨酸'
    assert decode_encoded_filename('丙氨酸生物合成') == '丙氨酸生物合成'


def test_discover_pathways_decodes_names(tmp_path):
    encoded_file = tmp_path / ('#U5f02#U4eae#U6c28#U9178#U751f#U7269#U5408#U6210-'
                               '#U82cf#U6c28#U9178.fasta')
    encoded_file.write_text('>gene1\nMKT\n', encoding='utf-8')

    stub = types.SimpleNamespace(
        PATHWAY_NAME_MAP=CultivationAnalyzer.PATHWAY_NAME_MAP
    )
    pathways = CultivationAnalyzer._discover_pathways(
        stub, str(tmp_path), 'Amino acid'
    )
    assert pathways == {
        'Isoleucine biosynthesis-threonine': str(encoded_file)
    }


def test_extract_missing_amino_acids_uses_analyzer():
    from metharct.cli.commands import _extract_missing_amino_acids

    class StubAnalyzer:
        def analyze_amino_acid_biosynthesis(self, results):
            return {
                'summary': {
                    'missing_list': ['Lysine (Lys)', 'Arginine (Arg)']
                }
            }

    assert _extract_missing_amino_acids({}, StubAnalyzer()) == 'Lysine, Arginine'


# ------------------------------------------------------- diamond summary

def _make_diamond_analyzer(reference_counts, subject_counts):
    analyzer = DiamondAnalyzer.__new__(DiamondAnalyzer)
    analyzer.logger = types.SimpleNamespace(
        warning=lambda *a, **k: None,
        debug=lambda *a, **k: None,
    )
    analyzer.pathway_names = {}
    analyzer.reference_counts = reference_counts
    analyzer._reference_subject_counts = subject_counts
    return analyzer


def test_sr_denominator_uses_unique_subjects():
    analyzer = _make_diamond_analyzer({'SR': 3}, {'SR': 4})
    assert analyzer._get_reference_count('SR') == 4


def test_zero_hit_pathways_are_kept_and_sr_cannot_reach_100(tmp_path):
    input_faa = tmp_path / 'genome.faa'
    input_faa.write_text('>p1\nMKT\n', encoding='utf-8')

    analyzer = _make_diamond_analyzer(
        {'SR': 3, 'DSR': 8}, {'SR': 4, 'DSR': 11}
    )

    search_results = {
        'SR': {
            'database': 'SR',
            'total_hits': 3,
            'high_threshold_hits': _hits(['K18367', 'K17219', 'K17220']),
            'low_threshold_hits': _hits(['K18367', 'K17219', 'K17220']),
            'output_file': '',
        },
        'DSR': {
            'database': 'DSR',
            'total_hits': 0,
            'high_threshold_hits': pd.DataFrame(),
            'low_threshold_hits': pd.DataFrame(),
            'output_file': '',
        },
    }

    results = analyzer._process_results(
        search_results, input_faa, databases=['SR', 'DSR', 'ASR']
    )
    pathways = results['pathway_results']

    # 3 of 4 unique subjects -> 75 %, never 100 %
    assert pathways['SR']['reference_sequences'] == 4
    assert pathways['SR']['low_completeness_percentage'] == 75.0

    # Zero-hit and not-searched databases still produce rows
    assert pathways['DSR']['low_completeness_percentage'] == 0.0
    assert pathways['DSR']['reference_sequences'] == 11
    assert pathways['ASR']['search_status'] == 'not_searched'

    summary = results['summary']
    assert summary['databases_requested'] == 3
    assert summary['databases_searched'] == 2
    assert summary['pathways_detected'] == 1
    assert summary['complete_pathways'] == 0


# ---------------------------------------------------------- antibiotic

def test_antibiotic_rule_labels_match_documentation():
    from metharct.core.antibiotic_analyzer import AntibioticAnalyzer

    rules = {rule['name']: rule for rule in AntibioticAnalyzer.ANTIBIOTIC_RULES}

    assert ('Methanobacterium.faa', 50.0, '>') in rules['Bacitracin']['conditions']
    assert rules['Bacitracin']['logic'] == 'any'

    assert rules['Tunicamycin']['logic'] == 'all'
    assert all(
        operator == '<'
        for _, _, operator in rules['Tunicamycin']['conditions']
    )


# --------------------------------------------------------------- wsl

def test_wsl_path_conversion():
    from metharct.utils.file_utils import FileUtils

    converted = FileUtils.to_wsl_path('C:\\data\\protein.faa')
    if os.name == 'nt':
        assert converted == '/mnt/c/data/protein.faa'
    else:
        assert converted == 'C:/data/protein.faa'
