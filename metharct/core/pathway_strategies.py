#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pathway evaluation strategies for MethArCT DiamondAnalyzer.

Each entry in ``PATHWAY_STRATEGIES`` describes *how* to evaluate a specific
pathway.  The ``evaluate_pathway`` dispatcher reads the strategy and returns
the matched gene count, eliminating hundreds of lines of ``if/elif`` in the
main analyser.
"""

from typing import Dict, Any, List, Optional
import pandas as pd


# ------------------------------------------------------------------
# Helper functions used inside strategy definitions
# ------------------------------------------------------------------

def _unique_sseqid(hits: pd.DataFrame) -> set:
    """Return unique subject-sequence IDs (column index 1)."""
    if hits.empty:
        return set()
    return set(hits.iloc[:, 1].unique())


def _strict_hits(hits: pd.DataFrame, gene: str,
                 min_bitscore: float, max_evalue: float) -> pd.DataFrame:
    """Filter *hits* for a single gene with strict quality thresholds."""
    mask = (
        (hits.iloc[:, 1] == gene)
        & (hits['bitscore'] > min_bitscore)
        & (hits['evalue'] <= max_evalue)
    )
    return hits[mask]


# ------------------------------------------------------------------
# Strategy definition
# ------------------------------------------------------------------
#
# Each strategy is a dict with:
#   key_genes      – list of gene names that must ALL be present (simple match)
#   strict_genes   – list of dicts {gene, min_bitscore, max_evalue} for strict matching
#   weighted_genes – list of dicts {genes: [...], weight: float} for weighted scoring
#   extra_logic    – optional callable(hits, ref_count, strategy) -> int | None
#                    that can override the matched count directly
#
# If *all* key/strict/weighted conditions are met the pathway is considered
# 100 % complete (matched = ref_count).  Otherwise the actual unique hit
# count is used (optionally adjusted by ``extra_logic``).
# ------------------------------------------------------------------

PATHWAY_STRATEGIES: Dict[str, Dict[str, Any]] = {
    # ── JIASUAN-CH4 (Formate methanogenesis) ──────────────────────
    'JIASUAN-CH4': {
        'weighted_genes': [
            {'genes': ['FwdF', 'FwdG', 'FwdH'], 'weight': 15},
            {'genes': ['Hmd', 'Mtd'], 'weight': 15},
        ],
        'strict_genes': [
            {'gene': 'FdhA', 'min_bitscore': 380, 'max_evalue': 1e-100, 'weight': 15},
            {'gene': 'FdhB', 'min_bitscore': 380, 'max_evalue': 1e-100, 'weight': 15},
        ],
        'base_weight': 40,
        'complete_score': 80,
    },
    # ── CO2-CH4 (CO2 reduction methanogenesis) ───────────────────
    'CO2-CH4': {
        'weighted_genes': [
            {'genes': ['FwdF', 'FwdG', 'FwdH'], 'weight': 20},
            {'genes': ['Hmd', 'Mtd'], 'weight': 20},
        ],
        'strict_genes': [],
        'base_weight': 60,
        'complete_score': 80,
    },
    # ── JIALIUCHUN-CH4 (Methanethiol methanogenesis) ────────────
    'JIALIUCHUN-CH4': {
        'strict_genes': [
            {'gene': 'MtsA1', 'min_bitscore': 200, 'max_evalue': 1e-100, 'weight': 0},
            {'gene': 'MtsA2', 'min_bitscore': 200, 'max_evalue': 1e-100, 'weight': 0},
        ],
        'mcr_genes': ['KYC55281.1', 'KYC55283.1', 'KYC55284.1', 'KYC55314.1'],
        'mcr_min_count': 3,
    },
    # ── Glycine betaine methanogenesis ───────────────────────────
    'Glycine betaine methanogenesis': {
        'key_genes': ['MtgB', 'dimethylamine_corrinoid_protein_3', 'MV10360'],
    },
    # ── Methylthiopropionate methanogenesis ──────────────────────
    'Methylthiopropionate methanogenesis': {
        'strict_genes': [
            {'gene': 'mtpA1', 'min_bitscore': 100, 'max_evalue': 1e-5, 'weight': 0},
            {'gene': 'mtsA1', 'min_bitscore': 100, 'max_evalue': 1e-5, 'weight': 0},
            {'gene': 'mtpA2', 'min_bitscore': 100, 'max_evalue': 1e-5, 'weight': 0},
            {'gene': 'mtsA2', 'min_bitscore': 100, 'max_evalue': 1e-5, 'weight': 0},
        ],
    },
    # ── Tetramethylammonium methanogenesis ───────────────────────
    'Tetramethylammonium methanogenesis': {
        'strict_genes': [
            {'gene': 'MtqA/MT2', 'min_bitscore': 200, 'max_evalue': 1e-100, 'weight': 0},
            {'gene': 'MtqB', 'min_bitscore': 200, 'max_evalue': 1e-100, 'weight': 0},
            {'gene': 'MtqC', 'min_bitscore': 200, 'max_evalue': 1e-100, 'weight': 0},
        ],
    },
    # ── Methanol dismutation methanogenesis ──────────────────────
    'Methanol dismutation methanogenesis': {
        'key_genes_any': ['MvhA', 'elpA'],  # at least one must be present
        'other_required': [
            'FwdA', 'FwdB', 'FwdC', 'FwdD', 'FwdE', 'FwdF', 'FwdG', 'FwdH',
            'Ftr', 'Mch', 'Hmd', 'Mtd', 'elpB', 'elpC',
        ],
    },
}


# ------------------------------------------------------------------
# Public dispatcher
# ------------------------------------------------------------------

def evaluate_pathway(
    db_name: str,
    low_threshold_hits: pd.DataFrame,
    reference_count: int,
    logger=None,
) -> int:
    """Evaluate a pathway and return the *matched gene count*.

    Parameters
    ----------
    db_name : str
        Pathway database key (e.g. ``'JIASUAN-CH4'``).
    low_threshold_hits : pd.DataFrame
        Diamond hits that passed the low-threshold filter.
    reference_count : int
        Number of reference sequences for this pathway.
    logger : logging.Logger, optional
        Logger for debug output.

    Returns
    -------
    int
        Number of matched genes (capped at *reference_count*).
    """
    strategy = PATHWAY_STRATEGIES.get(db_name)
    if strategy is None:
        # Default: count unique subject IDs
        return len(low_threshold_hits.drop_duplicates(subset=['sseqid']))

    unique_hits = _unique_sseqid(low_threshold_hits)

    # ── JIASUAN-CH4 / CO2-CH4: weighted scoring ─────────────────
    if 'base_weight' in strategy:
        return _evaluate_weighted(
            db_name, low_threshold_hits, unique_hits,
            reference_count, strategy, logger,
        )

    # ── JIALIUCHUN-CH4: MtsA + MCR genes ────────────────────────
    if 'mcr_genes' in strategy:
        return _evaluate_mcr(
            db_name, low_threshold_hits, unique_hits,
            reference_count, strategy, logger,
        )

    # ── Glycine betaine: simple key-gene presence ────────────────
    if 'key_genes' in strategy:
        all_present = all(g in unique_hits for g in strategy['key_genes'])
        if all_present:
            _debug(logger, f"{db_name}: All key genes satisfied")
            return reference_count
        _debug(logger, f"{db_name}: Key genes not met, using actual count {len(unique_hits)}")
        return len(unique_hits)

    # ── Methylthiopropionate / Tetramethylammonium: strict genes ─
    if 'strict_genes' in strategy and not strategy.get('weighted_genes'):
        all_present = all(
            len(_strict_hits(low_threshold_hits, sg['gene'],
                             sg['min_bitscore'], sg['max_evalue'])) > 0
            for sg in strategy['strict_genes']
        )
        if all_present:
            _debug(logger, f"{db_name}: All strict gene conditions satisfied")
            return reference_count
        _debug(logger, f"{db_name}: Strict gene conditions not met, using actual count {len(unique_hits)}")
        return len(unique_hits)

    # ── Methanol dismutation: other_required + any of key_genes_any ─
    if 'other_required' in strategy:
        other_count = sum(1 for g in strategy['other_required'] if g in unique_hits)
        all_other = other_count == len(strategy['other_required'])
        any_key = any(g in unique_hits for g in strategy.get('key_genes_any', []))
        if all_other and any_key:
            _debug(logger, f"{db_name}: All conditions satisfied")
            return reference_count
        actual = len(unique_hits)
        _debug(logger, f"{db_name}: Conditions not met, using actual count {actual}")
        return min(actual, reference_count)

    # Fallback
    return len(low_threshold_hits.drop_duplicates(subset=['sseqid']))


# ------------------------------------------------------------------
# Private helpers
# ------------------------------------------------------------------

def _evaluate_weighted(db_name, hits, unique_hits, ref_count, strategy, logger):
    """Weighted scoring for JIASUAN-CH4 and CO2-CH4."""
    base_weight = strategy['base_weight']
    complete_score = strategy.get('complete_score', 80)

    base_score = min(len(unique_hits) / ref_count * base_weight, base_weight)

    weighted_score = 0.0
    for wg in strategy.get('weighted_genes', []):
        if any(g in unique_hits for g in wg['genes']):
            weighted_score += wg['weight']

    strict_score = 0.0
    all_strict_met = True
    for sg in strategy.get('strict_genes', []):
        if len(_strict_hits(hits, sg['gene'], sg['min_bitscore'], sg['max_evalue'])) > 0:
            strict_score += sg.get('weight', 0)
        else:
            all_strict_met = False

    total_score = base_score + weighted_score + strict_score
    _debug(logger, f"{db_name}: total_score={total_score:.1f} (base={base_score:.1f}, "
                   f"weighted={weighted_score:.1f}, strict={strict_score:.1f})")

    # Check whether ALL key conditions are satisfied
    all_weighted_met = all(
        any(g in unique_hits for g in wg['genes'])
        for wg in strategy.get('weighted_genes', [])
    )
    all_conditions = all_weighted_met and all_strict_met

    if all_conditions:
        _debug(logger, f"{db_name}: All key conditions satisfied -> ref_count={ref_count}")
        return ref_count
    if total_score >= complete_score:
        _debug(logger, f"{db_name}: Score {total_score} >= {complete_score} -> ref_count={ref_count}")
        return ref_count
    _debug(logger, f"{db_name}: Score {total_score} < {complete_score} -> actual {len(unique_hits)}")
    return len(unique_hits)


def _evaluate_mcr(db_name, hits, unique_hits, ref_count, strategy, logger):
    """MCR-based evaluation for JIALIUCHUN-CH4."""
    # Strict gene checks
    strict_present = {}
    for sg in strategy.get('strict_genes', []):
        strict_present[sg['gene']] = (
            len(_strict_hits(hits, sg['gene'],
                             sg['min_bitscore'], sg['max_evalue'])) > 0
        )

    mcr_genes = strategy['mcr_genes']
    mcr_count = sum(1 for g in mcr_genes if g in unique_hits)
    mcr_ok = mcr_count >= strategy['mcr_min_count']

    _debug(logger, f"{db_name}: strict={strict_present}, mcr={mcr_count}/{len(mcr_genes)}")

    # MtsA1 is the essential strict gene
    mtsA1_ok = strict_present.get('MtsA1', False)
    mtsA2_ok = strict_present.get('MtsA2', False)

    if mtsA1_ok and mtsA2_ok and mcr_ok:
        _debug(logger, f"{db_name}: All conditions satisfied -> ref_count={ref_count}")
        return ref_count

    adjusted = len(unique_hits)
    if not mtsA1_ok:
        adjusted = max(0, adjusted - 1)
        _debug(logger, f"{db_name}: MtsA1 not strict-qualified, adjusted to {adjusted}")
    else:
        _debug(logger, f"{db_name}: Conditions not fully met, using actual {len(unique_hits)}")
    return adjusted


def _debug(logger, msg: str):
    """Log debug message if logger is available."""
    if logger is not None:
        logger.debug(msg)
