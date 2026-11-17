# -*- coding: utf-8 -*-
"""Shared fixtures for the AI-layer test suite."""

import os
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if REPO_ROOT not in sys.path:  # pragma: no cover - normally already importable
    sys.path.insert(0, REPO_ROOT)

FIXTURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "fixtures")


@pytest.fixture(scope="session")
def protein_fasta() -> str:
    """A small synthetic protein FASTA fixture (not a real genome)."""
    return os.path.abspath(os.path.join(FIXTURES_DIR, "protein_small.faa"))


@pytest.fixture(scope="session")
def nucleotide_fasta() -> str:
    """A small synthetic nucleotide FASTA fixture (not a real genome)."""
    return os.path.abspath(os.path.join(FIXTURES_DIR, "nucleotide_small.fna"))


@pytest.fixture()
def client(tmp_path):
    """A MethArCT client bound to a throw-away workspace."""
    from metharct.ai import MethArCTClient

    return MethArCTClient(workspace=str(tmp_path / "workspace"), max_parallel_jobs=1)
