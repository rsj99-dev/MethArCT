# -*- coding: utf-8 -*-
"""Tests for FASTA inspection and sequence type detection."""

import os

import pytest

from metharct.ai.errors import ErrorCode, MethArCTError
from metharct.ai.validators import assert_sequence_type_supported, inspect_fasta


def test_inspects_protein_fixture(protein_fasta):
    info = inspect_fasta(protein_fasta)
    assert info["file_exists"] is True
    assert info["format_detected"] == "fasta"
    assert info["sequence_type"] == "protein"
    assert info["sequence_type_confidence"] == "high"
    assert info["sequence_count"] == 3
    assert info["total_residues"] > 300
    assert info["length"]["min"] <= info["length"]["median"] <= info["length"]["max"]
    assert info["invalid_sequence_count"] == 0
    assert info["duplicate_id_count"] == 0
    assert info["valid"] is True
    assert len(info["sha256"]) == 64


def test_inspects_nucleotide_fixture(nucleotide_fasta):
    info = inspect_fasta(nucleotide_fasta)
    assert info["sequence_type"] == "nucleotide"
    assert info["sequence_count"] == 2
    assert info["valid"] is True


def test_missing_file_raises_input_file_not_found(tmp_path):
    with pytest.raises(MethArCTError) as excinfo:
        inspect_fasta(str(tmp_path / "missing.faa"))
    assert excinfo.value.code == ErrorCode.INPUT_FILE_NOT_FOUND


def test_empty_file_raises(tmp_path):
    empty = tmp_path / "empty.faa"
    empty.write_text("", encoding="utf-8")
    with pytest.raises(MethArCTError) as excinfo:
        inspect_fasta(str(empty))
    assert excinfo.value.code == ErrorCode.INPUT_FORMAT_INVALID


def test_file_without_header_raises(tmp_path):
    bad = tmp_path / "no_header.faa"
    bad.write_text("MKVLLAG\nMADKIL\n", encoding="utf-8")
    with pytest.raises(MethArCTError) as excinfo:
        inspect_fasta(str(bad))
    assert excinfo.value.code == ErrorCode.INPUT_FORMAT_INVALID


def test_duplicate_ids_and_invalid_characters_are_reported(tmp_path):
    path = tmp_path / "dirty.faa"
    path.write_text(
        ">dup protein one\nMKVLLAGTHAGK\n"
        ">dup protein two\nMADKILDFGGAW\n"
        ">third protein\nMKVLL\x01AGTHAGK\n",
        encoding="utf-8",
    )
    info = inspect_fasta(str(path))
    assert info["duplicate_id_count"] == 1
    assert info["invalid_sequence_count"] == 1
    assert info["valid"] is False
    assert any("duplicate" in warning for warning in info["warnings"])
    assert any("outside the protein" in error for error in info["errors"])


def test_oversized_file_is_rejected(tmp_path):
    path = tmp_path / "big.faa"
    path.write_text(">a\n" + "M" * 5000 + "\n", encoding="utf-8")
    with pytest.raises(MethArCTError) as excinfo:
        inspect_fasta(str(path), max_file_size_mb=0.001)
    assert excinfo.value.code == ErrorCode.INPUT_FILE_TOO_LARGE


def test_expected_type_mismatch_is_invalid(protein_fasta):
    info = inspect_fasta(protein_fasta, expected_type="nucleotide_fasta")
    assert info["expected_type"] == "nucleotide"
    assert info["valid"] is False
    assert any("Expected nucleotide" in error for error in info["errors"])


def test_unsupported_expected_type_raises(protein_fasta):
    with pytest.raises(MethArCTError) as excinfo:
        inspect_fasta(protein_fasta, expected_type="rna")
    assert excinfo.value.code == ErrorCode.INVALID_ARGUMENT


def test_compressed_input_is_rejected(tmp_path):
    path = tmp_path / "genome.faa.gz"
    path.write_bytes(b"\x1f\x8b\x08\x00" + b"\x00" * 32)
    with pytest.raises(MethArCTError) as excinfo:
        inspect_fasta(str(path))
    assert excinfo.value.code == ErrorCode.INPUT_FORMAT_INVALID


def test_allowed_roots_are_enforced(tmp_path, protein_fasta):
    with pytest.raises(MethArCTError) as excinfo:
        inspect_fasta(protein_fasta, allowed_roots=[str(tmp_path)])
    assert excinfo.value.code in (ErrorCode.INPUT_FILE_UNREADABLE, ErrorCode.ARTIFACT_ACCESS_DENIED)


def test_sequence_type_compatibility_helper(protein_fasta):
    info = inspect_fasta(protein_fasta)
    assert_sequence_type_supported(info, "pathway", ("protein",))
    with pytest.raises(MethArCTError) as excinfo:
        assert_sequence_type_supported(info, "pathway", ("nucleotide",))
    assert excinfo.value.code == ErrorCode.INPUT_MODULE_INCOMPATIBLE


def test_sha256_matches_hashlib(protein_fasta):
    import hashlib

    expected = hashlib.sha256(open(protein_fasta, "rb").read()).hexdigest()
    assert inspect_fasta(protein_fasta)["sha256"] == expected


def test_inspection_does_not_modify_input(protein_fasta):
    before = os.path.getsize(protein_fasta)
    inspect_fasta(protein_fasta)
    assert os.path.getsize(protein_fasta) == before
