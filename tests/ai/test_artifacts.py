# -*- coding: utf-8 -*-
"""Tests for artifact registration, path whitelisting and pagination."""

import json
import os

import pytest

from metharct.ai.artifacts import ArtifactRegistry, infer_module_from_filename
from metharct.ai.errors import ErrorCode, MethArCTError


@pytest.fixture()
def registry(tmp_path):
    job_dir = tmp_path / "job"
    (job_dir / "artifacts").mkdir(parents=True)
    return ArtifactRegistry(str(job_dir), str(job_dir / "artifacts"))


def _write(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def test_register_and_manifest(registry):
    _write(registry.artifact_dir / "genome_pH_Summary.tsv", "target\tvalue\nph_optimum\t7.1\n")
    _write(registry.artifact_dir / "genome_cultivation.csv", "a,b\n1,2\n")
    entries = registry.scan()

    assert len(entries) == 2
    assert [entry["artifact_id"] for entry in entries] == ["a001", "a002"]
    by_filename = {entry["filename"]: entry for entry in entries}
    assert by_filename["genome_pH_Summary.tsv"]["module"] == "ph"
    assert by_filename["genome_cultivation.csv"]["module"] == "cultivation"
    for entry in entries:
        assert entry["path"].startswith("artifacts/")
        assert len(entry["sha256"]) == 64
        assert entry["size_bytes"] > 0

    manifest_path = registry.write_manifest()
    with open(manifest_path, "r", encoding="utf-8") as stream:
        manifest = json.load(stream)
    assert manifest["count"] == 2
    assert manifest["artifact_root"] == "artifacts"


def test_row_count_for_tabular_artifacts(registry):
    _write(registry.artifact_dir / "genome_metabolic_summary.csv", "a,b\n1,2\n3,4\n")
    entries = registry.scan()
    assert entries[0]["row_count"] == 2


def test_register_file_outside_job_dir_is_denied(registry, tmp_path):
    outside = tmp_path / "outside.txt"
    outside.write_text("x", encoding="utf-8")
    with pytest.raises(MethArCTError) as excinfo:
        registry.register_file(str(outside))
    assert excinfo.value.code == ErrorCode.ARTIFACT_ACCESS_DENIED


def test_unknown_artifact_id_raises(registry):
    _write(registry.artifact_dir / "a.csv", "a\n1\n")
    registry.scan()
    registry.write_manifest()
    with pytest.raises(MethArCTError) as excinfo:
        registry.read("a999")
    assert excinfo.value.code == ErrorCode.ARTIFACT_NOT_FOUND


def test_read_artifact_without_manifest_raises(tmp_path):
    job_dir = tmp_path / "job"
    (job_dir / "artifacts").mkdir(parents=True)
    with pytest.raises(MethArCTError) as excinfo:
        ArtifactRegistry(str(job_dir)).read("a001")
    assert excinfo.value.code == ErrorCode.ARTIFACT_NOT_FOUND


def test_pagination_returns_next_offset(registry):
    rows = ["id\tvalue"] + ["r{0}\t{1}".format(i, i * 10) for i in range(10)]
    _write(registry.artifact_dir / "table.tsv", "\n".join(rows) + "\n")
    registry.scan()
    registry.write_manifest()

    first = registry.read("a001", offset=0, limit=4)
    assert first["total_rows"] == 10
    assert first["returned_rows"] == 4
    assert first["next_offset"] == 4
    assert first["columns"] == ["id", "value"]
    assert first["rows"][0]["id"] == "r0"

    last = registry.read("a001", offset=8, limit=4)
    assert last["returned_rows"] == 2
    assert last["next_offset"] is None


def test_filter_by_column(registry):
    rows = ["id\tstatus"] + ["a\tdetected", "b\tnot_detected", "c\tdetected"]
    _write(registry.artifact_dir / "table.tsv", "\n".join(rows) + "\n")
    registry.scan()
    registry.write_manifest()

    result = registry.read("a001", filter_column="status", filter_value="detected")
    assert result["matched_rows"] == 2
    assert [row["id"] for row in result["rows"]] == ["a", "c"]

    contains = registry.read("a001", filter_column="status", filter_value="not_det", contains=True)
    assert contains["matched_rows"] == 1


def test_unknown_filter_column_raises(registry):
    _write(registry.artifact_dir / "table.tsv", "id\tstatus\na\tdetected\n")
    registry.scan()
    registry.write_manifest()
    with pytest.raises(MethArCTError) as excinfo:
        registry.read("a001", filter_column="missing", filter_value="x")
    assert excinfo.value.code == ErrorCode.INVALID_ARGUMENT


def test_binary_artifact_is_not_read_as_text(registry):
    path = registry.artifact_dir / "result.xlsx"
    path.write_bytes(b"PK\x03\x04binary")
    registry.scan()
    registry.write_manifest()
    with pytest.raises(MethArCTError) as excinfo:
        registry.read("a001")
    assert excinfo.value.code == ErrorCode.ARTIFACT_ACCESS_DENIED


def test_json_artifact_is_paged_by_lines(registry):
    _write(registry.artifact_dir / "details.json", json.dumps({"a": 1, "b": 2}, indent=2) + "\n")
    registry.scan()
    registry.write_manifest()
    result = registry.read("a001", offset=0, limit=3)
    assert result["total_lines"] >= 3
    assert result["returned_lines"] == 3
    assert result["next_offset"] == 3


def test_manifest_path_escape_is_rejected(tmp_path):
    job_dir = tmp_path / "job"
    (job_dir / "artifacts").mkdir(parents=True)
    secret = tmp_path / "secret.txt"
    secret.write_text("do not read", encoding="utf-8")
    manifest = {
        "artifacts": [
            {
                "artifact_id": "a001",
                "path": "../secret.txt",
                "media_type": "text/plain",
                "size_bytes": 10,
                "sha256": "0" * 64,
                "row_count": None,
                "description": "escape",
                "module": None,
                "filename": "secret.txt",
            }
        ],
        "artifact_root": "artifacts",
        "count": 1,
    }
    (job_dir / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

    registry = ArtifactRegistry(str(job_dir))
    with pytest.raises(MethArCTError) as excinfo:
        registry.read("a001")
    assert excinfo.value.code == ErrorCode.ARTIFACT_ACCESS_DENIED


def test_absolute_paths_are_redacted_in_artifact_pages(tmp_path):
    """Agent-facing content must not expose host absolute paths."""
    job_dir = tmp_path / "job"
    (job_dir / "artifacts").mkdir(parents=True)
    input_file = tmp_path / "genomes" / "in.faa"
    input_file.parent.mkdir()
    input_file.write_text(">a\nMKT\n", encoding="utf-8")

    details = job_dir / "artifacts" / "details.json"
    details.write_text(
        json.dumps({"input_file": str(input_file), "workdir": str(job_dir)}, indent=2) + "\n",
        encoding="utf-8",
    )
    table = job_dir / "artifacts" / "table.tsv"
    table.write_text("source\n{0}\n".format(input_file), encoding="utf-8")

    registry = ArtifactRegistry(
        str(job_dir),
        redactions=[(str(input_file), "in.faa"), (str(job_dir), ".")],
    )
    registry.scan()
    registry.write_manifest()

    ids = {entry["filename"]: entry["artifact_id"] for entry in registry.entries}
    json_page = registry.read(ids["details.json"])
    assert str(input_file) not in json_page["content"]
    assert "in.faa" in json_page["content"]

    tsv_page = registry.read(ids["table.tsv"])
    assert str(input_file) not in str(tsv_page["rows"])
    assert tsv_page["rows"][0]["source"] == "in.faa"


def test_redaction_does_not_modify_files_on_disk(tmp_path):
    job_dir = tmp_path / "job"
    (job_dir / "artifacts").mkdir(parents=True)
    payload = "value\n{0}\n".format(tmp_path / "in.faa")
    target = job_dir / "artifacts" / "raw.txt"
    target.write_text(payload, encoding="utf-8")

    registry = ArtifactRegistry(str(job_dir), redactions=[(str(tmp_path / "in.faa"), "in.faa")])
    registry.scan()
    registry.write_manifest()
    registry.read("a001")
    assert target.read_text(encoding="utf-8") == payload


def test_infer_module_from_filename():
    assert infer_module_from_filename("x_metabolic_summary.csv") == "pathway"
    assert infer_module_from_filename("x_CO2-CH4_diamond.tsv") == "pathway"
    assert infer_module_from_filename("x_Huainanzi_Summary.tsv") == "temperature"
    assert infer_module_from_filename("x_SuSha_Summary.tsv") == "salinity"
    assert infer_module_from_filename("x_pH_Summary.tsv") == "ph"
    assert infer_module_from_filename("x_antibiotic_selection.tsv") == "antibiotic"
    assert infer_module_from_filename("x_checkm2_summary.csv") == "genome_quality"
    assert infer_module_from_filename("unknown.bin") is None


def test_registry_ignores_hidden_and_manifest_files(registry):
    _write(registry.artifact_dir / ".hidden", "x")
    _write(registry.artifact_dir / "manifest.json", "{}")
    _write(registry.artifact_dir / "visible.txt", "x")
    entries = registry.scan()
    assert [entry["filename"] for entry in entries] == ["visible.txt"]
    assert os.path.basename(entries[0]["path"]) == "visible.txt"
