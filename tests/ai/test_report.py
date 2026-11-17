# -*- coding: utf-8 -*-
"""Tests for the unified AI report and the bundled JSON schemas."""

import json
import os

import pytest

from metharct.ai.errors import ErrorCode, MethArCTError
from metharct.ai.report import (
    CAPABILITIES_SCHEMA_FILENAME,
    INPUT_SCHEMA_FILENAME,
    JOB_SCHEMA_FILENAME,
    REPORT_SCHEMA_FILENAME,
    assert_report_valid,
    build_report,
    load_schema,
    schemas_dir,
    validate_payload,
)
from metharct.ai.schema_check import validate_document
from metharct.ai.serializers import serialize_module

SCHEMA_FILES = (
    REPORT_SCHEMA_FILENAME,
    JOB_SCHEMA_FILENAME,
    CAPABILITIES_SCHEMA_FILENAME,
    INPUT_SCHEMA_FILENAME,
)


@pytest.mark.parametrize("filename", SCHEMA_FILES)
def test_schema_files_exist_and_are_valid_json(filename):
    path = os.path.join(schemas_dir(), filename)
    assert os.path.isfile(path), path
    with open(path, "r", encoding="utf-8") as stream:
        document = json.load(stream)
    assert document["$schema"].startswith("http://json-schema.org/draft-07")
    assert document["type"] == "object"


def _report_fixture():
    job = {
        "job_id": "job_20260101T000000_abcdef12",
        "created_at": "2026-01-01T00:00:00Z",
        "started_at": "2026-01-01T00:00:01Z",
        "finished_at": "2026-01-01T00:00:02Z",
        "status": "succeeded",
        "duration_seconds": 1.0,
        "input_filename": "fixture.faa",
        "input_sha256": "b" * 64,
        "requested_modules": ["temperature"],
        "warnings": [],
        "errors": [],
    }
    input_info = {
        "filename": "fixture.faa",
        "sha256": "b" * 64,
        "format_detected": "fasta",
        "sequence_type": "protein",
        "sequence_count": 2,
        "total_residues": 120,
    }
    module_results = {
        "temperature": serialize_module(
            "temperature",
            {
                "status": "success",
                "prediction": {"T_min": 12.0, "T_opt": 30.0, "T_max": 45.0},
                "summary": {"temperature_category": "Mesophilic"},
            },
            runtime_seconds=0.5,
        )
    }
    artifacts = [
        {
            "artifact_id": "a001",
            "path": "artifacts/fixture_Huainanzi_Summary.tsv",
            "media_type": "text/tab-separated-values",
            "size_bytes": 120,
            "sha256": "c" * 64,
            "row_count": 1,
            "description": "temperature output",
            "module": "temperature",
            "filename": "fixture_Huainanzi_Summary.tsv",
        }
    ]
    return build_report(job, input_info, module_results, artifacts)


def test_report_fixture_matches_schema():
    report = _report_fixture()
    problems = validate_payload(report, REPORT_SCHEMA_FILENAME)
    assert problems == [], problems[:5]
    assert report["schema_version"] == "1.0"
    assert report["input"]["type"] == "protein"


def test_report_preserves_null_for_missing_input_stats():
    job = {"job_id": "j", "requested_modules": [], "warnings": [], "errors": []}
    info = {"filename": "x.faa", "format_detected": "fasta", "sequence_type": "protein"}
    report = build_report(job, info, {}, [])
    assert report["input"]["sequence_count"] is None
    assert report["input"]["total_residues"] is None
    assert report["input"]["sha256"] is None


def test_assert_report_valid_raises_with_problems():
    report = _report_fixture()
    report["modules"]["temperature"]["status"] = "mystery"
    with pytest.raises(MethArCTError) as excinfo:
        assert_report_valid(report)
    assert excinfo.value.code == ErrorCode.SCHEMA_VALIDATION_FAILED
    assert excinfo.value.details["problem_count"] >= 1


def test_missing_schema_file_raises(tmp_path, monkeypatch):
    with pytest.raises(MethArCTError) as excinfo:
        load_schema("does-not-exist.schema.json")
    assert excinfo.value.code == ErrorCode.SCHEMA_VALIDATION_FAILED


def test_validator_reports_missing_required_property():
    problems = validate_payload({"schema_version": "1.0"}, REPORT_SCHEMA_FILENAME)
    assert any("required" in problem for problem in problems)


def test_validator_rejects_wrong_enum_value():
    report = _report_fixture()
    report["run"]["status"] = "done"
    problems = validate_payload(report, REPORT_SCHEMA_FILENAME)
    assert any("not in" in problem for problem in problems)


def test_artifact_path_cannot_be_absolute():
    report = _report_fixture()
    report["artifacts"][0]["path"] = "/etc/passwd"
    problems = validate_payload(report, REPORT_SCHEMA_FILENAME)
    assert problems


def test_minimal_validator_supports_refs_and_types():
    document = {
        "type": "object",
        "required": ["a"],
        "properties": {"a": {"$ref": "#/definitions/positive"}},
        "definitions": {
            "positive": {"type": "integer", "minimum": 1},
        },
    }
    assert validate_document({"a": 2}, document) == []
    assert validate_document({"a": 0}, document)
    assert validate_document({}, document)
    assert validate_document({"a": "x"}, document)


def test_schemas_are_cross_checked_with_jsonschema_when_available():
    jsonschema = pytest.importorskip("jsonschema")
    report = _report_fixture()
    schema = load_schema(REPORT_SCHEMA_FILENAME)
    jsonschema.validate(report, schema)
    problems = validate_payload(report, REPORT_SCHEMA_FILENAME)
    assert problems == []
