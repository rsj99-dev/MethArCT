# -*- coding: utf-8 -*-
"""Tests for result normalisation into the unified report shape."""

import math

import pytest

from metharct.ai.errors import ErrorCode, MethArCTError
from metharct.ai.serializers import json_safe, normalize_status, serialize_module


def test_json_safe_handles_nan_and_nested_values():
    payload = {
        "nan": float("nan"),
        "inf": float("inf"),
        "ok": 1.5,
        "nested": [{"a": float("nan")}],
        "flag": True,
    }
    result = json_safe(payload)
    assert result["nan"] is None
    assert result["inf"] is None
    assert result["ok"] == 1.5
    assert result["nested"][0]["a"] is None
    assert result["flag"] is True


def test_json_safe_converts_numpy_like_objects():
    class FakeScalar(object):
        def item(self):
            return 3

    class FakeArray(object):
        def tolist(self):
            return [1, 2, 3]

    assert json_safe(FakeScalar()) == 3
    assert json_safe(FakeArray()) == [1, 2, 3]


def test_normalize_status_mapping():
    assert normalize_status("success") == "success"
    assert normalize_status("failed") == "failed"
    assert normalize_status("skipped") == "not_available"
    assert normalize_status(None) == "success"


def test_unknown_module_cannot_be_serialized():
    with pytest.raises(MethArCTError) as excinfo:
        serialize_module("metabolomics", {"status": "success"})
    assert excinfo.value.code == ErrorCode.MODULE_UNKNOWN


def test_none_result_is_not_available():
    result = serialize_module("temperature", None)
    assert result["status"] == "not_available"
    assert result["values"] == {}
    assert result["summary"] == "The module produced no result."


def test_failed_module_keeps_error_message():
    result = serialize_module("temperature", {"status": "failed", "error": "boom"})
    assert result["status"] == "failed"
    assert result["error"] == "boom"
    assert result["values"] == {}


def test_skipped_module_becomes_not_available():
    result = serialize_module(
        "genome_quality", {"status": "skipped", "message": "CheckM2 executable is not available"}
    )
    assert result["status"] == "not_available"
    assert "CheckM2" in result["summary"]


def test_temperature_result_carries_units_and_model():
    native = {
        "status": "success",
        "prediction": {"T_min": 10.0, "T_opt": 25.0, "T_max": 40.0},
        "summary": {"temperature_category": "Mesophilic"},
    }
    result = serialize_module("temperature", native, runtime_seconds=1.25)
    assert result["status"] == "success"
    assert result["values"]["T_opt"] == 25.0
    assert result["units"]["T_opt"] == "degC"
    assert result["model"]["name"] == "MLHuaiNanzi"
    assert result["runtime_seconds"] == 1.25
    assert result["limitations"]


def test_prediction_without_values_is_partial():
    native = {
        "status": "success",
        "prediction": {
            "ph_optimum": {"value": None, "warning": "genome missing features"},
            "ph_max": {"value": None},
            "ph_min": {"value": None},
        },
    }
    result = serialize_module("ph", native)
    assert result["status"] == "partial"
    assert result["values"]["ph_optimum"] is None
    assert "genome missing features" in result["summary"]
    assert result["evidence"][0]["kind"] == "warning"


def test_pathway_detection_status_never_conflates_absent_and_undetected():
    native = {
        "status": "success",
        "summary": {"databases_requested": 3, "databases_searched": 2, "pathways_detected": 1},
        "pathway_results": {
            "CO2-CH4": {
                "pathway_name": "CO2 reduction",
                "search_status": "searched",
                "low_threshold_hits": 5,
                "reference_sequences": 15,
                "low_completeness_percentage": 33.33,
                "best_hits_low": [
                    {"qseqid": "q1", "sseqid": "K00399", "pident": 72.5, "evalue": 1e-30, "bitscore": 210.0}
                ],
            },
            "DSR": {
                "pathway_name": "Dissimilatory sulfate reduction",
                "search_status": "searched",
                "low_threshold_hits": 0,
                "reference_sequences": 11,
                "low_completeness_percentage": 0.0,
            },
            "NIT": {
                "pathway_name": "Nitrification",
                "search_status": "not_searched",
                "low_threshold_hits": 0,
                "reference_sequences": 4,
                "low_completeness_percentage": 0.0,
            },
        },
    }
    result = serialize_module("pathway", native)
    assert result["status"] == "success"
    statuses = {entry["database"]: entry["detection_status"] for entry in result["values"]["pathways"]}
    assert statuses["CO2-CH4"] == "detected"
    assert statuses["DSR"] == "not_detected"
    assert statuses["NIT"] == "not_assessed"

    assert result["evidence"][0]["reference_id"] == "K00399"
    assert any("not mean absent" in note for note in result["limitations"])


def test_salinity_result_maps_top3_evidence():
    native = {
        "status": "success",
        "prediction": {"salinity_label": "Halotolerant", "salinity_index": 1, "confidence": 0.62},
        "top3_predictions": [
            {"rank": 1, "label": "Halotolerant", "probability": 0.62},
            {"rank": 2, "label": "Salt-sensitive", "probability": 0.30},
        ],
    }
    result = serialize_module("salinity", native)
    assert result["values"]["salinity_label"] == "Halotolerant"
    assert result["units"]["confidence"] == "model_probability"
    assert len(result["evidence"]) == 2


def test_cultivation_result_surfaces_amino_acid_gaps():
    native = {
        "status": "success",
        "pathway_results": {
            "Lysine biosynthesis-DAP pathway": {
                "completeness": 0.5,
                "found_genes": 2,
                "total_genes": 4,
            }
        },
        "cultivability_assessment": {
            "overall_cultivability": "medium",
            "average_pathway_completeness": 0.55,
        },
        "amino_acid_biosynthesis": {
            "summary": {"present": 12, "missing_or_partial": 8, "missing_list": ["Lysine (Lys)"]}
        },
    }
    result = serialize_module("cultivation", native)
    assert result["values"]["overall_cultivability"] == "medium"
    assert result["values"]["amino_acid_summary"]["missing_list"] == ["Lysine (Lys)"]
    assert result["values"]["pathways"][0]["detection_status"] == "not_detected"


def test_source_artifacts_use_manifest_ids():
    native = {"status": "success", "output_files": {"tsv": "/tmp/job/artifacts/x_pH_Summary.tsv"}}
    result = serialize_module("ph", native, artifact_resolver=lambda path: "a004")
    assert result["source_artifacts"] == [{"artifact_id": "a004", "kind": "tsv"}]


def test_serializers_do_not_invent_values():
    native = {"status": "success", "prediction": {}}
    result = serialize_module("temperature", native)
    assert result["values"]["T_min"] is None
    assert result["values"]["T_opt"] is None
    assert result["values"]["T_max"] is None
    assert not any(isinstance(value, float) and math.isnan(value) for value in result["values"].values())
