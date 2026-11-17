# -*- coding: utf-8 -*-
"""Tests for capability discovery and the module registry."""

import os

import pytest

from metharct.ai.capabilities import (
    MODULE_NAMES,
    MODULE_REGISTRY,
    compatible_modules,
    get_capabilities,
    get_module,
    normalize_module_list,
    probe_module,
)
from metharct.ai.errors import ErrorCode, MethArCTError
from metharct.ai.report import CAPABILITIES_SCHEMA_FILENAME, validate_payload


def test_registry_covers_documented_modules():
    assert set(MODULE_NAMES) == {
        "pathway",
        "cultivation",
        "temperature",
        "salinity",
        "ph",
        "antibiotic",
        "genome_quality",
    }


def test_registry_entries_are_self_consistent():
    for module in MODULE_REGISTRY:
        assert module.input_types, module.name
        assert module.description, module.name
        assert module.output_artifacts, module.name
        for option in module.options:
            assert option.type in ("integer", "number", "string", "array[string]"), option.name


def test_unknown_module_raises():
    with pytest.raises(MethArCTError) as excinfo:
        get_module("proteomics")
    assert excinfo.value.code == ErrorCode.MODULE_UNKNOWN


def test_normalize_module_list_defaults_to_all_and_deduplicates():
    assert normalize_module_list(None) == list(MODULE_NAMES)
    assert normalize_module_list([]) == list(MODULE_NAMES)
    assert normalize_module_list(["ph", "PH", "temperature"]) == ["ph", "temperature"]


def test_normalize_module_list_rejects_typos():
    with pytest.raises(MethArCTError) as excinfo:
        normalize_module_list(["temperture"])
    assert excinfo.value.code == ErrorCode.MODULE_UNKNOWN


def test_compatible_modules_by_sequence_type():
    compatible, incompatible = compatible_modules("protein")
    assert set(compatible) == set(MODULE_NAMES)
    assert incompatible == []

    compatible, incompatible = compatible_modules("nucleotide")
    assert compatible == []
    assert set(incompatible) == set(MODULE_NAMES)


def test_capabilities_match_schema():
    payload = get_capabilities()
    problems = validate_payload(payload, CAPABILITIES_SCHEMA_FILENAME)
    assert problems == [], problems[:5]


def test_capabilities_without_details_omit_checks():
    payload = get_capabilities(include_details=False)
    assert all("checks" not in module for module in payload["modules"])


def test_probe_reports_status_from_real_checks():
    for module in MODULE_REGISTRY:
        entry = probe_module(module)
        assert entry["status"] in ("available", "unavailable", "misconfigured", "not_required")
        assert isinstance(entry["enabled"], bool)
        assert entry["enabled"] == (entry["status"] == "available")
        if not entry["enabled"]:
            assert entry["reasons"], module.name


def test_missing_checks_are_reported_for_unavailable_modules():
    """A disabled module must explain which real dependency is missing."""
    for module in MODULE_REGISTRY:
        entry = probe_module(module)
        if entry["enabled"]:
            continue
        assert entry["missing"], module.name
        for item in entry["missing"]:
            assert item["status"] in ("unavailable", "misconfigured")


def test_capability_manifest_contains_no_host_paths():
    """The capability manifest is agent-facing and must stay path-free."""
    import json

    payload = get_capabilities(include_details=True)
    text = json.dumps(payload, ensure_ascii=False)
    assert ":\\" not in text
    assert ":/" not in text.replace("metharct://", "")
    assert "/" not in payload["environment"]["python_executable"]


def test_doctor_keeps_full_paths_for_diagnosis(tmp_path):
    """Doctor must still report the real configuration it found."""
    from metharct.ai.doctor import run_doctor

    report = run_doctor(workspace=str(tmp_path), check_databases=False)
    assert report["environment"]["python_executable"]
    assert os.path.isabs(report["environment"]["python_executable"])


def test_temperature_models_are_shipped():
    """The embedded MLHuaiNanzi models must be present in this repository."""
    entry = probe_module(get_module("temperature"))
    model_checks = [check for check in entry["checks"] if check["name"].startswith("glob:")]
    assert model_checks, entry["checks"]
    assert any(check["status"] == "available" for check in model_checks)
