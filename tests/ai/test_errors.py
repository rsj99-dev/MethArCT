# -*- coding: utf-8 -*-
"""Tests for stable error codes and envelope formatting."""

import subprocess

from metharct.ai.errors import (
    ALL_ERROR_CODES,
    ErrorCode,
    MethArCTError,
    error,
    to_metharct_error,
)

#: Error codes promised by docs/AI_REPORT_SCHEMA.md and docs/MCP_TOOLS.md.
CONTRACT_CODES = {
    "INPUT_FILE_NOT_FOUND",
    "INPUT_FILE_UNREADABLE",
    "INPUT_FORMAT_INVALID",
    "INPUT_TYPE_AMBIGUOUS",
    "INPUT_MODULE_INCOMPATIBLE",
    "MODULE_UNKNOWN",
    "MODULE_UNAVAILABLE",
    "DEPENDENCY_MISSING",
    "DATABASE_MISSING",
    "DATABASE_INVALID",
    "ANALYSIS_TIMEOUT",
    "ANALYSIS_FAILED",
    "JOB_NOT_FOUND",
    "ARTIFACT_NOT_FOUND",
    "ARTIFACT_ACCESS_DENIED",
    "SCHEMA_VALIDATION_FAILED",
    "INTERNAL_ERROR",
}


def test_contract_error_codes_are_present():
    assert CONTRACT_CODES.issubset(set(ALL_ERROR_CODES))


def test_error_envelope_shape():
    exc = error(ErrorCode.INPUT_FILE_NOT_FOUND, "missing", details={"path": "x"})
    payload = exc.to_dict()
    assert set(payload) == {"code", "message", "retryable", "details"}
    assert payload["code"] == "INPUT_FILE_NOT_FOUND"
    assert payload["retryable"] is False
    assert payload["details"] == {"path": "x"}
    assert str(exc).startswith("INPUT_FILE_NOT_FOUND:")


def test_mapping_of_common_exceptions():
    assert to_metharct_error(FileNotFoundError("x")).code == ErrorCode.INPUT_FILE_NOT_FOUND
    assert to_metharct_error(PermissionError("x")).code == ErrorCode.INPUT_FILE_UNREADABLE
    assert to_metharct_error(ValueError("x")).code == ErrorCode.INVALID_ARGUMENT
    timeout = subprocess.TimeoutExpired(cmd="diamond", timeout=1)
    mapped = to_metharct_error(timeout)
    assert mapped.code == ErrorCode.ANALYSIS_TIMEOUT
    assert mapped.retryable is True


def test_mapping_preserves_original_exception_details():
    mapped = to_metharct_error(RuntimeError("boom"))
    assert mapped.code == ErrorCode.INTERNAL_ERROR
    assert mapped.details["exception_type"] == "RuntimeError"
    assert mapped.details["exception_message"] == "boom"


def test_mapping_is_idempotent():
    original = MethArCTError(ErrorCode.JOB_NOT_FOUND, "nope")
    assert to_metharct_error(original) is original
