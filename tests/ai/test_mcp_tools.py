# -*- coding: utf-8 -*-
"""Tests for the MCP tool layer (no MCP SDK required)."""

import json

import pytest

from metharct.ai.errors import ErrorCode, MethArCTError
from metharct.ai.mcp_tools import TOOL_NAMES, call_tool, describe_tools


def test_eight_documented_tools_are_registered():
    assert TOOL_NAMES == (
        "metharct_capabilities",
        "metharct_doctor",
        "metharct_inspect_input",
        "metharct_run_analysis",
        "metharct_get_job_status",
        "metharct_get_result_summary",
        "metharct_read_artifact",
        "metharct_cancel_job",
    )


def test_every_tool_has_a_json_schema():
    tools = describe_tools()
    assert len(tools) == 8
    for tool in tools:
        assert tool["description"]
        schema = tool["input_schema"]
        assert schema["type"] == "object"
        assert isinstance(schema.get("properties"), dict)
        assert schema.get("additionalProperties") is False
        for required in schema.get("required", []):
            assert required in schema["properties"], tool["name"]


def test_unknown_tool_returns_validation_error(client):
    payload = call_tool(client, "metharct_drop_database", {})
    assert payload["ok"] is False
    assert payload["error"]["code"] == ErrorCode.INVALID_ARGUMENT
    assert "metharct_capabilities" in payload["error"]["details"]["available_tools"]


def test_unexpected_argument_is_rejected(client):
    payload = call_tool(client, "metharct_capabilities", {"command": "rm -rf /"})
    assert payload["ok"] is False
    assert payload["error"]["code"] == ErrorCode.INVALID_ARGUMENT
    assert payload["error"]["details"]["unexpected"] == ["command"]


def test_missing_required_argument_is_rejected(client):
    payload = call_tool(client, "metharct_inspect_input", {})
    assert payload["ok"] is False
    assert payload["error"]["details"]["missing"] == ["path"]


def test_wrong_argument_type_is_rejected(client):
    payload = call_tool(client, "metharct_capabilities", {"include_details": "yes"})
    assert payload["ok"] is False
    assert payload["error"]["code"] == ErrorCode.INVALID_ARGUMENT


def test_capabilities_tool_returns_envelope(client):
    payload = call_tool(client, "metharct_capabilities", {"include_details": False})
    assert payload["ok"] is True
    assert payload["error"] is None
    assert payload["data"]["software"]["name"] == "MethArCT"
    assert all("checks" not in module for module in payload["data"]["modules"])
    json.dumps(payload)  # must be JSON serialisable


def test_inspect_tool_reports_fixture(client, protein_fasta):
    payload = call_tool(client, "metharct_inspect_input", {"path": protein_fasta})
    assert payload["ok"] is True
    assert payload["data"]["sequence_type"] == "protein"
    assert payload["data"]["valid"] is True


def test_inspect_tool_missing_file_returns_stable_error(client, tmp_path):
    payload = call_tool(client, "metharct_inspect_input", {"path": str(tmp_path / "nope.faa")})
    assert payload["ok"] is False
    assert payload["error"]["code"] == ErrorCode.INPUT_FILE_NOT_FOUND


def test_result_tool_unknown_job_returns_stable_error(client):
    payload = call_tool(client, "metharct_get_result_summary", {"job_id": "job_20260101T000000_deadbeef"})
    assert payload["ok"] is False
    assert payload["error"]["code"] == ErrorCode.JOB_NOT_FOUND


def test_read_artifact_tool_requires_job(client):
    payload = call_tool(client, "metharct_read_artifact", {"job_id": "job_20260101T000000_deadbeef"})
    assert payload["ok"] is False
    assert payload["error"]["code"] == ErrorCode.JOB_NOT_FOUND


def test_read_artifact_tool_rejects_bad_paging(client):
    payload = call_tool(
        client,
        "metharct_read_artifact",
        {"job_id": "job_x", "artifact_id": "a001", "limit": 0},
    )
    assert payload["ok"] is False
    assert payload["error"]["code"] == ErrorCode.INVALID_ARGUMENT


def test_run_analysis_tool_requires_input_path(client):
    payload = call_tool(client, "metharct_run_analysis", {"input_path": ""})
    assert payload["ok"] is False
    assert payload["error"]["code"] == ErrorCode.INVALID_ARGUMENT


def test_mcp_server_module_imports_without_sdk():
    """The server module must be importable even when ``mcp`` is absent."""
    from metharct.ai import mcp_server

    assert callable(mcp_server.main)
    assert callable(mcp_server.build_server)


def test_build_server_reports_missing_sdk():
    """Without the optional SDK the server must fail with a clear, stable error."""
    try:
        import mcp  # noqa: F401
    except ImportError:
        from metharct.ai.mcp_server import build_server

        with pytest.raises(MethArCTError) as excinfo:
            build_server()
        assert excinfo.value.code == ErrorCode.DEPENDENCY_MISSING
        assert "mcp" in excinfo.value.details["install"]
    else:  # pragma: no cover - only when the optional SDK is installed
        from metharct.ai.mcp_server import build_server

        server = build_server(workspace=".mcp_test_workspace")
        assert server is not None
