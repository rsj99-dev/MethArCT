# -*- coding: utf-8 -*-
"""End-to-end MCP test: a real client session against ``metharct-mcp``.

The MCP Python SDK requires Python >= 3.10, while MethArCT itself supports
3.8+. This module therefore skips (with an explicit reason) on older
interpreters instead of silently passing.

Run it in a 3.10+ environment::

    pip install ".[mcp]" pytest
    python -m pytest tests/ai/test_mcp_integration.py -q
"""

import asyncio
import json
import os
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _require_sdk():
    if sys.version_info < (3, 10):
        pytest.skip(
            "The MCP Python SDK requires Python >= 3.10; this interpreter is "
            "{0}.{1}. Run this test on a newer interpreter.".format(*sys.version_info[:2])
        )
    try:
        import mcp  # noqa: F401
    except ImportError:
        pytest.skip('MCP SDK not installed; run: pip install "metharct[mcp]"')
    from mcp import ClientSession, StdioServerParameters

    return ClientSession, StdioServerParameters


def _server_params(workspace, allowed_roots):
    _client_session, StdioServerParameters = _require_sdk()
    env = dict(os.environ)
    existing = env.get("PYTHONPATH")
    env["PYTHONPATH"] = REPO_ROOT + os.pathsep + existing if existing else REPO_ROOT
    env["PYTHONUNBUFFERED"] = "1"
    env["PYTHONIOENCODING"] = "utf-8"

    args = ["-m", "metharct.ai.mcp_server", "--workspace", str(workspace)]
    if allowed_roots:
        args.extend(["--allowed-roots", os.pathsep.join(allowed_roots)])
    return StdioServerParameters(
        command=sys.executable,
        args=args,
        env=env,
        cwd=REPO_ROOT,
    )


def _tool_schema(tool):
    """Return a tool's input schema across MCP SDK 1.x (inputSchema) and 2.x."""
    schema = getattr(tool, "input_schema", None)
    if schema is None:
        schema = getattr(tool, "inputSchema", None)
    return schema or {}


def _payload(result):
    """Decode a ``CallToolResult`` payload produced by the MethArCT tools."""
    assert result.is_error is False, getattr(result, "content", None)
    texts = [
        block.text
        for block in result.content
        if getattr(block, "type", None) == "text" and getattr(block, "text", None)
    ]
    assert texts, "tool returned no text content"
    return json.loads("".join(texts))


async def _open_session(workspace, allowed_roots, body):
    ClientSession, _params = _require_sdk()
    from mcp.client.stdio import stdio_client

    params = _server_params(workspace, allowed_roots)
    async with stdio_client(params) as (read_stream, write_stream):
        async with ClientSession(read_stream, write_stream) as session:
            init = await session.initialize()
            return await body(session, init)


def test_mcp_session_lists_tools_and_resources(tmp_path):
    """A real MCP client must see the eight documented tools and four resources."""
    from metharct.ai.mcp_tools import TOOL_NAMES

    workspace = tmp_path / "workspace"

    async def body(session, init):
        server_name = getattr(getattr(init, "server_info", None), "name", "") or ""
        tools = await session.list_tools()
        resources = await session.list_resources()
        limitations = await session.read_resource("metharct://docs/limitations")
        schema = await session.read_resource("metharct://docs/result-schema")
        return server_name, tools, resources, limitations, schema

    server_name, tools, resources, limitations, schema = asyncio.run(
        _open_session(str(workspace), [REPO_ROOT], body)
    )

    assert "metharct" in server_name.lower()
    names = [tool.name for tool in tools.tools]
    assert names == list(TOOL_NAMES), names
    for tool in tools.tools:
        assert tool.description
        assert _tool_schema(tool)["type"] == "object"

    uris = sorted(str(resource.uri) for resource in resources.resources)
    assert uris == sorted(
        [
            "metharct://capabilities",
            "metharct://docs/analysis-modules",
            "metharct://docs/result-schema",
            "metharct://docs/limitations",
        ]
    )

    limit_text = limitations.contents[0].text
    assert "not evidence" in limit_text or "not `absent`" in limit_text
    schema_payload = json.loads(schema.contents[0].text)
    assert schema_payload["title"]
    assert schema_payload["properties"]["schema_version"]["const"] == "1.0"


def test_mcp_session_runs_a_real_analysis(tmp_path, protein_fasta):
    """Submit a job through MCP, poll it and read the structured result."""

    async def body(session, init):
        capabilities = _payload(await session.call_tool("metharct_capabilities", {"include_details": False}))
        assert capabilities["ok"] is True
        enabled = [module["name"] for module in capabilities["data"]["modules"] if module["enabled"]]
        if not enabled:
            return {"skip": "no analysis module is available in this interpreter"}

        module = "temperature" if "temperature" in enabled else enabled[0]

        inspection = _payload(
            await session.call_tool("metharct_inspect_input", {"path": protein_fasta})
        )
        assert inspection["ok"] is True
        assert inspection["data"]["sequence_type"] == "protein"

        submitted = _payload(
            await session.call_tool(
                "metharct_run_analysis",
                {
                    "input_path": protein_fasta,
                    "modules": [module],
                    "wait": True,
                },
            )
        )
        job_id = submitted["data"]["job_id"]

        status = _payload(await session.call_tool("metharct_get_job_status", {"job_id": job_id}))
        summary = _payload(
            await session.call_tool(
                "metharct_get_result_summary", {"job_id": job_id, "modules": [module]}
            )
        )
        listing = _payload(await session.call_tool("metharct_read_artifact", {"job_id": job_id}))
        artifacts = listing["data"]["artifacts"]
        page = None
        if artifacts:
            tabular = next(
                (
                    item
                    for item in artifacts
                    if item["media_type"] in ("text/tab-separated-values", "text/csv")
                ),
                artifacts[0],
            )
            page = _payload(
                await session.call_tool(
                    "metharct_read_artifact",
                    {"job_id": job_id, "artifact_id": tabular["artifact_id"], "limit": 5},
                )
            )

        cancelled = _payload(await session.call_tool("metharct_cancel_job", {"job_id": job_id}))
        missing = _payload(
            await session.call_tool(
                "metharct_get_job_status", {"job_id": "job_20260101T000000_deadbeef"}
            )
        )
        return {
            "module": module,
            "submitted": submitted,
            "status": status,
            "summary": summary,
            "artifacts": artifacts,
            "page": page,
            "cancelled": cancelled,
            "missing": missing,
        }

    workspace = tmp_path / "workspace"
    outcome = asyncio.run(_open_session(str(workspace), [REPO_ROOT], body))
    if "skip" in outcome:
        pytest.skip(outcome["skip"])

    submitted = outcome["submitted"]
    assert submitted["ok"] is True
    assert submitted["data"]["job_id"].startswith("job_")
    assert submitted["data"]["status"] in ("succeeded", "partially_succeeded"), submitted["data"]

    status = outcome["status"]
    assert status["ok"] is True
    assert status["data"]["status"] in ("succeeded", "partially_succeeded")
    assert status["data"]["started_at"] and status["data"]["finished_at"]

    summary = outcome["summary"]
    assert summary["ok"] is True
    assert summary["data"]["report_available"] is True
    module_result = summary["data"]["modules"][outcome["module"]]
    assert module_result["status"] in ("success", "partial")
    assert isinstance(module_result["values"], dict)
    assert module_result["limitations"]

    artifacts = outcome["artifacts"]
    assert artifacts, "a completed job must register at least one artifact"
    for entry in artifacts:
        assert entry["path"].startswith("artifacts/")

    page = outcome["page"]
    assert page["ok"] is True
    assert page["data"]["job_id"]

    cancelled = outcome["cancelled"]
    assert cancelled["ok"] is True
    assert cancelled["data"]["accepted"] is False  # the job already finished
    assert cancelled["data"]["cancellation_supported"] is True

    missing = outcome["missing"]
    assert missing["ok"] is False
    assert missing["error"]["code"] == "JOB_NOT_FOUND"


def test_mcp_session_rejects_inputs_outside_allowed_roots(tmp_path):
    """The MCP server must honour ``--allowed-roots`` for input paths."""
    outside = tmp_path / "outside.faa"
    outside.write_text(">a\nMKVLLAGTHAGK\n", encoding="utf-8")

    async def body(session, init):
        return _payload(
            await session.call_tool("metharct_inspect_input", {"path": str(outside)})
        )

    result = asyncio.run(_open_session(str(tmp_path / "workspace"), [REPO_ROOT], body))
    assert result["ok"] is False
    assert result["error"]["code"] in ("INPUT_FILE_UNREADABLE", "ARTIFACT_ACCESS_DENIED")


def test_mcp_transport_drops_undeclared_arguments(tmp_path):
    """Fields such as ``command`` can never reach MethArCT over the MCP transport.

    The MCP SDK validates each tool against its input schema before dispatch, so
    an undeclared field is dropped and the tool body never sees it. The Python
    tool layer additionally rejects such fields with ``INVALID_ARGUMENT``
    (covered by ``tests/ai/test_mcp_tools.py``), which keeps the guarantee for
    direct/in-process callers as well.
    """

    async def body(session, init):
        declared = _payload(
            await session.call_tool("metharct_capabilities", {"include_details": False})
        )
        injected = _payload(
            await session.call_tool(
                "metharct_capabilities",
                {"include_details": False, "command": "rm -rf /"},
            )
        )
        return declared, injected

    declared, injected = asyncio.run(
        _open_session(str(tmp_path / "workspace"), [REPO_ROOT], body)
    )

    assert declared["ok"] is True
    assert injected["ok"] is True, injected.get("error")
    # The injected field must have no effect at all: same software block, no echo.
    assert injected["data"]["software"] == declared["data"]["software"]
    assert "rm -rf" not in json.dumps(injected, ensure_ascii=False)
