# -*- coding: utf-8 -*-
"""Optional MCP server exposing MethArCT to agent hosts.

The server is an *optional* entry point: the MCP Python SDK is not a hard
dependency of MethArCT. Install it with::

    pip install "metharct[mcp]"

Run it with::

    python -m metharct.ai.mcp_server --workspace ./metharct_workspace

The default transport is local ``stdio``; no network port is opened. Every tool
delegates to :class:`metharct.ai.api.MethArCTClient` and returns a compact JSON
envelope -- never a whole result table.
"""

import argparse
import json
import os
import sys
from typing import Any, Dict, List, Optional

from .constants import DEFAULT_WORKSPACE, ENV_ALLOWED_ROOTS, ENV_WORKSPACE
from .errors import ErrorCode, MethArCTError
from .mcp_tools import TOOL_DEFINITIONS, call_tool

#: Guidance handed to the agent host when the server starts.
SERVER_INSTRUCTIONS = (
    "MethArCT predicts methanogen metabolic pathways, growth temperature range, "
    "salinity adaptation, growth pH preference, antibiotic selection and genome "
    "quality from a protein FASTA file. Call metharct_capabilities first: only "
    "modules reported as available can run on this machine. Then call "
    "metharct_inspect_input, submit work with metharct_run_analysis, poll with "
    "metharct_get_job_status, and read results with metharct_get_result_summary. "
    "Predictions are model outputs, not experimental measurements. A pathway "
    "reported as not detected is not evidence that it is absent."
)

#: Static limitation text exposed as an MCP resource.
LIMITATIONS_TEXT = """# MethArCT limitations

- Temperature, salinity and pH values are model predictions derived from
  genome sequence composition. They are not experimental measurements.
- Model probabilities are model outputs and are not calibrated experimental
  confidence values.
- Pathway detection is a homology search. `not_detected` is not `absent`.
- Cultivability and amino-acid auxotrophy are inferred from biosynthetic
  pathway completeness using fixed rule thresholds.
- Antibiotic recommendations follow fixed AAI rules against five reference
  methanogen proteomes and are not a susceptibility test.
- Genome quality depends on the external CheckM2 database, which MethArCT does
  not bundle.
"""


def _missing_sdk_error() -> MethArCTError:
    return MethArCTError(
        ErrorCode.DEPENDENCY_MISSING,
        "The MCP Python SDK is not available for this interpreter.",
        details={
            "install": 'pip install "metharct[mcp]"',
            "module": "mcp",
            "python_version": "{0}.{1}.{2}".format(*sys.version_info[:3]),
            "python_requires": ">=3.10 (MCP Python SDK constraint)",
            "documentation": "https://github.com/modelcontextprotocol/python-sdk",
        },
    )


def load_server_class():
    """Return the SDK server class for both MCP SDK generations.

    * MCP SDK 2.x: ``mcp.server.mcpserver.MCPServer`` (FastMCP was renamed).
    * MCP SDK 1.x: ``mcp.server.fastmcp.FastMCP``.

    Both expose the ``tool`` / ``resource`` decorators and ``run(transport)``
    used below, so a single registration path works on either.
    """
    try:
        from mcp.server.mcpserver import MCPServer  # noqa: F401 - SDK >= 2.0

        return MCPServer
    except ImportError:
        pass
    try:
        from mcp.server.fastmcp import FastMCP  # noqa: F401 - SDK 1.x

        return FastMCP
    except ImportError:
        raise _missing_sdk_error()


def sdk_version() -> Optional[str]:
    """Return the installed MCP SDK version, or ``None`` when absent."""
    try:
        from importlib import metadata

        return metadata.version("mcp")
    except Exception:  # pragma: no cover - metadata may be unavailable
        return None


def _client(workspace: Optional[str] = None, allowed_roots: Optional[List[str]] = None, config_file: Optional[str] = None):
    from .api import MethArCTClient

    return MethArCTClient(
        workspace=workspace or os.environ.get(ENV_WORKSPACE) or DEFAULT_WORKSPACE,
        allowed_roots=allowed_roots,
        config_file=config_file,
    )


def _json(payload: Dict[str, Any]) -> str:
    return json.dumps(payload, ensure_ascii=False)


def build_server(
    workspace: Optional[str] = None,
    allowed_roots: Optional[List[str]] = None,
    config_file: Optional[str] = None,
):
    """Create the MCP server with all MethArCT tools registered.

    Raises:
        MethArCTError: when the optional ``mcp`` dependency is unavailable.
    """
    server_class = load_server_class()

    client = _client(workspace, allowed_roots, config_file)
    server = server_class("metharct", instructions=SERVER_INSTRUCTIONS)

    @server.tool(
        name="metharct_capabilities",
        description=TOOL_DEFINITIONS[0]["description"],
    )
    def metharct_capabilities(include_details: bool = True) -> str:
        return _json(call_tool(client, "metharct_capabilities", {"include_details": include_details}))

    @server.tool(
        name="metharct_doctor",
        description=TOOL_DEFINITIONS[1]["description"],
    )
    def metharct_doctor(check_databases: bool = True) -> str:
        return _json(call_tool(client, "metharct_doctor", {"check_databases": check_databases}))

    @server.tool(
        name="metharct_inspect_input",
        description=TOOL_DEFINITIONS[2]["description"],
    )
    def metharct_inspect_input(path: str, expected_type: str = "auto") -> str:
        return _json(
            call_tool(client, "metharct_inspect_input", {"path": path, "expected_type": expected_type})
        )

    @server.tool(
        name="metharct_run_analysis",
        description=TOOL_DEFINITIONS[3]["description"],
    )
    def metharct_run_analysis(
        input_path: str,
        input_type: str = "auto",
        modules: Optional[List[str]] = None,
        options: Optional[Dict[str, Any]] = None,
        wait: bool = False,
        timeout_seconds: Optional[float] = None,
        idempotency_key: Optional[str] = None,
    ) -> str:
        arguments: Dict[str, Any] = {
            "input_path": input_path,
            "input_type": input_type,
            "wait": wait,
        }
        if modules is not None:
            arguments["modules"] = modules
        if options is not None:
            arguments["options"] = options
        if timeout_seconds is not None:
            arguments["timeout_seconds"] = timeout_seconds
        if idempotency_key is not None:
            arguments["idempotency_key"] = idempotency_key
        return _json(call_tool(client, "metharct_run_analysis", arguments))

    @server.tool(
        name="metharct_get_job_status",
        description=TOOL_DEFINITIONS[4]["description"],
    )
    def metharct_get_job_status(job_id: str) -> str:
        return _json(call_tool(client, "metharct_get_job_status", {"job_id": job_id}))

    @server.tool(
        name="metharct_get_result_summary",
        description=TOOL_DEFINITIONS[5]["description"],
    )
    def metharct_get_result_summary(job_id: str, modules: Optional[List[str]] = None) -> str:
        arguments: Dict[str, Any] = {"job_id": job_id}
        if modules is not None:
            arguments["modules"] = modules
        return _json(call_tool(client, "metharct_get_result_summary", arguments))

    @server.tool(
        name="metharct_read_artifact",
        description=TOOL_DEFINITIONS[6]["description"],
    )
    def metharct_read_artifact(
        job_id: str,
        artifact_id: str = "",
        offset: int = 0,
        limit: int = 50,
        filter: Optional[Dict[str, Any]] = None,
    ) -> str:
        arguments: Dict[str, Any] = {
            "job_id": job_id,
            "artifact_id": artifact_id,
            "offset": offset,
            "limit": limit,
        }
        if filter is not None:
            arguments["filter"] = filter
        return _json(call_tool(client, "metharct_read_artifact", arguments))

    @server.tool(
        name="metharct_cancel_job",
        description=TOOL_DEFINITIONS[7]["description"],
    )
    def metharct_cancel_job(job_id: str) -> str:
        return _json(call_tool(client, "metharct_cancel_job", {"job_id": job_id}))

    @server.resource("metharct://capabilities")
    def capabilities_resource() -> str:
        """Live capability manifest (modules, dependencies, options)."""
        return _json(client.get_capabilities(include_details=True))

    @server.resource("metharct://docs/analysis-modules")
    def analysis_modules_resource() -> str:
        """Registered analysis modules and their scientific limitations."""
        from .capabilities import MODULE_REGISTRY

        lines = ["# MethArCT analysis modules", ""]
        for module in MODULE_REGISTRY:
            lines.append("## {0} (`{1}`)".format(module.display_name, module.name))
            lines.append("")
            lines.append(module.description)
            lines.append("")
            for limitation in module.limitations:
                lines.append("- {0}".format(limitation))
            lines.append("")
        return "\n".join(lines)

    @server.resource("metharct://docs/result-schema")
    def result_schema_resource() -> str:
        """The JSON schema of metharct_ai_report.json."""
        from .report import load_schema, REPORT_SCHEMA_FILENAME

        return json.dumps(load_schema(REPORT_SCHEMA_FILENAME), ensure_ascii=False, indent=2)

    @server.resource("metharct://docs/limitations")
    def limitations_resource() -> str:
        """Scientific limitations that must accompany every interpretation."""
        return LIMITATIONS_TEXT

    return server


def main(argv: Optional[List[str]] = None) -> int:
    """Entry point for ``metharct-mcp`` / ``python -m metharct.ai.mcp_server``."""
    parser = argparse.ArgumentParser(
        prog="metharct-mcp",
        description="Run the MethArCT MCP server over stdio (local, no network port).",
    )
    parser.add_argument(
        "--workspace",
        default=os.environ.get(ENV_WORKSPACE) or DEFAULT_WORKSPACE,
        help="Workspace directory holding jobs/ (default: %(default)s).",
    )
    parser.add_argument(
        "--allowed-roots",
        default=os.environ.get(ENV_ALLOWED_ROOTS) or "",
        help="Path-list separator separated directories the server may read "
        "(default: %s)." % ENV_ALLOWED_ROOTS,
    )
    parser.add_argument("--config", default=None, help="Optional MethArCT configuration file.")
    parser.add_argument(
        "--list-tools",
        action="store_true",
        help="Print the tool manifest and exit (does not start the server).",
    )
    args = parser.parse_args(argv)

    allowed_roots = [part for part in args.allowed_roots.split(os.pathsep) if part.strip()]

    if args.list_tools:
        from .mcp_tools import describe_tools

        print(json.dumps({"tools": describe_tools()}, ensure_ascii=False, indent=2))
        return 0

    try:
        server = build_server(
            workspace=args.workspace,
            allowed_roots=allowed_roots or None,
            config_file=args.config,
        )
    except MethArCTError as exc:
        print("metharct-mcp: {0}".format(exc.message), file=sys.stderr)
        if exc.details:
            print(json.dumps(exc.details, ensure_ascii=False, indent=2), file=sys.stderr)
        return 2

    try:
        server.run(transport="stdio")
    except TypeError:  # pragma: no cover - older SDK signatures
        server.run()
    return 0


if __name__ == "__main__":  # pragma: no cover - process entry point
    sys.exit(main())
