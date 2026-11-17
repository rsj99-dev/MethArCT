# -*- coding: utf-8 -*-
"""MCP tool layer: schemas, argument validation and JSON envelopes.

Keeping the tool definitions here (instead of inside ``mcp_server.py``) has two
benefits:

* the tool contracts are testable without installing the MCP SDK;
* the MCP server, the CLI and the documentation all read the same definitions.

Every handler only touches :class:`metharct.ai.api.MethArCTClient`.
"""

import uuid
from typing import Any, Callable, Dict, List, Optional

from .errors import ErrorCode, MethArCTError, to_metharct_error

#: JSON-Schema fragments shared by several tools.
_JOB_ID_SCHEMA = {"type": "string", "description": "Job identifier returned by metharct_run_analysis."}
_MODULES_SCHEMA = {
    "type": "array",
    "items": {"type": "string"},
    "description": "Module identifiers; omit to use every registered module.",
}


def _request_id() -> str:
    return "mcp_" + uuid.uuid4().hex[:12]


def envelope_ok(data: Any, warnings: Optional[List[str]] = None) -> Dict[str, Any]:
    """Build the success envelope shared by every MCP tool."""
    return {
        "ok": True,
        "request_id": _request_id(),
        "data": data,
        "warnings": list(warnings or []),
        "error": None,
    }


def envelope_error(exc: MethArCTError) -> Dict[str, Any]:
    """Build the failure envelope shared by every MCP tool."""
    return {
        "ok": False,
        "request_id": _request_id(),
        "data": None,
        "warnings": [],
        "error": exc.to_dict(),
    }


# ---------------------------------------------------------------------------
# Tool handlers
# ---------------------------------------------------------------------------

def _tool_capabilities(client, arguments: Dict[str, Any]) -> Dict[str, Any]:
    include_details = arguments.get("include_details", True)
    if not isinstance(include_details, bool):
        raise MethArCTError(
            ErrorCode.INVALID_ARGUMENT,
            "include_details must be a boolean.",
            details={"value": include_details},
        )
    return envelope_ok(client.get_capabilities(include_details=include_details))


def _tool_doctor(client, arguments: Dict[str, Any]) -> Dict[str, Any]:
    check_databases = arguments.get("check_databases", True)
    if not isinstance(check_databases, bool):
        raise MethArCTError(
            ErrorCode.INVALID_ARGUMENT,
            "check_databases must be a boolean.",
            details={"value": check_databases},
        )
    data = client.doctor(check_databases=check_databases)
    warnings = [issue.get("detail") for issue in data.get("blocking_issues", []) if issue.get("detail")]
    return envelope_ok(data, warnings=warnings)


def _tool_inspect_input(client, arguments: Dict[str, Any]) -> Dict[str, Any]:
    path = arguments.get("path")
    if not path or not isinstance(path, str):
        raise MethArCTError(
            ErrorCode.INVALID_ARGUMENT,
            "path is required and must be a string.",
        )
    expected_type = arguments.get("expected_type") or "auto"
    data = client.inspect_input(path, expected_type=expected_type)
    return envelope_ok(data, warnings=data.get("warnings"))


def _tool_run_analysis(client, arguments: Dict[str, Any]) -> Dict[str, Any]:
    input_path = arguments.get("input_path")
    if not input_path or not isinstance(input_path, str):
        raise MethArCTError(
            ErrorCode.INVALID_ARGUMENT,
            "input_path is required and must be a string.",
        )
    modules = arguments.get("modules")
    if modules is not None and (
        not isinstance(modules, list) or not all(isinstance(item, str) for item in modules)
    ):
        raise MethArCTError(
            ErrorCode.INVALID_ARGUMENT,
            "modules must be an array of strings.",
            details={"value": modules},
        )
    options = arguments.get("options") or {}
    if not isinstance(options, dict):
        raise MethArCTError(
            ErrorCode.INVALID_ARGUMENT,
            "options must be an object keyed by module name.",
            details={"value": options},
        )

    job = client.run_analysis(
        input_file=input_path,
        input_type=arguments.get("input_type") or "auto",
        modules=modules,
        options=options,
        timeout_seconds=arguments.get("timeout_seconds"),
        idempotency_key=arguments.get("idempotency_key"),
        wait=bool(arguments.get("wait", False)),
        poll_interval=0.5,
    )
    return envelope_ok(job, warnings=job.get("warnings"))


def _tool_get_job_status(client, arguments: Dict[str, Any]) -> Dict[str, Any]:
    job_id = arguments.get("job_id")
    if not job_id or not isinstance(job_id, str):
        raise MethArCTError(ErrorCode.INVALID_ARGUMENT, "job_id is required and must be a string.")
    return envelope_ok(client.get_job_status(job_id))


def _tool_get_result_summary(client, arguments: Dict[str, Any]) -> Dict[str, Any]:
    job_id = arguments.get("job_id")
    if not job_id or not isinstance(job_id, str):
        raise MethArCTError(ErrorCode.INVALID_ARGUMENT, "job_id is required and must be a string.")
    modules = arguments.get("modules")
    if modules is not None and not isinstance(modules, list):
        raise MethArCTError(
            ErrorCode.INVALID_ARGUMENT,
            "modules must be an array of module identifiers.",
            details={"value": modules},
        )
    data = client.get_result_summary(job_id, modules=modules)
    return envelope_ok(data, warnings=data.get("warnings"))


def _tool_read_artifact(client, arguments: Dict[str, Any]) -> Dict[str, Any]:
    job_id = arguments.get("job_id")
    if not job_id or not isinstance(job_id, str):
        raise MethArCTError(ErrorCode.INVALID_ARGUMENT, "job_id is required and must be a string.")

    artifact_id = arguments.get("artifact_id") or ""
    if artifact_id and not isinstance(artifact_id, str):
        raise MethArCTError(
            ErrorCode.INVALID_ARGUMENT,
            "artifact_id must be a string.",
            details={"value": artifact_id},
        )
    if not artifact_id:
        return envelope_ok({"job_id": job_id, "artifacts": client.list_artifacts(job_id)})

    offset = arguments.get("offset", 0)
    limit = arguments.get("limit", 50)
    if isinstance(offset, bool) or not isinstance(offset, int) or offset < 0:
        raise MethArCTError(
            ErrorCode.INVALID_ARGUMENT,
            "offset must be a non-negative integer.",
            details={"value": offset},
        )
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 1:
        raise MethArCTError(
            ErrorCode.INVALID_ARGUMENT,
            "limit must be a positive integer.",
            details={"value": limit},
        )

    filter_spec = arguments.get("filter") or {}
    if filter_spec and not isinstance(filter_spec, dict):
        raise MethArCTError(
            ErrorCode.INVALID_ARGUMENT,
            "filter must be an object with column/value keys.",
            details={"value": filter_spec},
        )
    column = filter_spec.get("column") if filter_spec else None
    value = filter_spec.get("value") if filter_spec else None
    contains = bool(filter_spec.get("contains", False)) if filter_spec else False
    if column is not None and not isinstance(column, str):
        raise MethArCTError(ErrorCode.INVALID_ARGUMENT, "filter.column must be a string.")
    if value is not None and not isinstance(value, str):
        raise MethArCTError(ErrorCode.INVALID_ARGUMENT, "filter.value must be a string.")

    data = client.read_artifact(
        job_id,
        artifact_id,
        offset=offset,
        limit=limit,
        filter_column=column,
        filter_value=value,
        contains=contains,
    )
    return envelope_ok(data)


def _tool_cancel_job(client, arguments: Dict[str, Any]) -> Dict[str, Any]:
    job_id = arguments.get("job_id")
    if not job_id or not isinstance(job_id, str):
        raise MethArCTError(ErrorCode.INVALID_ARGUMENT, "job_id is required and must be a string.")
    return envelope_ok(client.cancel_job(job_id))


#: The eight MCP tools, in the order used by the documentation.
TOOL_DEFINITIONS: List[Dict[str, Any]] = [
    {
        "name": "metharct_capabilities",
        "description": (
            "Discover the installed MethArCT version, the analysis modules, their "
            "options and their real dependency status. Call this first."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "include_details": {
                    "type": "boolean",
                    "default": True,
                    "description": "Include per-dependency check details.",
                }
            },
            "additionalProperties": False,
        },
        "handler": _tool_capabilities,
    },
    {
        "name": "metharct_doctor",
        "description": (
            "Check the local environment (Python packages, DIAMOND, databases, "
            "embedded models, writable directories). Never installs anything."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "check_databases": {
                    "type": "boolean",
                    "default": True,
                    "description": "Also verify database and model files.",
                }
            },
            "additionalProperties": False,
        },
        "handler": _tool_doctor,
    },
    {
        "name": "metharct_inspect_input",
        "description": (
            "Validate a FASTA file before analysis: sequence type, sequence count, "
            "residues, duplicate ids, invalid characters and compatible modules."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "path": {"type": "string", "description": "Path to the FASTA file."},
                "expected_type": {
                    "type": "string",
                    "enum": ["auto", "protein_fasta", "nucleotide_fasta"],
                    "default": "auto",
                },
            },
            "required": ["path"],
            "additionalProperties": False,
        },
        "handler": _tool_inspect_input,
    },
    {
        "name": "metharct_run_analysis",
        "description": (
            "Submit an analysis job. Returns a job_id immediately unless wait=true. "
            "Poll with metharct_get_job_status and read results with "
            "metharct_get_result_summary."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "input_path": {"type": "string", "description": "Path to the protein FASTA file."},
                "input_type": {
                    "type": "string",
                    "enum": ["auto", "protein_fasta", "nucleotide_fasta"],
                    "default": "auto",
                },
                "modules": _MODULES_SCHEMA,
                "options": {
                    "type": "object",
                    "description": "Per-module options, e.g. {\"pathway\": {\"threads\": 8}}.",
                },
                "wait": {
                    "type": "boolean",
                    "default": False,
                    "description": "Block until the job finishes (may exceed the host timeout).",
                },
                "timeout_seconds": {"type": "number", "minimum": 0},
                "idempotency_key": {"type": "string"},
            },
            "required": ["input_path"],
            "additionalProperties": False,
        },
        "handler": _tool_run_analysis,
    },
    {
        "name": "metharct_get_job_status",
        "description": "Return the current job status, phase, progress and error summary.",
        "input_schema": {
            "type": "object",
            "properties": {"job_id": _JOB_ID_SCHEMA},
            "required": ["job_id"],
            "additionalProperties": False,
        },
        "handler": _tool_get_job_status,
    },
    {
        "name": "metharct_get_result_summary",
        "description": (
            "Return the unified report summary for a job, optionally restricted to "
            "specific modules. Never returns whole result tables."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "job_id": _JOB_ID_SCHEMA,
                "modules": _MODULES_SCHEMA,
            },
            "required": ["job_id"],
            "additionalProperties": False,
        },
        "handler": _tool_get_result_summary,
    },
    {
        "name": "metharct_read_artifact",
        "description": (
            "List the artifacts of a job (omit artifact_id) or read one bounded page "
            "of a CSV/TSV/JSON artifact registered in the job manifest."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "job_id": _JOB_ID_SCHEMA,
                "artifact_id": {
                    "type": "string",
                    "description": "Artifact identifier such as 'a001'. Omit to list artifacts.",
                },
                "offset": {"type": "integer", "minimum": 0, "default": 0},
                "limit": {"type": "integer", "minimum": 1, "maximum": 1000, "default": 50},
                "filter": {
                    "type": "object",
                    "description": "Optional row filter: {column, value, contains}.",
                    "properties": {
                        "column": {"type": "string"},
                        "value": {"type": "string"},
                        "contains": {"type": "boolean", "default": False},
                    },
                },
            },
            "required": ["job_id"],
            "additionalProperties": False,
        },
        "handler": _tool_read_artifact,
    },
    {
        "name": "metharct_cancel_job",
        "description": (
            "Request cancellation of a running job. Returns whether the worker could "
            "be stopped and the resulting job status."
        ),
        "input_schema": {
            "type": "object",
            "properties": {"job_id": _JOB_ID_SCHEMA},
            "required": ["job_id"],
            "additionalProperties": False,
        },
        "handler": _tool_cancel_job,
    },
]

TOOL_NAMES = tuple(definition["name"] for definition in TOOL_DEFINITIONS)
_TOOLS_BY_NAME = {definition["name"]: definition for definition in TOOL_DEFINITIONS}


def describe_tools() -> List[Dict[str, Any]]:
    """Return the public tool manifest (name, description, input schema)."""
    return [
        {
            "name": definition["name"],
            "description": definition["description"],
            "input_schema": definition["input_schema"],
        }
        for definition in TOOL_DEFINITIONS
    ]


def call_tool(client, tool_name: str, arguments: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Execute one tool by name.

    Unknown tools and unexpected arguments are rejected explicitly; nothing is
    silently ignored.
    """
    definition = _TOOLS_BY_NAME.get(tool_name)
    if definition is None:
        return envelope_error(
            MethArCTError(
                ErrorCode.INVALID_ARGUMENT,
                "Unknown tool name.",
                details={"tool": tool_name, "available_tools": list(TOOL_NAMES)},
            )
        )

    arguments = dict(arguments or {})
    schema = definition["input_schema"]
    allowed = set((schema.get("properties") or {}).keys())
    unexpected = sorted(set(arguments) - allowed)
    if unexpected:
        return envelope_error(
            MethArCTError(
                ErrorCode.INVALID_ARGUMENT,
                "Unexpected argument(s) for tool '{0}'.".format(tool_name),
                details={"unexpected": unexpected, "allowed": sorted(allowed)},
            )
        )
    missing = [name for name in schema.get("required", []) if name not in arguments]
    if missing:
        return envelope_error(
            MethArCTError(
                ErrorCode.INVALID_ARGUMENT,
                "Missing required argument(s) for tool '{0}'.".format(tool_name),
                details={"missing": missing},
            )
        )

    handler: Callable[[Any, Dict[str, Any]], Dict[str, Any]] = definition["handler"]
    try:
        return handler(client, arguments)
    except MethArCTError as exc:
        return envelope_error(exc)
    except Exception as exc:  # noqa: BLE001 - exposed as a stable error code
        return envelope_error(to_metharct_error(exc))
