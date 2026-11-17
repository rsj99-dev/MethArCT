# -*- coding: utf-8 -*-
"""Machine-readable CLI commands for the MethArCT AI layer.

These subcommands are additive: the legacy ``comprehensive`` / ``diamond`` /
``susha`` / ... commands keep their original arguments and text output.

In ``--format json`` mode stdout carries exactly one JSON object and every
progress or diagnostic message goes to stderr.
"""

import argparse
import json
import sys
import uuid
from typing import Any, Dict, List, Optional

from ..ai.constants import DEFAULT_WORKSPACE
from ..ai.errors import ErrorCode, MethArCTError, to_metharct_error
from ..ai.paths import atomic_write_json

#: Subcommands implemented by this module (used by ``metharct.cli.main``).
AI_COMMANDS = frozenset(
    ["capabilities", "doctor", "inspect", "run", "result", "artifact", "cancel", "jobs"]
)

#: Exit codes documented in ``docs/AI_API.md``.
EXIT_OK = 0
EXIT_ERROR = 1
EXIT_USAGE = 2
EXIT_JOB_FAILED = 3


def add_ai_subparsers(subparsers) -> None:
    """Register every AI subcommand on the argparse *subparsers* action."""

    def _common(parser: argparse.ArgumentParser, workspace: bool = True) -> None:
        parser.add_argument(
            "--format",
            choices=("text", "json"),
            default="text",
            help="Output format (default: text).",
        )
        parser.add_argument(
            "--json-out",
            default=None,
            help="Write the JSON payload atomically to this path.",
        )
        if workspace:
            parser.add_argument(
                "--workspace",
                default=DEFAULT_WORKSPACE,
                help="Workspace directory holding jobs/ (default: %(default)s).",
            )

    capabilities = subparsers.add_parser(
        "capabilities",
        help="Show dynamically probed modules, dependencies and options",
        description="Report MethArCT version, available modules and their real dependency status.",
    )
    _common(capabilities)
    capabilities.add_argument(
        "--no-details",
        action="store_true",
        help="Omit the per-check dependency details.",
    )

    doctor = subparsers.add_parser(
        "doctor",
        help="Check the local environment (never installs anything)",
        description="Verify Python, DIAMOND, databases, embedded models and writable directories.",
    )
    _common(doctor)
    doctor.add_argument(
        "--skip-databases",
        action="store_true",
        help="Skip database and model file checks.",
    )

    inspect = subparsers.add_parser(
        "inspect",
        help="Inspect a FASTA input file and report module compatibility",
        description="Streaming FASTA validation: sequence type, counts, duplicates, warnings.",
    )
    _common(inspect)
    inspect.add_argument("path", help="Path to the FASTA file")
    inspect.add_argument(
        "--expected-type",
        default="auto",
        choices=("auto", "protein_fasta", "nucleotide_fasta", "protein", "nucleotide"),
        help="Expected sequence type (default: auto).",
    )
    inspect.add_argument(
        "--max-file-size-mb",
        type=int,
        default=500,
        help="Reject files larger than this size (default: 500).",
    )

    run = subparsers.add_parser(
        "run",
        help="Submit an analysis job and optionally wait for it",
        description="Run one or more analysis modules and write a unified AI report.",
    )
    _common(run)
    run.add_argument("input", help="Protein FASTA file")
    run.add_argument(
        "-m",
        "--modules",
        default="",
        help="Comma separated module list (default: all registered modules).",
    )
    run.add_argument(
        "-o",
        "--output",
        default=None,
        help="Requested output directory (recorded in job metadata).",
    )
    run.add_argument(
        "--input-type",
        default="auto",
        choices=("auto", "protein_fasta", "nucleotide_fasta", "protein", "nucleotide"),
        help="Declared input type (default: auto).",
    )
    run.add_argument(
        "--options-json",
        default=None,
        help="Per-module options as a JSON string, e.g. '{\"pathway\": {\"threads\": 8}}'.",
    )
    run.add_argument(
        "--timeout",
        type=float,
        default=None,
        help="Soft job timeout in seconds.",
    )
    run.add_argument(
        "--no-wait",
        action="store_true",
        help="Return immediately after submitting the job.",
    )
    run.add_argument(
        "--idempotency-key",
        default=None,
        help="Reuse an existing job with the same key instead of starting a new one.",
    )

    result = subparsers.add_parser(
        "result",
        help="Read the unified AI report of a job",
        description="Print a compact report summary, or one module result with --module.",
    )
    _common(result)
    result.add_argument("job_id", help="Job identifier")
    result.add_argument("--module", default=None, help="Return only this module result")

    artifact = subparsers.add_parser(
        "artifact",
        help="Read one page of a job artifact",
        description="List artifacts when no artifact id is given; otherwise read a bounded page.",
    )
    _common(artifact)
    artifact.add_argument("job_id", help="Job identifier")
    artifact.add_argument("artifact_id", nargs="?", default=None, help="Artifact identifier (a001, ...)")
    artifact.add_argument("--offset", type=int, default=0, help="First row/line offset (default: 0).")
    artifact.add_argument("--limit", type=int, default=50, help="Maximum rows/lines (default: 50).")
    artifact.add_argument("--filter-column", default=None, help="Column to filter on (CSV/TSV only).")
    artifact.add_argument("--filter-value", default=None, help="Value the filter column must equal.")
    artifact.add_argument(
        "--contains",
        action="store_true",
        help="Use substring matching for --filter-value.",
    )

    cancel = subparsers.add_parser(
        "cancel",
        help="Request cancellation of a running job",
        description="Flag the job for cancellation and stop the worker process.",
    )
    _common(cancel)
    cancel.add_argument("job_id", help="Job identifier")

    jobs = subparsers.add_parser(
        "jobs",
        help="List jobs in the workspace",
        description="List job identifiers and their current status.",
    )
    _common(jobs)
    jobs.add_argument("--limit", type=int, default=20, help="Maximum number of jobs (default: 20).")


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def _request_id() -> str:
    return "req_" + uuid.uuid4().hex[:12]


def _write_json_stdout(payload: Dict[str, Any]) -> None:
    """Write exactly one JSON object to stdout as UTF-8 bytes.

    Writing to ``sys.stdout.buffer`` keeps the stream valid JSON even when the
    console code page is not UTF-8 (default on Chinese Windows).
    """
    text = json.dumps(payload, ensure_ascii=False, indent=2) + "\n"
    buffer = getattr(sys.stdout, "buffer", None)
    if buffer is not None:
        buffer.write(text.encode("utf-8"))
        buffer.flush()
    else:  # pragma: no cover - exotic stream replacement
        sys.stdout.write(text)
        sys.stdout.flush()


def _emit(payload: Dict[str, Any], args: argparse.Namespace, exit_code: int = EXIT_OK) -> int:
    """Write the payload to stdout (and ``--json-out``) and return the exit code."""
    if getattr(args, "json_out", None):
        try:
            atomic_write_json(args.json_out, payload)
        except OSError as exc:
            print("error: cannot write --json-out: {0}".format(exc), file=sys.stderr)
            return EXIT_ERROR

    if getattr(args, "format", "text") == "json":
        _write_json_stdout(payload)
        return exit_code

    _print_text(payload, args)
    return exit_code


def _print_text(payload: Dict[str, Any], args: argparse.Namespace) -> None:
    if not payload.get("ok"):
        error = payload.get("error") or {}
        print("error [{0}]: {1}".format(error.get("code"), error.get("message")), file=sys.stderr)
        return

    data = payload.get("data") or {}
    command = getattr(args, "command", None)

    if command == "capabilities":
        software = data.get("software", {})
        print("MethArCT {0} (AI API {1})".format(software.get("version"), software.get("ai_api_version")))
        print("Modules:")
        for module in data.get("modules", []):
            state = "available" if module.get("enabled") else module.get("status")
            print("  - {0:<16} {1}".format(module.get("name"), state))
            for reason in module.get("reasons", []):
                print("      {0}".format(reason))
        return

    if command == "doctor":
        print("Environment check: {0}".format("OK" if data.get("ok") else "ISSUES FOUND"))
        for check in data.get("checks", []):
            print("  [{0}] {1} - {2}".format(check.get("status"), check.get("id"), check.get("detail")))
        return

    if command == "inspect":
        print("File: {0}".format(data.get("filename")))
        print("Sequences: {0}   Residues: {1}".format(data.get("sequence_count"), data.get("total_residues")))
        print("Type: {0} (confidence: {1})".format(data.get("sequence_type"), data.get("sequence_type_confidence")))
        print("Valid: {0}".format(data.get("valid")))
        for warning in data.get("warnings", []):
            print("  warning: {0}".format(warning))
        for error in data.get("errors", []):
            print("  error: {0}".format(error))
        return

    if command == "run":
        print("Job: {0}".format(data.get("job_id")))
        print("Status: {0}".format(data.get("status")))
        print("Modules: {0}".format(", ".join(data.get("requested_modules", []))))
        for item in data.get("unavailable_modules", []):
            print("  unavailable: {0}".format(item.get("module")))
        return

    if command == "result":
        print("Job: {0}   status: {1}".format(data.get("job_id"), data.get("status")))
        for name, module in (data.get("modules") or {}).items():
            print("  - {0}: {1}".format(name, module.get("status")))
            if module.get("summary"):
                print("      {0}".format(module.get("summary")))
        return

    if command == "artifact":
        if "columns" in data:
            print("\t".join(data.get("columns") or []))
            for row in data.get("rows", []):
                print("\t".join(str(row.get(column)) for column in data.get("columns", [])))
            print(
                "# rows {0}/{1}, next offset: {2}".format(
                    data.get("returned_rows"), data.get("matched_rows"), data.get("next_offset")
                ),
                file=sys.stderr,
            )
        elif "content" in data:
            print(data.get("content"))
        else:
            for entry in data.get("artifacts", []):
                print(
                    "{0}\t{1}\t{2}\t{3} bytes\t{4}".format(
                        entry.get("artifact_id"),
                        entry.get("path"),
                        entry.get("media_type"),
                        entry.get("size_bytes"),
                        entry.get("module"),
                    )
                )
        return

    if command == "cancel":
        print("Job {0}: {1} (accepted: {2})".format(data.get("job_id"), data.get("status"), data.get("accepted")))
        return

    if command == "jobs":
        for job in data.get("jobs", []):
            print("{0}\t{1}\t{2}".format(job.get("job_id"), job.get("status"), job.get("input_filename")))
        return

    print(json.dumps(data, ensure_ascii=False, indent=2))


def _fail(args: argparse.Namespace, exc: MethArCTError) -> int:
    payload = {
        "ok": False,
        "request_id": _request_id(),
        "data": None,
        "warnings": [],
        "error": exc.to_dict(),
    }
    return _emit(payload, args, EXIT_ERROR)


def _client(args: argparse.Namespace):
    from ..ai.api import MethArCTClient

    return MethArCTClient(
        workspace=getattr(args, "workspace", DEFAULT_WORKSPACE),
        config_file=getattr(args, "config", None),
    )


# ---------------------------------------------------------------------------
# Command handlers
# ---------------------------------------------------------------------------

def _handle_capabilities(args: argparse.Namespace) -> int:
    client = _client(args)
    data = client.get_capabilities(include_details=not args.no_details)
    return _emit(
        {"ok": True, "request_id": _request_id(), "data": data, "warnings": [], "error": None},
        args,
    )


def _handle_doctor(args: argparse.Namespace) -> int:
    client = _client(args)
    data = client.doctor(check_databases=not args.skip_databases)
    warnings = [issue["detail"] for issue in data.get("blocking_issues", [])]
    return _emit(
        {
            "ok": True,
            "request_id": _request_id(),
            "data": data,
            "warnings": warnings,
            "error": None,
        },
        args,
    )


def _handle_inspect(args: argparse.Namespace) -> int:
    client = _client(args)
    data = client.inspect_input(
        args.path,
        expected_type=args.expected_type,
        max_file_size_mb=args.max_file_size_mb,
    )
    return _emit(
        {
            "ok": True,
            "request_id": _request_id(),
            "data": data,
            "warnings": data.get("warnings", []),
            "error": None,
        },
        args,
        EXIT_OK if data.get("valid") else EXIT_JOB_FAILED,
    )


def _parse_options(raw: Optional[str]) -> Dict[str, Any]:
    if not raw:
        return {}
    try:
        payload = json.loads(raw)
    except ValueError as exc:
        raise MethArCTError(
            ErrorCode.INVALID_ARGUMENT,
            "--options-json is not valid JSON.",
            details={"error": str(exc)},
        )
    if not isinstance(payload, dict):
        raise MethArCTError(
            ErrorCode.INVALID_ARGUMENT,
            "--options-json must decode to an object.",
            details={"type": type(payload).__name__},
        )
    return payload


def _handle_run(args: argparse.Namespace) -> int:
    client = _client(args)
    modules: List[str] = [item.strip() for item in (args.modules or "").split(",") if item.strip()]
    options = _parse_options(args.options_json)
    job = client.run_analysis(
        input_file=args.input,
        input_type=args.input_type,
        modules=modules or None,
        output_dir=args.output,
        options=options,
        timeout_seconds=args.timeout,
        idempotency_key=args.idempotency_key,
        wait=not args.no_wait,
    )
    status = job.get("status")
    exit_code = EXIT_OK
    if status in ("failed", "cancelled"):
        exit_code = EXIT_JOB_FAILED
    elif status == "partially_succeeded":
        exit_code = EXIT_JOB_FAILED
    elif args.no_wait:
        exit_code = EXIT_OK
    return _emit(
        {
            "ok": True,
            "request_id": _request_id(),
            "data": job,
            "warnings": job.get("warnings", []),
            "error": None,
        },
        args,
        exit_code,
    )


def _handle_result(args: argparse.Namespace) -> int:
    client = _client(args)
    if args.module:
        data = client.get_module_result(args.job_id, args.module)
        data = {"job_id": args.job_id, "modules": {args.module: data}, "status": client.get_job_status(args.job_id).get("status")}
    else:
        data = client.get_result_summary(args.job_id)
    return _emit(
        {
            "ok": True,
            "request_id": _request_id(),
            "data": data,
            "warnings": data.get("warnings", []),
            "error": None,
        },
        args,
    )


def _handle_artifact(args: argparse.Namespace) -> int:
    client = _client(args)
    if not args.artifact_id:
        data = {"job_id": args.job_id, "artifacts": client.list_artifacts(args.job_id)}
    else:
        data = client.read_artifact(
            args.job_id,
            args.artifact_id,
            offset=args.offset,
            limit=args.limit,
            filter_column=args.filter_column,
            filter_value=args.filter_value,
            contains=args.contains,
        )
    return _emit(
        {"ok": True, "request_id": _request_id(), "data": data, "warnings": [], "error": None},
        args,
    )


def _handle_cancel(args: argparse.Namespace) -> int:
    client = _client(args)
    data = client.cancel_job(args.job_id)
    return _emit(
        {"ok": True, "request_id": _request_id(), "data": data, "warnings": [], "error": None},
        args,
    )


def _handle_jobs(args: argparse.Namespace) -> int:
    client = _client(args)
    data = {"workspace": str(client.workspace), "jobs": client.list_jobs(limit=args.limit)}
    return _emit(
        {"ok": True, "request_id": _request_id(), "data": data, "warnings": [], "error": None},
        args,
    )


_HANDLERS = {
    "capabilities": _handle_capabilities,
    "doctor": _handle_doctor,
    "inspect": _handle_inspect,
    "run": _handle_run,
    "result": _handle_result,
    "artifact": _handle_artifact,
    "cancel": _handle_cancel,
    "jobs": _handle_jobs,
}


def run_ai_command(args: argparse.Namespace) -> int:
    """Dispatch an AI subcommand and return the process exit code."""
    handler = _HANDLERS.get(args.command)
    if handler is None:  # pragma: no cover - argparse prevents this
        return EXIT_USAGE
    try:
        return handler(args)
    except MethArCTError as exc:
        return _fail(args, exc)
    except KeyboardInterrupt:  # pragma: no cover - interactive only
        print("interrupted", file=sys.stderr)
        return EXIT_ERROR
    except Exception as exc:  # noqa: BLE001 - mapped to a stable error code
        if getattr(args, "verbose", False):
            import traceback

            traceback.print_exc()
        return _fail(args, to_metharct_error(exc))
