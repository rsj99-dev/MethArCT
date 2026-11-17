# -*- coding: utf-8 -*-
"""Background worker that executes one MethArCT job.

The worker is started as ``python -m metharct.ai.worker --workspace <ws>
--job-id <id>`` in a separate process so a long DIAMOND run cannot block the
API server, the CLI or an MCP host.  It updates the job record after every
module, so an interrupted worker still leaves a truthful state behind.
"""

import argparse
import os
import sys
import time
import traceback
from typing import Any, Dict, List, Optional

from .artifacts import ArtifactRegistry
from .constants import (
    AI_REPORT_FILENAME,
    JOB_STATUS_CANCELLED,
    JOB_STATUS_FAILED,
    JOB_STATUS_PARTIALLY_SUCCEEDED,
    JOB_STATUS_RUNNING,
    JOB_STATUS_SUCCEEDED,
    JOB_STATUS_VALIDATING,
    MODULE_STATUS_FAILED,
    MODULE_STATUS_SUCCESS,
)
from .errors import ErrorCode, MethArCTError, to_metharct_error
from .jobs import JobManager
from .report import assert_report_valid, build_report, write_report
from .runner import build_output_prefix, execute_module, ordered_modules
from .serializers import serialize_module
from .validators import inspect_fasta


def _log(message: str) -> None:
    """Write a progress line to stdout (captured in ``logs/stdout.log``)."""
    print("[metharct-worker] {0}".format(message), flush=True)


def _final_status(module_results: Dict[str, Any]) -> str:
    """Derive the job status from the individual module outcomes."""
    statuses = [result.get("status") for result in module_results.values()]
    successes = [status for status in statuses if status == MODULE_STATUS_SUCCESS]
    failures = [status for status in statuses if status == MODULE_STATUS_FAILED]
    if failures and successes:
        return JOB_STATUS_PARTIALLY_SUCCEEDED
    if failures:
        return JOB_STATUS_FAILED
    if not successes:
        return JOB_STATUS_FAILED
    return JOB_STATUS_SUCCEEDED


def run_job(workspace: str, job_id: str, config_file: Optional[str] = None) -> Dict[str, Any]:
    """Execute *job_id* and write the unified report.

    Returns the report dictionary that was written to disk.
    """
    from ..utils.config import Config

    manager = JobManager(workspace, config_file=config_file)
    job = manager.load(job_id)

    if job.get("status") in (JOB_STATUS_SUCCEEDED, JOB_STATUS_PARTIALLY_SUCCEEDED, JOB_STATUS_CANCELLED):
        _log("job {0} is already finished ({1})".format(job_id, job.get("status")))
        return {}

    directory = str(manager.job_dir(job_id))
    artifact_dir = os.path.join(directory, "artifacts")
    os.makedirs(artifact_dir, exist_ok=True)

    job = manager.transition(job_id, JOB_STATUS_VALIDATING, phase="validating")
    manager.update(job_id, pid=os.getpid())

    try:
        inspection = inspect_fasta(
            job["input_file"],
            expected_type="auto",
            max_file_size_mb=int(job.get("input_evidence", {}).get("max_file_size_mb") or 500),
        )
    except MethArCTError as exc:
        manager.append_error(job_id, exc.code, exc.message)
        manager.transition(job_id, JOB_STATUS_FAILED, phase="validation_failed")
        return {}

    recorded_sha = job.get("input_sha256")
    if recorded_sha and inspection.get("sha256") and recorded_sha != inspection.get("sha256"):
        message = "Input file changed between inspection and analysis; refusing to run."
        manager.append_error(job_id, ErrorCode.INPUT_FORMAT_INVALID, message)
        manager.transition(job_id, JOB_STATUS_FAILED, phase="validation_failed")
        return {}

    job = manager.transition(job_id, JOB_STATUS_RUNNING, phase="running")
    _log("job {0} started (pid {1})".format(job_id, os.getpid()))

    config = Config(config_file or job.get("config_file"))
    config.set("output.base_dir", artifact_dir)

    requested = ordered_modules(job.get("requested_modules"))
    options_by_module = job.get("options") or {}
    deadline = None
    if job.get("timeout_seconds"):
        deadline = time.time() + float(job["timeout_seconds"])

    native_results: Dict[str, Any] = {}
    completed: List[str] = []
    failed: List[str] = []
    skipped: List[str] = []
    cancel_hit = False

    for index, module_name in enumerate(requested, start=1):
        current = manager.load(job_id)
        if current.get("cancel_requested"):
            cancel_hit = True
            manager.append_warning(job_id, "Cancellation was requested; remaining modules were not started.")
            break
        if deadline is not None and time.time() > deadline:
            manager.append_error(
                job_id,
                ErrorCode.ANALYSIS_TIMEOUT,
                "Job timeout reached before module '{0}' started.".format(module_name),
            )
            failed.append(module_name)
            break

        manager.update(
            job_id,
            progress={"current_module": module_name, "completed": len(completed), "total": len(requested)},
        )
        _log("running module {0} ({1}/{2})".format(module_name, index, len(requested)))

        native = execute_module(
            module_name=module_name,
            config=config,
            input_file=job["input_file"],
            output_prefix=build_output_prefix(artifact_dir, job["input_file"]),
            options=options_by_module.get(module_name),
        )
        native_results[module_name] = native

        status = native.get("status")
        if status == "failed":
            failed.append(module_name)
            manager.append_error(
                job_id,
                native.get("error_code") or ErrorCode.ANALYSIS_FAILED,
                "{0}: {1}".format(module_name, native.get("error") or "analysis failed"),
            )
        elif status in ("skipped", "not_available"):
            skipped.append(module_name)
            manager.append_warning(
                job_id,
                "{0} was skipped: {1}".format(module_name, native.get("message") or "unavailable"),
            )
        else:
            completed.append(module_name)

        manager.update(
            job_id,
            completed_modules=completed,
            failed_modules=failed,
            skipped_modules=skipped,
            progress={"current_module": None, "completed": len(completed), "total": len(requested)},
        )
        _log("module {0} finished with status {1}".format(module_name, status))

    registry = ArtifactRegistry(
        directory,
        artifact_dir,
        redactions=[
            (job["input_file"], os.path.basename(job["input_file"])),
            (directory, "."),
            (str(manager.workspace), "<workspace>"),
        ],
    )
    registry.scan()
    registry.write_manifest()

    def _resolve_artifact(path: str) -> Optional[str]:
        entry = registry.find_by_path(path)
        return entry["artifact_id"] if entry else None

    module_results = {}
    for module_name in requested:
        native = native_results.get(module_name)
        module_results[module_name] = serialize_module(
            module_name,
            native,
            runtime_seconds=(native or {}).get("runtime_seconds"),
            artifact_resolver=_resolve_artifact,
        )

    for entry in registry.entries:
        module = entry.get("module")
        if not module or module not in module_results:
            continue
        result = module_results[module]
        sources = result.setdefault("source_artifacts", [])
        known = {item.get("artifact_id") for item in sources if isinstance(item, dict)}
        if entry["artifact_id"] not in known:
            sources.append({"artifact_id": entry["artifact_id"], "kind": "output"})

    status = JOB_STATUS_CANCELLED if cancel_hit else _final_status(module_results)
    job = manager.update(
        job_id,
        artifacts=registry.entries,
        completed_modules=completed,
        failed_modules=failed,
        skipped_modules=skipped,
        input_sha256=inspection.get("sha256") or job.get("input_sha256"),
        input_evidence={
            "sequence_count": inspection.get("sequence_count"),
            "total_residues": inspection.get("total_residues"),
            "sequence_type": inspection.get("sequence_type"),
            "sequence_type_confidence": inspection.get("sequence_type_confidence"),
            "file_size_bytes": inspection.get("file_size_bytes"),
            "max_file_size_mb": 500,
        },
    )
    job = manager.transition(job_id, status, phase=status)

    report = build_report(
        job=job,
        input_info=inspection,
        module_results=module_results,
        artifacts=registry.entries,
        config=config,
    )
    assert_report_valid(report)
    write_report(directory, report)

    manager.update(job_id, report_path=AI_REPORT_FILENAME)
    _log("job {0} finished with status {1}".format(job_id, job.get("status")))
    return report


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="metharct.ai.worker",
        description="Execute one MethArCT AI job (internal entry point).",
    )
    parser.add_argument("--workspace", required=True, help="Workspace directory containing jobs/")
    parser.add_argument("--job-id", required=True, help="Job identifier to execute")
    parser.add_argument("--config", default=None, help="Optional MethArCT configuration file")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entry point used by the process launcher."""
    args = _build_parser().parse_args(argv)
    try:
        run_job(args.workspace, args.job_id, config_file=args.config)
    except MethArCTError as exc:
        print("[metharct-worker] error: {0}".format(exc.to_dict()), file=sys.stderr, flush=True)
        return 1
    except Exception as exc:  # noqa: BLE001 - last-resort reporting
        traceback.print_exc()
        mapped = to_metharct_error(exc)
        print("[metharct-worker] unhandled error: {0}".format(mapped.to_dict()), file=sys.stderr, flush=True)
        try:
            manager = JobManager(args.workspace)
            manager.append_error(args.job_id, mapped.code, mapped.message)
            manager.transition(args.job_id, JOB_STATUS_FAILED, phase="crashed")
        except Exception:  # pragma: no cover - the job record may be gone
            pass
        return 1
    return 0


if __name__ == "__main__":  # pragma: no cover - process entry point
    sys.exit(main())
