# -*- coding: utf-8 -*-
"""Public MethArCT AI API: :class:`MethArCTClient`.

This is the only entry point agents should use.  It wraps the existing
analyzer classes; it never re-implements scientific logic.
"""

import os
import signal
import subprocess
import sys
import time
from typing import Any, Dict, Iterable, List, Optional

from .artifacts import ArtifactRegistry
from .capabilities import (
    MODULE_NAMES,
    compatible_modules,
    get_capabilities,
    get_module,
    normalize_module_list,
)
from .constants import (
    AI_API_VERSION,
    DEFAULT_ARTIFACT_LIMIT,
    DEFAULT_MAX_FILE_SIZE_MB,
    DEFAULT_WORKSPACE,
    ENV_ALLOWED_ROOTS,
    ENV_WORKSPACE,
    JOB_STATUS_CANCELLED,
    JOB_STATUS_FAILED,
    SCHEMA_VERSION,
    TERMINAL_JOB_STATUSES,
)
from .doctor import run_doctor
from .errors import ErrorCode, MethArCTError
from .jobs import JobManager, summarize_status
from .paths import enforce_allowed_roots, is_within, normalize
from .report import validate_payload, REPORT_SCHEMA_FILENAME
from .validators import inspect_fasta

#: Accepted ``input_type`` values (``auto`` is resolved during inspection).
_INPUT_TYPES = {
    "auto": "auto",
    "protein": "protein_fasta",
    "protein_fasta": "protein_fasta",
    "nucleotide": "nucleotide_fasta",
    "nucleotide_fasta": "nucleotide_fasta",
}


def _package_root() -> str:
    """Return the repository root that contains the ``metharct`` package."""
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _env_allowed_roots() -> Optional[List[str]]:
    raw = os.environ.get(ENV_ALLOWED_ROOTS)
    if not raw:
        return None
    return [part for part in raw.split(os.pathsep) if part.strip()]


def _process_alive(pid: Optional[int]) -> bool:
    """Return True when process *pid* is still running.

    Uses ``os.kill(pid, 0)`` on POSIX and ``OpenProcess`` on Windows; no
    external dependency (``psutil``) is required.
    """
    if not pid:
        return False
    pid = int(pid)
    if os.name == "nt":
        try:
            import ctypes

            query_limited_information = 0x1000
            handle = ctypes.windll.kernel32.OpenProcess(  # type: ignore[attr-defined]
                query_limited_information, False, pid
            )
            if not handle:
                return False
            ctypes.windll.kernel32.CloseHandle(handle)  # type: ignore[attr-defined]
            return True
        except Exception:  # pragma: no cover - unusual Windows environments
            return True
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:  # pragma: no cover - process owned by another user
        return True
    except OSError:  # pragma: no cover - defensive
        return True
    return True


class MethArCTClient(object):
    """Stable, synchronised facade over the MethArCT analysis pipeline."""

    def __init__(
        self,
        workspace: Optional[str] = None,
        allowed_roots: Optional[Iterable[str]] = None,
        config_file: Optional[str] = None,
        max_parallel_jobs: int = 1,
        max_file_size_mb: int = DEFAULT_MAX_FILE_SIZE_MB,
        default_timeout_seconds: Optional[float] = None,
    ):
        workspace = workspace or os.environ.get(ENV_WORKSPACE) or DEFAULT_WORKSPACE
        self.workspace = normalize(workspace)
        self.workspace.mkdir(parents=True, exist_ok=True)
        self.config_file = config_file
        self.max_parallel_jobs = max(int(max_parallel_jobs or 1), 1)
        self.max_file_size_mb = int(max_file_size_mb or DEFAULT_MAX_FILE_SIZE_MB)
        self.default_timeout_seconds = default_timeout_seconds
        self.allowed_roots = list(allowed_roots) if allowed_roots is not None else _env_allowed_roots()
        self._manager = JobManager(str(self.workspace), config_file=config_file)

    # -- metadata ----------------------------------------------------------
    @property
    def manager(self) -> JobManager:
        return self._manager

    def get_version_info(self) -> Dict[str, Any]:
        from .. import __version__

        return {
            "software": "MethArCT",
            "version": __version__,
            "ai_api_version": AI_API_VERSION,
            "schema_version": SCHEMA_VERSION,
            "workspace": str(self.workspace),
            "modules": list(MODULE_NAMES),
        }

    def _load_config(self):
        from ..utils.config import Config

        return Config(self.config_file)

    # -- discovery ---------------------------------------------------------
    def get_capabilities(self, include_details: bool = True) -> Dict[str, Any]:
        """Return the dynamically probed capability description."""
        return get_capabilities(self._load_config(), include_details=include_details)

    def doctor(self, check_databases: bool = True) -> Dict[str, Any]:
        """Return an environment self-check report (never mutates the system)."""
        return run_doctor(
            self._load_config(),
            workspace=str(self.workspace),
            check_databases=check_databases,
        )

    # -- inputs ------------------------------------------------------------
    def inspect_input(
        self,
        path: str,
        expected_type: str = "auto",
        max_file_size_mb: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Inspect a FASTA file and annotate module compatibility."""
        resolved_type = self._resolve_input_type(expected_type)
        info = inspect_fasta(
            enforce_allowed_roots(path, self.allowed_roots),
            expected_type=resolved_type,
            max_file_size_mb=int(max_file_size_mb or self.max_file_size_mb),
            allowed_roots=self.allowed_roots,
        )
        compatible, incompatible = compatible_modules(info.get("sequence_type"))
        info["compatible_modules"] = compatible
        info["incompatible_modules"] = incompatible
        info["max_file_size_mb"] = int(max_file_size_mb or self.max_file_size_mb)
        return info

    @staticmethod
    def _resolve_input_type(input_type: Optional[str]) -> str:
        key = (input_type or "auto").strip().lower()
        if key not in _INPUT_TYPES:
            raise MethArCTError(
                ErrorCode.INVALID_ARGUMENT,
                "Unsupported input_type value.",
                details={"input_type": input_type, "allowed": sorted(_INPUT_TYPES)},
            )
        return _INPUT_TYPES[key]

    def _validate_module_compatibility(self, modules: List[str], sequence_type: str) -> None:
        compatible, incompatible = compatible_modules(sequence_type)
        blocked = [name for name in modules if name in incompatible and name not in compatible]
        if blocked:
            raise MethArCTError(
                ErrorCode.INPUT_MODULE_INCOMPATIBLE,
                "Requested modules cannot analyse this input sequence type.",
                details={
                    "sequence_type": sequence_type,
                    "modules": blocked,
                    "compatible_modules": compatible,
                },
            )

    def _availability(self) -> Dict[str, Dict[str, Any]]:
        capabilities = get_capabilities(self._load_config(), include_details=True)
        return {entry["name"]: entry for entry in capabilities["modules"]}

    # -- analysis ----------------------------------------------------------
    def run_analysis(
        self,
        input_file: str,
        input_type: str = "auto",
        modules: Optional[Iterable[str]] = None,
        output_dir: Optional[str] = None,
        options: Optional[Dict[str, Dict[str, Any]]] = None,
        timeout_seconds: Optional[float] = None,
        metadata: Optional[Dict[str, Any]] = None,
        idempotency_key: Optional[str] = None,
        wait: bool = True,
        poll_interval: float = 1.0,
    ) -> Dict[str, Any]:
        """Submit (and by default wait for) one analysis job.

        Args:
            input_file: Protein FASTA path.
            input_type: ``auto`` / ``protein_fasta`` / ``nucleotide_fasta``.
            modules: Module identifiers; ``None`` means every registered module.
            output_dir: Optional caller-preferred output directory inside the
                workspace.  MethArCT always writes into the job directory; the
                value is recorded in the job metadata for traceability.
            options: Per-module option dictionaries (``{module: {option: value}}``).
            timeout_seconds: Soft job timeout enforced between modules.
            metadata: Free-form caller notes (never used in calculations).
            idempotency_key: Reuse an existing job with the same key.
            wait: Block until the job reaches a terminal state.
            poll_interval: Polling interval in seconds while waiting.

        Returns:
            A job view dictionary (no host absolute paths).
        """
        resolved_type = self._resolve_input_type(input_type)
        module_list = normalize_module_list(modules)
        resolved_options = self._normalize_options(module_list, options)
        resolved_timeout = (
            float(timeout_seconds) if timeout_seconds is not None else self.default_timeout_seconds
        )

        inspection = self.inspect_input(input_file, expected_type=resolved_type)
        if not inspection["valid"]:
            raise MethArCTError(
                ErrorCode.INPUT_FORMAT_INVALID,
                "Input file failed validation.",
                details={"errors": inspection["errors"], "warnings": inspection["warnings"]},
            )
        self._validate_module_compatibility(module_list, inspection["sequence_type"])

        availability = self._availability()
        unavailable = [
            {"module": name, "reasons": availability.get(name, {}).get("reasons", [])}
            for name in module_list
            if not availability.get(name, {}).get("enabled", False)
        ]
        if len(unavailable) == len(module_list):
            raise MethArCTError(
                ErrorCode.MODULE_UNAVAILABLE,
                "None of the requested modules is available in this environment.",
                details={"modules": unavailable},
                retryable=False,
            )

        self._assert_capacity()

        if idempotency_key:
            existing = self._manager.find_by_idempotency_key(idempotency_key)
            if existing is not None:
                view = self._job_view(existing)
                view["deduplicated"] = True
                if wait and existing.get("status") not in TERMINAL_JOB_STATUSES:
                    existing = self._wait_for_job(existing["job_id"], poll_interval)
                    view = self._job_view(existing)
                    view["deduplicated"] = True
                return view

        job_metadata = dict(metadata or {})
        if output_dir:
            job_metadata["requested_output_dir"] = str(
                enforce_allowed_roots(output_dir, self.allowed_roots, code=ErrorCode.INVALID_ARGUMENT)
            )

        job = self._manager.create_job(
            input_file=str(normalize(input_file)),
            input_type=resolved_type,
            modules=module_list,
            options=resolved_options,
            metadata=job_metadata,
            timeout_seconds=resolved_timeout,
            idempotency_key=idempotency_key,
            input_sha256=inspection.get("sha256"),
            input_evidence={
                "sequence_count": inspection.get("sequence_count"),
                "total_residues": inspection.get("total_residues"),
                "sequence_type": inspection.get("sequence_type"),
                "sequence_type_confidence": inspection.get("sequence_type_confidence"),
                "file_size_bytes": inspection.get("file_size_bytes"),
                "max_file_size_mb": inspection.get("max_file_size_mb"),
            },
        )
        self._launch_worker(job["job_id"])

        if wait:
            job = self._wait_for_job(job["job_id"], poll_interval, timeout=resolved_timeout)

        view = self._job_view(job)
        view["unavailable_modules"] = unavailable
        return view

    def submit_analysis(self, input_file: str, **kwargs: Any) -> Dict[str, Any]:
        """Non-blocking convenience wrapper around :meth:`run_analysis`."""
        kwargs["wait"] = False
        return self.run_analysis(input_file, **kwargs)

    @staticmethod
    def _normalize_options(
        modules: List[str], options: Optional[Dict[str, Dict[str, Any]]]
    ) -> Dict[str, Dict[str, Any]]:
        if not options:
            return {}
        if not isinstance(options, dict):
            raise MethArCTError(
                ErrorCode.INVALID_ARGUMENT,
                "options must be a mapping of module name to option dictionary.",
                details={"type": type(options).__name__},
            )
        unknown_modules = [key for key in options if key not in MODULE_NAMES]
        if unknown_modules:
            raise MethArCTError(
                ErrorCode.MODULE_UNKNOWN,
                "Options were supplied for unknown modules.",
                details={"modules": unknown_modules, "available_modules": list(MODULE_NAMES)},
            )
        for name in modules:
            if name in options and not isinstance(options[name], dict):
                raise MethArCTError(
                    ErrorCode.INVALID_ARGUMENT,
                    "Module options must be a dictionary.",
                    details={"module": name, "type": type(options[name]).__name__},
                )
        return {name: dict(payload) for name, payload in options.items() if name in modules}

    def _assert_capacity(self) -> None:
        active = self._manager.active_jobs()
        if len(active) >= self.max_parallel_jobs:
            raise MethArCTError(
                ErrorCode.INVALID_ARGUMENT,
                "The configured parallel job limit has been reached.",
                retryable=True,
                details={
                    "max_parallel_jobs": self.max_parallel_jobs,
                    "active_jobs": [job.get("job_id") for job in active],
                    "reason": "parallel_job_limit",
                },
            )

    def _launch_worker(self, job_id: str) -> subprocess.Popen:
        """Start the worker process with an argument array (never a shell)."""
        job_dir = str(self._manager.job_dir(job_id))
        log_dir = os.path.join(job_dir, "logs")
        os.makedirs(log_dir, exist_ok=True)
        stdout_path = os.path.join(log_dir, "stdout.log")
        stderr_path = os.path.join(log_dir, "stderr.log")

        command = [
            sys.executable,
            "-m",
            "metharct.ai.worker",
            "--workspace",
            str(self.workspace),
            "--job-id",
            job_id,
        ]
        if self.config_file:
            command.extend(["--config", str(self.config_file)])

        env = os.environ.copy()
        env["PYTHONIOENCODING"] = "utf-8"
        env["PYTHONUNBUFFERED"] = "1"
        package_root = _package_root()
        existing_path = env.get("PYTHONPATH")
        env["PYTHONPATH"] = (
            package_root + os.pathsep + existing_path if existing_path else package_root
        )

        stdout_handle = open(stdout_path, "wb")
        stderr_handle = open(stderr_path, "wb")
        try:
            process = subprocess.Popen(
                command,
                cwd=package_root,
                stdin=subprocess.DEVNULL,
                stdout=stdout_handle,
                stderr=stderr_handle,
                shell=False,
                env=env,
            )
        finally:
            stdout_handle.close()
            stderr_handle.close()

        self._manager.update(job_id, pid=process.pid)
        return process

    def _wait_for_job(
        self,
        job_id: str,
        poll_interval: float = 1.0,
        timeout: Optional[float] = None,
    ) -> Dict[str, Any]:
        deadline = None if not timeout else time.time() + float(timeout) + 5.0
        interval = max(0.05, float(poll_interval or 0.2))
        dead_checks = 0
        while True:
            job = self._manager.load(job_id)
            if job.get("status") in TERMINAL_JOB_STATUSES:
                return job
            if deadline is not None and time.time() > deadline:
                self._terminate_worker(job_id)
                job = self._manager.transition(
                    job_id,
                    JOB_STATUS_FAILED,
                    phase="timeout",
                    error={
                        "code": ErrorCode.ANALYSIS_TIMEOUT,
                        "message": "Job did not finish within the configured timeout.",
                    },
                )
                return job
            if _process_alive(job.get("pid")):
                dead_checks = 0
            else:
                # Give the worker a moment to persist its final state before
                # declaring the process lost.
                dead_checks += 1
                if dead_checks >= 3:
                    return self._manager.transition(
                        job_id,
                        JOB_STATUS_FAILED,
                        phase="worker_lost",
                        error={
                            "code": ErrorCode.INTERNAL_ERROR,
                            "message": "Worker process exited without writing a final job state.",
                        },
                    )
            time.sleep(interval)

    # -- job access --------------------------------------------------------
    def _job_view(self, job: Dict[str, Any]) -> Dict[str, Any]:
        view = JobManager.public_view(job)
        view["status_tool"] = "metharct_get_job_status"
        view["result_tool"] = "metharct_get_result_summary"
        view["workspace_job_dir"] = "jobs/{0}".format(job.get("job_id"))
        return view

    def get_job(self, job_id: str) -> Dict[str, Any]:
        return self._job_view(self._manager.load(job_id))

    def get_job_status(self, job_id: str) -> Dict[str, Any]:
        return summarize_status(self._manager.load(job_id))

    def list_jobs(self, limit: Optional[int] = None) -> List[Dict[str, Any]]:
        return [self._job_view(job) for job in self._manager.list_jobs(limit=limit)]

    def cancel_job(self, job_id: str) -> Dict[str, Any]:
        """Request cancellation and terminate the worker when possible."""
        job = self._manager.load(job_id)
        if job.get("status") in TERMINAL_JOB_STATUSES:
            return {
                "job_id": job_id,
                "cancellation_supported": True,
                "accepted": False,
                "status": job.get("status"),
                "message": "Job already finished.",
            }

        self._manager.request_cancel(job_id)
        terminated = self._terminate_worker(job_id)
        final = self._manager.load(job_id)
        if final.get("status") not in TERMINAL_JOB_STATUSES:
            final = self._manager.transition(
                job_id,
                JOB_STATUS_CANCELLED,
                phase="cancelled",
                error={"code": "CANCELLED", "message": "Cancellation requested by the caller."},
            )
        return {
            "job_id": job_id,
            "cancellation_supported": True,
            "accepted": True,
            "worker_terminated": terminated,
            "status": final.get("status"),
            "message": "Cancellation accepted; the worker was stopped."
            if terminated
            else "Cancellation flag set; the worker will stop before the next module.",
        }

    def _terminate_worker(self, job_id: str) -> bool:
        """Best-effort termination of the worker process.  Never uses a shell."""
        job = self._manager.load(job_id)
        pid = job.get("pid")
        if not pid:
            return False
        try:
            if os.name == "nt":
                subprocess.run(
                    ["taskkill", "/PID", str(pid), "/T", "/F"],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    shell=False,
                    timeout=20,
                )
            else:
                os.kill(int(pid), signal.SIGTERM)
        except (OSError, subprocess.SubprocessError):
            return False
        return True

    # -- results -----------------------------------------------------------
    def _report_path(self, job: Dict[str, Any]) -> str:
        job_dir = str(self._manager.job_dir(job["job_id"]))
        return os.path.join(job_dir, job.get("report_path") or "metharct_ai_report.json")

    def get_result(self, job_id: str, validate: bool = False) -> Dict[str, Any]:
        """Return the full unified report for *job_id*."""
        job = self._manager.load(job_id)
        path = self._report_path(job)
        from .paths import read_json

        report = read_json(path)
        if not isinstance(report, dict):
            raise MethArCTError(
                ErrorCode.ARTIFACT_NOT_FOUND,
                "The AI report has not been written yet.",
                details={"job_id": job_id, "job_status": job.get("status")},
            )
        if validate:
            problems = validate_payload(report, REPORT_SCHEMA_FILENAME)
            if problems:
                raise MethArCTError(
                    ErrorCode.SCHEMA_VALIDATION_FAILED,
                    "Stored report does not satisfy the report schema.",
                    details={"problems": problems[:10]},
                )
        return report

    def get_result_summary(self, job_id: str, modules: Optional[Iterable[str]] = None) -> Dict[str, Any]:
        """Return a compact summary suitable for an LLM context window."""
        job = self._manager.load(job_id)
        try:
            report = self.get_result(job_id)
        except MethArCTError:
            summary = summarize_status(job)
            summary["report_available"] = False
            summary["modules"] = {}
            summary["input"] = {
                "filename": job.get("input_filename"),
                "type": job.get("input_type"),
                "sha256": job.get("input_sha256"),
            }
            summary["artifacts"] = job.get("artifacts", [])
            return summary

        requested = normalize_module_list(list(modules)) if modules else list(report.get("modules", {}))
        selected = {}
        for name in report.get("modules", {}):
            if name in requested:
                selected[name] = report["modules"][name]

        return {
            "job_id": job_id,
            "status": report.get("run", {}).get("status"),
            "report_available": True,
            "input": report.get("input", {}),
            "requested_modules": report.get("requested_modules", []),
            "modules": selected,
            "warnings": report.get("warnings", []),
            "errors": report.get("errors", []),
            "artifacts": report.get("artifacts", []),
            "provenance": report.get("provenance", {}),
        }

    def get_module_result(self, job_id: str, module: str) -> Dict[str, Any]:
        """Return one module result from the report."""
        get_module(module)  # validates the identifier
        report = self.get_result(job_id)
        module_result = report.get("modules", {}).get(module)
        if module_result is None:
            raise MethArCTError(
                ErrorCode.ARTIFACT_NOT_FOUND,
                "The module has no result in this job.",
                details={"job_id": job_id, "module": module, "modules": list(report.get("modules", {}))},
            )
        return module_result

    # -- artifacts ---------------------------------------------------------
    def _registry(self, job_id: str) -> ArtifactRegistry:
        job = self._manager.load(job_id)
        job_dir = str(self._manager.job_dir(job_id))
        input_file = job.get("input_file") or ""
        redactions = [
            (input_file, os.path.basename(input_file)) if input_file else ("", ""),
            (job_dir, "."),
            (str(self.workspace), "<workspace>"),
        ]
        return ArtifactRegistry(job_dir, os.path.join(job_dir, "artifacts"), redactions=redactions)

    def list_artifacts(self, job_id: str) -> List[Dict[str, Any]]:
        registry = self._registry(job_id)
        return registry.load_manifest()

    def read_artifact(
        self,
        job_id: str,
        artifact_id: str,
        offset: int = 0,
        limit: int = DEFAULT_ARTIFACT_LIMIT,
        filter_column: Optional[str] = None,
        filter_value: Optional[str] = None,
        contains: bool = False,
    ) -> Dict[str, Any]:
        """Read a bounded page of one artifact (manifest whitelist only)."""
        registry = self._registry(job_id)
        payload = registry.read(
            artifact_id,
            offset=offset,
            limit=limit,
            filter_column=filter_column,
            filter_value=filter_value,
            contains=contains,
        )
        payload["job_id"] = job_id
        return payload


def job_directory(client: MethArCTClient, job_id: str) -> str:
    """Return the absolute job directory (internal helper for tests/CLI)."""
    path = str(client.manager.job_dir(job_id))
    if not is_within(path, client.workspace):  # pragma: no cover - defensive
        raise MethArCTError(ErrorCode.ARTIFACT_ACCESS_DENIED, "Job directory escapes the workspace.")
    return path
