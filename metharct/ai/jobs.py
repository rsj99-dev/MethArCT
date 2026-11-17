# -*- coding: utf-8 -*-
"""Job lifecycle management (create / inspect / cancel).

A job owns a directory under ``<workspace>/jobs/<job_id>/`` containing
``job.json``, the unified AI report, the captured console logs and the
``artifacts/`` directory with the original result files.
"""

import os
import uuid
from typing import Any, Dict, Iterable, List, Optional

from .clock import parse_iso, utc_now as utc_now_dt, utc_now_iso
from .constants import (
    AI_REPORT_FILENAME,
    JOB_STATUS_CANCELLED,
    JOB_STATUS_FAILED,
    JOB_STATUS_PARTIALLY_SUCCEEDED,
    JOB_STATUS_QUEUED,
    JOB_STATUS_RUNNING,
    JOB_STATUS_SUCCEEDED,
    JOB_STATUS_VALIDATING,
    SCHEMA_VERSION,
    TERMINAL_JOB_STATUSES,
)
from .errors import ErrorCode, MethArCTError
from .paths import atomic_write_json, ensure_dir, normalize, read_json

JOB_RECORD_FILENAME = "job.json"

#: Allowed forward transitions of the job state machine.
_ALLOWED_TRANSITIONS = {
    JOB_STATUS_QUEUED: {JOB_STATUS_VALIDATING, JOB_STATUS_RUNNING, JOB_STATUS_FAILED, JOB_STATUS_CANCELLED},
    JOB_STATUS_VALIDATING: {JOB_STATUS_RUNNING, JOB_STATUS_FAILED, JOB_STATUS_CANCELLED},
    JOB_STATUS_RUNNING: {
        JOB_STATUS_SUCCEEDED,
        JOB_STATUS_PARTIALLY_SUCCEEDED,
        JOB_STATUS_FAILED,
        JOB_STATUS_CANCELLED,
    },
    JOB_STATUS_SUCCEEDED: set(),
    JOB_STATUS_PARTIALLY_SUCCEEDED: set(),
    JOB_STATUS_FAILED: set(),
    JOB_STATUS_CANCELLED: set(),
}


def utc_now() -> str:
    """Return the current UTC time as an ISO-8601 string with ``Z`` suffix."""
    return utc_now_iso()


def _duration_seconds(started_at: Optional[str], finished_at: Optional[str]) -> Optional[float]:
    start = parse_iso(started_at)
    end = parse_iso(finished_at)
    if start is None or end is None:
        return None
    return round((end - start).total_seconds(), 3)


class JobManager(object):
    """Filesystem-backed job store.

    The store is intentionally simple (one JSON document per job, atomic
    replacement) so an interrupted process can never leave a half-written job
    record behind.
    """

    def __init__(self, workspace: str, config_file: Optional[str] = None):
        self.workspace = normalize(workspace)
        self.config_file = config_file
        self.jobs_dir = self.workspace / "jobs"

    # -- layout ------------------------------------------------------------
    def ensure_layout(self) -> None:
        ensure_dir(self.jobs_dir)

    def job_dir(self, job_id: str) -> "os.PathLike[str]":
        if not job_id or os.sep in job_id or "/" in job_id or job_id.startswith("."):
            raise MethArCTError(
                ErrorCode.JOB_NOT_FOUND,
                "Invalid job id.",
                details={"job_id": job_id},
            )
        return self.jobs_dir / job_id

    def job_record_path(self, job_id: str) -> "os.PathLike[str]":
        return os.path.join(str(self.job_dir(job_id)), JOB_RECORD_FILENAME)

    # -- creation ----------------------------------------------------------
    @staticmethod
    def new_job_id() -> str:
        stamp = utc_now_dt().strftime("%Y%m%dT%H%M%S")
        return "job_{0}_{1}".format(stamp, uuid.uuid4().hex[:8])

    def create_job(
        self,
        input_file: str,
        input_type: str,
        modules: List[str],
        options: Optional[Dict[str, Any]] = None,
        metadata: Optional[Dict[str, Any]] = None,
        timeout_seconds: Optional[float] = None,
        idempotency_key: Optional[str] = None,
        input_sha256: Optional[str] = None,
        input_evidence: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """Create a new job directory and persist its ``job.json``."""
        self.ensure_layout()
        job_id = self.new_job_id()
        directory = ensure_dir(self.job_dir(job_id))
        ensure_dir(os.path.join(str(directory), "artifacts"))
        ensure_dir(os.path.join(str(directory), "logs"))

        job = {
            "schema_version": SCHEMA_VERSION,
            "job_id": job_id,
            "status": JOB_STATUS_QUEUED,
            "phase": "queued",
            "created_at": utc_now(),
            "started_at": None,
            "finished_at": None,
            "duration_seconds": None,
            "input_file": str(normalize(input_file)),
            "input_filename": os.path.basename(str(input_file)),
            "input_type": input_type,
            "input_sha256": input_sha256,
            "input_evidence": input_evidence or {},
            "requested_modules": list(modules),
            "completed_modules": [],
            "failed_modules": [],
            "skipped_modules": [],
            "options": options or {},
            "metadata": metadata or {},
            "timeout_seconds": timeout_seconds,
            "idempotency_key": idempotency_key,
            "report_path": AI_REPORT_FILENAME,
            "artifacts": [],
            "warnings": [],
            "errors": [],
            "error_code": None,
            "pid": None,
            "cancel_requested": False,
            "workspace": str(self.workspace),
            "config_file": self.config_file,
            "progress": {"current_module": None, "completed": 0, "total": len(modules)},
        }
        self.save(job)
        return job

    # -- persistence -------------------------------------------------------
    def save(self, job: Dict[str, Any]) -> Dict[str, Any]:
        job["duration_seconds"] = _duration_seconds(job.get("started_at"), job.get("finished_at"))
        atomic_write_json(self.job_record_path(job["job_id"]), job)
        return job

    def load(self, job_id: str) -> Dict[str, Any]:
        path = self.job_record_path(job_id)
        payload = read_json(path)
        if not isinstance(payload, dict):
            raise MethArCTError(
                ErrorCode.JOB_NOT_FOUND,
                "Job does not exist in this workspace.",
                details={"job_id": job_id, "workspace": str(self.workspace)},
            )
        return payload

    def exists(self, job_id: str) -> bool:
        try:
            return os.path.isfile(self.job_record_path(job_id))
        except MethArCTError:
            return False

    def list_jobs(self, limit: Optional[int] = None) -> List[Dict[str, Any]]:
        """Return job records sorted by creation time (newest first)."""
        if not os.path.isdir(str(self.jobs_dir)):
            return []
        jobs = []
        for name in os.listdir(str(self.jobs_dir)):
            record = read_json(os.path.join(str(self.jobs_dir), name, JOB_RECORD_FILENAME))
            if isinstance(record, dict):
                jobs.append(record)
        jobs.sort(
            key=lambda item: (item.get("created_at") or "", item.get("job_id") or ""),
            reverse=True,
        )
        if limit is not None:
            jobs = jobs[: int(limit)]
        return jobs

    # -- state machine -----------------------------------------------------
    def transition(
        self,
        job_id: str,
        status: str,
        phase: Optional[str] = None,
        error: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """Atomically move a job to *status* following the allowed transitions."""
        job = self.load(job_id)
        current = job.get("status", JOB_STATUS_QUEUED)
        if status == current:
            if phase:
                job["phase"] = phase
            return self.save(job)
        allowed = _ALLOWED_TRANSITIONS.get(current, set())
        if status not in allowed:
            raise MethArCTError(
                ErrorCode.JOB_ALREADY_FINISHED,
                "Illegal job state transition.",
                details={"job_id": job_id, "from": current, "to": status},
            )
        job["status"] = status
        job["phase"] = phase or status
        if status == JOB_STATUS_RUNNING and not job.get("started_at"):
            job["started_at"] = utc_now()
        if status in TERMINAL_JOB_STATUSES:
            job["finished_at"] = utc_now()
        if error:
            job.setdefault("errors", []).append(error)
            job["error_code"] = error.get("code")
        return self.save(job)

    def update(self, job_id: str, **fields: Any) -> Dict[str, Any]:
        """Merge *fields* into the job record and persist it."""
        job = self.load(job_id)
        job.update(fields)
        return self.save(job)

    def append_error(self, job_id: str, code: str, message: str) -> Dict[str, Any]:
        job = self.load(job_id)
        job.setdefault("errors", []).append({"code": code, "message": message})
        job["error_code"] = code
        return self.save(job)

    def append_warning(self, job_id: str, message: str) -> Dict[str, Any]:
        job = self.load(job_id)
        job.setdefault("warnings", []).append(message)
        return self.save(job)

    def request_cancel(self, job_id: str) -> Dict[str, Any]:
        """Flag a job for cancellation (the caller terminates the worker)."""
        job = self.load(job_id)
        if job.get("status") in TERMINAL_JOB_STATUSES:
            raise MethArCTError(
                ErrorCode.JOB_ALREADY_FINISHED,
                "Job already finished; nothing to cancel.",
                details={"job_id": job_id, "status": job.get("status")},
            )
        job["cancel_requested"] = True
        return self.save(job)

    def active_jobs(self) -> List[Dict[str, Any]]:
        return [
            job
            for job in self.list_jobs()
            if job.get("status") not in TERMINAL_JOB_STATUSES
        ]

    def find_by_idempotency_key(self, key: str) -> Optional[Dict[str, Any]]:
        """Return an existing unfinished-or-successful job for *key*."""
        if not key:
            return None
        for job in self.list_jobs():
            if job.get("idempotency_key") != key:
                continue
            if job.get("status") in (JOB_STATUS_FAILED, JOB_STATUS_CANCELLED):
                continue
            return job
        return None

    # -- presentation ------------------------------------------------------
    @staticmethod
    def public_view(job: Dict[str, Any]) -> Dict[str, Any]:
        """Return a job record without host absolute paths.

        Agents must never receive arbitrary absolute paths from the user's
        machine; relative paths inside the job directory are enough.
        """
        view = dict(job)
        input_file = view.pop("input_file", None)
        view.pop("workspace", None)
        if input_file:
            view["input_filename"] = os.path.basename(str(input_file))
        report_path = view.get("report_path")
        if report_path:
            view["report_path"] = str(report_path).replace(os.sep, "/")
        return view


def summarize_status(job: Dict[str, Any]) -> Dict[str, Any]:
    """Compact status payload shared by the CLI and MCP ``get_job_status``."""
    return {
        "job_id": job.get("job_id"),
        "status": job.get("status"),
        "phase": job.get("phase"),
        "created_at": job.get("created_at"),
        "started_at": job.get("started_at"),
        "finished_at": job.get("finished_at"),
        "duration_seconds": job.get("duration_seconds"),
        "requested_modules": job.get("requested_modules", []),
        "completed_modules": job.get("completed_modules", []),
        "failed_modules": job.get("failed_modules", []),
        "skipped_modules": job.get("skipped_modules", []),
        "progress": job.get("progress", {}),
        "error_code": job.get("error_code"),
        "errors": job.get("errors", []),
        "warnings": job.get("warnings", []),
        "cancel_requested": job.get("cancel_requested", False),
    }


def merge_unique(existing: Iterable[str], extra: Iterable[str]) -> List[str]:
    """Merge two name lists, preserving order and removing duplicates."""
    merged: List[str] = []
    for name in list(existing) + list(extra):
        if name not in merged:
            merged.append(name)
    return merged
