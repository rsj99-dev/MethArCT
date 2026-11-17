# -*- coding: utf-8 -*-
"""Stable constants for the MethArCT AI-ready adaptation layer.

These values are part of the public contract of the AI layer.  They must not
change without a documented compatibility note (see ``docs/AI_API.md``).
"""

#: Version of the unified AI report structure (``metharct_ai_report.json``).
SCHEMA_VERSION = "1.0"

#: Version of the public Python API (``MethArCTClient``).
AI_API_VERSION = "1.0.0"

#: Default upper bound for input file inspection (megabytes).
DEFAULT_MAX_FILE_SIZE_MB = 500

#: Default page size when reading tabular artifacts.
DEFAULT_ARTIFACT_LIMIT = 50

#: Hard upper bound for rows returned by a single artifact read.
MAX_ARTIFACT_LIMIT = 1000

#: Hard upper bound (bytes) for a single artifact read response.
MAX_ARTIFACT_BYTES = 2 * 1024 * 1024

#: Default workspace directory name used by the CLI.
DEFAULT_WORKSPACE = "metharct_workspace"

#: Environment variable that restricts accessible roots for agents.
ENV_ALLOWED_ROOTS = "METHARCT_ALLOWED_ROOTS"

#: Environment variable pointing at a MethArCT workspace.
ENV_WORKSPACE = "METHARCT_WORKSPACE"

#: Module result status values (see AI report schema).
MODULE_STATUS_SUCCESS = "success"
MODULE_STATUS_PARTIAL = "partial"
MODULE_STATUS_FAILED = "failed"
MODULE_STATUS_SKIPPED = "skipped"
MODULE_STATUS_NOT_AVAILABLE = "not_available"
MODULE_STATUS_NOT_APPLICABLE = "not_applicable"

MODULE_STATUSES = (
    MODULE_STATUS_SUCCESS,
    MODULE_STATUS_PARTIAL,
    MODULE_STATUS_FAILED,
    MODULE_STATUS_SKIPPED,
    MODULE_STATUS_NOT_AVAILABLE,
    MODULE_STATUS_NOT_APPLICABLE,
)

#: Job status values (see job schema).
JOB_STATUS_QUEUED = "queued"
JOB_STATUS_VALIDATING = "validating"
JOB_STATUS_RUNNING = "running"
JOB_STATUS_SUCCEEDED = "succeeded"
JOB_STATUS_PARTIALLY_SUCCEEDED = "partially_succeeded"
JOB_STATUS_FAILED = "failed"
JOB_STATUS_CANCELLED = "cancelled"

JOB_STATUSES = (
    JOB_STATUS_QUEUED,
    JOB_STATUS_VALIDATING,
    JOB_STATUS_RUNNING,
    JOB_STATUS_SUCCEEDED,
    JOB_STATUS_PARTIALLY_SUCCEEDED,
    JOB_STATUS_FAILED,
    JOB_STATUS_CANCELLED,
)

#: Terminal job statuses (no further state transition is expected).
TERMINAL_JOB_STATUSES = frozenset(
    [
        JOB_STATUS_SUCCEEDED,
        JOB_STATUS_PARTIALLY_SUCCEEDED,
        JOB_STATUS_FAILED,
        JOB_STATUS_CANCELLED,
    ]
)

#: Dependency / environment status values used by doctor and capabilities.
STATUS_AVAILABLE = "available"
STATUS_UNAVAILABLE = "unavailable"
STATUS_MISCONFIGURED = "misconfigured"
STATUS_NOT_REQUIRED = "not_required"

#: Name of the unified AI report file written for every job.
AI_REPORT_FILENAME = "metharct_ai_report.json"
