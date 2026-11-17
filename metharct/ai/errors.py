# -*- coding: utf-8 -*-
"""Stable error codes and exceptions for the MethArCT AI layer.

Error codes are part of the public contract.  Once released they stay stable;
extra context belongs in ``details``.
"""

import subprocess
from typing import Any, Dict, Optional


class ErrorCode(object):
    """Container for the stable error code strings."""

    INPUT_FILE_NOT_FOUND = "INPUT_FILE_NOT_FOUND"
    INPUT_FILE_UNREADABLE = "INPUT_FILE_UNREADABLE"
    INPUT_FILE_TOO_LARGE = "INPUT_FILE_TOO_LARGE"
    INPUT_FORMAT_INVALID = "INPUT_FORMAT_INVALID"
    INPUT_TYPE_AMBIGUOUS = "INPUT_TYPE_AMBIGUOUS"
    INPUT_MODULE_INCOMPATIBLE = "INPUT_MODULE_INCOMPATIBLE"

    MODULE_UNKNOWN = "MODULE_UNKNOWN"
    MODULE_UNAVAILABLE = "MODULE_UNAVAILABLE"

    DEPENDENCY_MISSING = "DEPENDENCY_MISSING"
    DATABASE_MISSING = "DATABASE_MISSING"
    DATABASE_INVALID = "DATABASE_INVALID"

    ANALYSIS_TIMEOUT = "ANALYSIS_TIMEOUT"
    ANALYSIS_FAILED = "ANALYSIS_FAILED"

    JOB_NOT_FOUND = "JOB_NOT_FOUND"
    JOB_NOT_RUNNING = "JOB_NOT_RUNNING"
    JOB_ALREADY_FINISHED = "JOB_ALREADY_FINISHED"

    ARTIFACT_NOT_FOUND = "ARTIFACT_NOT_FOUND"
    ARTIFACT_ACCESS_DENIED = "ARTIFACT_ACCESS_DENIED"

    SCHEMA_VALIDATION_FAILED = "SCHEMA_VALIDATION_FAILED"
    INVALID_ARGUMENT = "INVALID_ARGUMENT"
    INTERNAL_ERROR = "INTERNAL_ERROR"


ALL_ERROR_CODES = tuple(
    value
    for name, value in sorted(vars(ErrorCode).items())
    if not name.startswith("_") and isinstance(value, str)
)


class MethArCTError(Exception):
    """Exception carrying a stable error code plus machine-readable details."""

    def __init__(
        self,
        code: str,
        message: str,
        retryable: bool = False,
        details: Optional[Dict[str, Any]] = None,
    ):
        super(MethArCTError, self).__init__(message)
        self.code = code
        self.message = message
        self.retryable = bool(retryable)
        self.details = details or {}

    def to_dict(self) -> Dict[str, Any]:
        """Return the error as the JSON envelope used by MCP / CLI JSON mode."""
        return {
            "code": self.code,
            "message": self.message,
            "retryable": self.retryable,
            "details": self.details,
        }

    def __str__(self) -> str:  # pragma: no cover - trivial
        return "{0}: {1}".format(self.code, self.message)


def error(
    code: str,
    message: str,
    retryable: bool = False,
    details: Optional[Dict[str, Any]] = None,
) -> MethArCTError:
    """Convenience factory for :class:`MethArCTError`."""
    return MethArCTError(code, message, retryable=retryable, details=details)


def to_metharct_error(exc: BaseException) -> MethArCTError:
    """Best-effort mapping of arbitrary exceptions onto stable error codes.

    The original exception type and message are preserved in ``details`` so no
    diagnostic information is lost while the agent-facing surface stays stable.
    """
    if isinstance(exc, MethArCTError):
        return exc

    details = {"exception_type": type(exc).__name__}
    text = str(exc)
    if text:
        details["exception_message"] = text

    if isinstance(exc, FileNotFoundError):
        return MethArCTError(ErrorCode.INPUT_FILE_NOT_FOUND, text or "File not found.", details=details)
    if isinstance(exc, PermissionError):
        return MethArCTError(ErrorCode.INPUT_FILE_UNREADABLE, text or "File not readable.", details=details)
    if isinstance(exc, subprocess.TimeoutExpired):
        return MethArCTError(
            ErrorCode.ANALYSIS_TIMEOUT,
            "Analysis exceeded the configured timeout.",
            retryable=True,
            details=details,
        )
    if isinstance(exc, (ValueError, TypeError)):
        return MethArCTError(ErrorCode.INVALID_ARGUMENT, text or "Invalid argument.", details=details)

    return MethArCTError(ErrorCode.INTERNAL_ERROR, text or "Unclassified internal error.", details=details)
