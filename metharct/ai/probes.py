# -*- coding: utf-8 -*-
"""Low-level environment probes shared by capabilities and doctor.

Every probe returns a plain dictionary with a stable ``status`` value
(``available`` / ``unavailable`` / ``misconfigured`` / ``not_required``) so no
diagnostic detail is lost behind a boolean.
"""

import importlib
import importlib.util
import os
import shlex
import subprocess
import sys
from typing import Any, Dict, List, Optional

from .constants import (
    STATUS_AVAILABLE,
    STATUS_MISCONFIGURED,
    STATUS_NOT_REQUIRED,
    STATUS_UNAVAILABLE,
)


def build_command(value: Any) -> List[str]:
    """Turn a configured command into an argv list without using a shell.

    ``"diamond"`` -> ``["diamond"]``; ``"wsl diamond"`` -> ``["wsl", "diamond"]``.
    An existing file path is kept as a single argument even if it contains
    spaces.
    """
    if isinstance(value, (list, tuple)):
        return [str(item) for item in value]
    text = str(value).strip()
    if not text:
        return []
    if os.path.isfile(text):
        return [text]
    try:
        if os.name == "nt":
            parts = shlex.split(text, posix=False)
        else:
            parts = shlex.split(text)
    except ValueError:  # pragma: no cover - unbalanced quotes
        parts = text.split()
    cleaned = [part.strip("\"'") for part in parts if part.strip("\"'")]
    return cleaned or [text]


def probe_python_module(module_name: str) -> Dict[str, Any]:
    """Check whether a Python module can be imported (without importing it)."""
    try:
        spec = importlib.util.find_spec(module_name)
    except (ImportError, ValueError, AttributeError):
        spec = None
    if spec is None:
        return {
            "name": "python_module:{0}".format(module_name),
            "status": STATUS_UNAVAILABLE,
            "detail": "Python module '{0}' is not installed.".format(module_name),
            "hint": "pip install {0}".format(module_name.split(".")[0]),
        }

    version = None
    try:
        module = importlib.import_module(module_name)
        version = getattr(module, "__version__", None)
    except Exception as exc:  # pragma: no cover - broken installs
        return {
            "name": "python_module:{0}".format(module_name),
            "status": STATUS_MISCONFIGURED,
            "detail": "Python module '{0}' failed to import: {1}".format(module_name, exc),
        }
    return {
        "name": "python_module:{0}".format(module_name),
        "status": STATUS_AVAILABLE,
        "detail": "Python module '{0}' is importable.".format(module_name),
        "version": version,
    }


def probe_package_version(distribution: str) -> Dict[str, Any]:
    """Read an installed distribution version (PEP 566 metadata)."""
    try:
        from importlib import metadata
    except ImportError:  # pragma: no cover - Python < 3.8
        import importlib_metadata as metadata  # type: ignore

    try:
        version = metadata.version(distribution)
    except Exception:
        return {
            "name": "package:{0}".format(distribution),
            "status": STATUS_UNAVAILABLE,
            "detail": "Distribution '{0}' is not installed.".format(distribution),
        }
    return {
        "name": "package:{0}".format(distribution),
        "status": STATUS_AVAILABLE,
        "detail": "Distribution '{0}' is installed.".format(distribution),
        "version": version,
    }


def probe_executable(command: Any, version_args: Optional[List[str]] = None, timeout: int = 20) -> Dict[str, Any]:
    """Check that an executable can be started and capture its version string."""
    argv = build_command(command)
    if not argv:
        return {
            "name": "executable",
            "status": STATUS_UNAVAILABLE,
            "detail": "No executable configured.",
        }

    # Only the program name (never the full configured path) enters the check
    # name, so the capability manifest stays path-free.
    program = os.path.basename(argv[0])
    check_name = "executable:{0}".format(program)
    args = list(version_args or ["--version"])
    try:
        result = subprocess.run(
            argv + args,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            shell=False,
            universal_newlines=True,
        )
    except FileNotFoundError:
        return {
            "name": check_name,
            "status": STATUS_UNAVAILABLE,
            "detail": "Executable not found on PATH.",
            "hint": "Install the tool or set its absolute path in metharct_config.yaml.",
            "command": argv[0],
        }
    except subprocess.TimeoutExpired:
        return {
            "name": check_name,
            "status": STATUS_MISCONFIGURED,
            "detail": "Executable did not respond within {0}s.".format(timeout),
            "command": argv[0],
        }
    except OSError as exc:
        return {
            "name": check_name,
            "status": STATUS_UNAVAILABLE,
            "detail": "Executable could not be started: {0}".format(exc),
            "command": argv[0],
        }

    output = (result.stdout or "").strip() or (result.stderr or "").strip()
    first_line = output.splitlines()[0].strip() if output else ""
    if result.returncode == 0:
        return {
            "name": check_name,
            "status": STATUS_AVAILABLE,
            "detail": "Executable responded successfully.",
            "version": first_line or None,
            "command": argv[0],
        }
    return {
        "name": check_name,
        "status": STATUS_MISCONFIGURED,
        "detail": "Executable returned exit code {0}: {1}".format(result.returncode, first_line),
        "command": argv[0],
    }


def probe_path(path: Optional[str], kind: str = "any") -> Dict[str, Any]:
    """Check that a filesystem path exists and matches *kind*.

    The check ``name`` intentionally carries only the basename so the
    capability manifest stays free of host absolute paths; the full path is
    available under ``path`` for the doctor report.
    """
    label = os.path.basename(str(path).rstrip("/\\")) if path else "<unset>"
    if not path:
        return {
            "name": "path:{0}".format(label),
            "status": STATUS_UNAVAILABLE,
            "detail": "Path is not configured.",
        }
    if kind == "dir":
        exists = os.path.isdir(path)
    elif kind == "file":
        exists = os.path.isfile(path)
    else:
        exists = os.path.exists(path)
    payload = {
        "name": "path:{0}".format(label),
        "status": STATUS_AVAILABLE if exists else STATUS_UNAVAILABLE,
        "detail": "Path exists." if exists else "Path does not exist.",
        "path": str(path),
    }
    return payload


def probe_glob(directory: Optional[str], pattern: str, minimum: int = 1) -> Dict[str, Any]:
    """Check that *directory* contains at least *minimum* files matching *pattern*."""
    import glob

    label = os.path.basename(str(directory).rstrip("/\\")) if directory else "<unset>"
    check_name = "glob:{0}/{1}".format(label, pattern)
    if not directory or not os.path.isdir(directory):
        return {
            "name": check_name,
            "status": STATUS_UNAVAILABLE,
            "detail": "Directory is missing.",
            "path": str(directory) if directory else None,
        }
    matches = sorted(glob.glob(os.path.join(directory, pattern)))
    if len(matches) >= minimum:
        return {
            "name": check_name,
            "status": STATUS_AVAILABLE,
            "detail": "Found {0} matching file(s).".format(len(matches)),
            "count": len(matches),
            "path": str(directory),
        }
    return {
        "name": check_name,
        "status": STATUS_MISCONFIGURED,
        "detail": "Found {0} matching file(s), expected at least {1}.".format(len(matches), minimum),
        "count": len(matches),
        "path": str(directory),
    }


def probe_writable_dir(path: str) -> Dict[str, Any]:
    """Check that a directory is writable (creating it when necessary)."""
    from .paths import is_writable_dir

    label = os.path.basename(str(path).rstrip("/\\")) or str(path)
    if is_writable_dir(path):
        return {
            "name": "writable:{0}".format(label),
            "status": STATUS_AVAILABLE,
            "detail": "Directory is writable.",
            "path": str(path),
        }
    return {
        "name": "writable:{0}".format(label),
        "status": STATUS_MISCONFIGURED,
        "detail": "Directory is not writable.",
        "path": str(path),
    }


def python_environment() -> Dict[str, Any]:
    """Describe the running interpreter."""
    return {
        "python_version": sys.version.split()[0],
        "python_executable": sys.executable,
        "platform": sys.platform,
        "implementation": sys.implementation.name if hasattr(sys, "implementation") else None,
    }


def aggregate_status(checks: List[Dict[str, Any]], required_names: Optional[List[str]] = None) -> str:
    """Derive an overall status from individual probe results.

    When *required_names* is provided only those checks influence the result,
    which lets a module stay "available" although an unrelated optional tool is
    missing.
    """
    relevant = checks
    if required_names:
        wanted = set(required_names)
        relevant = [check for check in checks if check.get("name") in wanted] or checks
    if not relevant:
        return STATUS_NOT_REQUIRED
    statuses = [check.get("status") for check in relevant]
    if all(status == STATUS_NOT_REQUIRED for status in statuses):
        return STATUS_NOT_REQUIRED
    if any(status == STATUS_MISCONFIGURED for status in statuses):
        return STATUS_MISCONFIGURED
    if all(status == STATUS_AVAILABLE for status in statuses):
        return STATUS_AVAILABLE
    return STATUS_UNAVAILABLE
