# -*- coding: utf-8 -*-
"""Path handling, hashing and atomic file writes for the AI layer.

Every user-supplied path goes through :func:`normalize` so that comparisons and
containment checks operate on absolute, case-normalised paths on both Windows
and POSIX systems.
"""

import hashlib
import json
import os
import tempfile
import time
from pathlib import Path
from typing import Any, Iterable, List, Optional, Union

from .errors import ErrorCode, MethArCTError

PathLike = Union[str, "os.PathLike[str]", Path]

#: Windows refuses ``os.replace`` while another process has the target open,
#: so atomic writes and reads retry briefly instead of failing spuriously.
_RETRY_ATTEMPTS = 10
_RETRY_DELAY_SECONDS = 0.02


def normalize(path: PathLike, base: Optional[PathLike] = None) -> Path:
    """Return an absolute path without resolving symlinks.

    Relative paths are resolved against *base* (default: current working
    directory).  ``~`` is expanded.
    """
    if path is None:
        raise MethArCTError(ErrorCode.INVALID_ARGUMENT, "Path must not be None.")
    text = os.fspath(path)
    if isinstance(text, bytes):  # pragma: no cover - defensive
        text = text.decode("utf-8", "replace")
    text = os.path.expanduser(text)
    if not os.path.isabs(text):
        root = os.fspath(base) if base is not None else os.getcwd()
        text = os.path.join(os.path.expanduser(str(root)), text)
    return Path(os.path.abspath(text))


def _compare_key(path: PathLike) -> str:
    """Normalised key for containment checks (case-insensitive on Windows)."""
    return os.path.normcase(os.path.realpath(os.fspath(path)))


def is_within(child: PathLike, root: PathLike) -> bool:
    """Return True when *child* is *root* or lives underneath it."""
    child_key = _compare_key(child)
    root_key = _compare_key(root)
    if child_key == root_key:
        return True
    if not root_key.endswith(os.sep):
        root_key = root_key + os.sep
    return child_key.startswith(root_key)


def normalize_roots(allowed_roots: Optional[Iterable[PathLike]]) -> Optional[List[Path]]:
    """Normalise an ``allowed_roots`` iterable, returning ``None`` when unset."""
    if allowed_roots is None:
        return None
    roots = []
    for root in allowed_roots:
        if root is None or str(root).strip() == "":
            continue
        roots.append(normalize(root))
    return roots or None


def enforce_allowed_roots(
    path: PathLike,
    allowed_roots: Optional[Iterable[PathLike]],
    code: str = ErrorCode.ARTIFACT_ACCESS_DENIED,
) -> Path:
    """Validate that *path* is inside one of *allowed_roots*.

    ``allowed_roots is None`` disables the restriction (the default for the
    plain Python API).  The MCP server always configures explicit roots.
    """
    resolved = normalize(path)
    roots = normalize_roots(allowed_roots)
    if roots is None:
        return resolved
    for root in roots:
        if is_within(resolved, root):
            return resolved
    raise MethArCTError(
        code,
        "Path is outside the configured allowed roots.",
        details={"path": str(resolved), "allowed_roots": [str(r) for r in roots]},
    )


def safe_join(root: PathLike, *parts: str) -> Path:
    """Join *parts* onto *root*, rejecting any escape from the root directory."""
    root_path = normalize(root)
    candidate = root_path
    for part in parts:
        part = str(part)
        if part in ("", "."):
            continue
        if os.path.isabs(part) or part.startswith("..") or "/.." in part or "\\.." in part:
            raise MethArCTError(
                ErrorCode.ARTIFACT_ACCESS_DENIED,
                "Refusing to resolve a path that escapes the job directory.",
                details={"part": part},
            )
        candidate = candidate / part
    if not is_within(candidate, root_path):
        raise MethArCTError(
            ErrorCode.ARTIFACT_ACCESS_DENIED,
            "Refusing to resolve a path that escapes the job directory.",
            details={"path": str(candidate)},
        )
    return candidate


def relative_posix(path: PathLike, base: PathLike) -> str:
    """Return *path* relative to *base* using forward slashes."""
    rel = os.path.relpath(os.fspath(normalize(path)), os.fspath(normalize(base)))
    return rel.replace(os.sep, "/")


def sha256_file(path: PathLike, chunk_size: int = 1024 * 1024) -> str:
    """Compute the SHA-256 digest of a file using chunked reads."""
    digest = hashlib.sha256()
    with open(os.fspath(path), "rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    """Compute the SHA-256 digest of a UTF-8 string."""
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _replace_with_retry(source: str, target: str) -> None:
    """Atomically replace *target* with *source*, retrying transient locks."""
    last_error = None
    for attempt in range(_RETRY_ATTEMPTS):
        try:
            os.replace(source, target)
            return
        except PermissionError as exc:  # Windows: reader holds the file open
            last_error = exc
        except OSError as exc:  # pragma: no cover - unusual filesystems
            last_error = exc
        time.sleep(_RETRY_DELAY_SECONDS * (attempt + 1))
    if last_error is not None:
        raise last_error


def atomic_write_text(path: PathLike, text: str, encoding: str = "utf-8") -> Path:
    """Write *text* atomically (temporary file + ``os.replace``)."""
    target = normalize(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    handle, temp_path = tempfile.mkstemp(
        prefix=".tmp_", suffix=target.suffix or ".tmp", dir=str(target.parent)
    )
    try:
        with os.fdopen(handle, "w", encoding=encoding, newline="") as stream:
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        _replace_with_retry(temp_path, str(target))
    except Exception:
        try:
            os.unlink(temp_path)
        except OSError:
            pass
        raise
    return target


def atomic_write_json(path: PathLike, payload: Any, indent: int = 2) -> Path:
    """Serialise *payload* to JSON and write it atomically."""
    text = json.dumps(payload, indent=indent, ensure_ascii=False, sort_keys=False)
    return atomic_write_text(path, text + "\n")


def read_json(path: PathLike, default: Any = None) -> Any:
    """Read a JSON file, returning *default* when it is missing or unreadable.

    Transient sharing violations (a concurrent atomic replace on Windows) are
    retried before giving up, so a polling reader never sees a spurious
    "missing" file.
    """
    target = os.fspath(path)
    for attempt in range(_RETRY_ATTEMPTS):
        try:
            with open(target, "r", encoding="utf-8") as stream:
                return json.load(stream)
        except FileNotFoundError:
            return default
        except (IOError, OSError, ValueError):
            time.sleep(_RETRY_DELAY_SECONDS * (attempt + 1))
    return default


def ensure_dir(path: PathLike) -> Path:
    """Create *path* (and parents) if needed and return it."""
    resolved = normalize(path)
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved


def is_writable_dir(path: PathLike) -> bool:
    """Return True when *path* exists (or can be created) and is writable."""
    try:
        resolved = normalize(path)
        resolved.mkdir(parents=True, exist_ok=True)
        probe = resolved / ".metharct_write_probe"
        with open(str(probe), "w", encoding="utf-8") as stream:
            stream.write("probe")
        os.unlink(str(probe))
        return True
    except (IOError, OSError):
        return False
