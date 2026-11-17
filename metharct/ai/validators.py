# -*- coding: utf-8 -*-
"""Input file inspection and validation.

The inspection is deliberately streaming: multi-hundred-megabyte proteomes are
never loaded into memory as a whole, and sequence type detection uses explicit
alphabet evidence with a rule-level confidence label (never a "calibrated
probability").
"""

import os
import statistics
import unicodedata
from typing import Any, Dict, List, Optional

from .constants import DEFAULT_MAX_FILE_SIZE_MB, SCHEMA_VERSION
from .errors import ErrorCode, MethArCTError
from .paths import enforce_allowed_roots, normalize, sha256_file

#: Canonical amino-acid alphabet plus ambiguity/extended codes.
PROTEIN_LETTERS = frozenset("ACDEFGHIKLMNPQRSTVWYBXZJUO")

#: IUPAC nucleotide alphabet (DNA + RNA + ambiguity codes).
NUCLEOTIDE_LETTERS = frozenset("ACGTUNRYSWKMBDHV")

#: Letters that can only be amino acids (never nucleotide ambiguity codes).
PROTEIN_ONLY_LETTERS = PROTEIN_LETTERS - NUCLEOTIDE_LETTERS

#: Characters tolerated inside sequence lines without being counted as residues.
IGNORED_SEQUENCE_CHARS = frozenset("*-.~")

#: Text encodings tried in order when reading a FASTA file.
_ENCODINGS = ("utf-8-sig", "utf-16", "latin-1")

SUPPORTED_SEQUENCE_TYPES = ("protein", "nucleotide", "ambiguous", "unknown")


def _open_text(path: str):
    """Open *path* for text reading, trying common encodings."""
    last_error = None
    for encoding in _ENCODINGS:
        try:
            stream = open(path, "r", encoding=encoding, newline="")
            # Force a decode of the first block so encoding errors surface here.
            stream.read(0)
        except (UnicodeError, LookupError) as exc:  # pragma: no cover - rare
            last_error = exc
            try:
                stream.close()
            except Exception:
                pass
            continue
        return stream
    raise MethArCTError(
        ErrorCode.INPUT_FORMAT_INVALID,
        "Unable to decode the input file as text.",
        details={"path": path, "last_error": str(last_error)},
    )


def _looks_binary(path: str, sample_size: int = 4096) -> bool:
    """Return True when the file contains NUL bytes in its first block."""
    with open(path, "rb") as handle:
        sample = handle.read(sample_size)
    return b"\x00" in sample and not sample.startswith(b"\xff\xfe") and not sample.startswith(b"\xfe\xff")


def _is_probably_compressed(path: str) -> bool:
    """Detect gzip / bzip2 / xz / zip magic numbers."""
    magic = (b"\x1f\x8b", b"BZh", b"\xfd7zXZ", b"PK\x03\x04")
    with open(path, "rb") as handle:
        head = handle.read(4)
    return any(head.startswith(prefix) for prefix in magic)


def _classify_sequence_type(protein_only: int, nucleotide_letters: int, total_letters: int) -> Dict[str, Any]:
    """Classify the sequence alphabet from residue counts.

    Returns ``{"sequence_type": ..., "sequence_type_confidence": ...}``.
    Confidence is a rule level (``high``/``medium``/``low``), never a
    probabilistically calibrated score.
    """
    if total_letters <= 0:
        return {"sequence_type": "unknown", "sequence_type_confidence": "none"}

    protein_only_fraction = protein_only / float(total_letters)
    nucleotide_fraction = nucleotide_letters / float(total_letters)

    if protein_only_fraction >= 0.10:
        return {"sequence_type": "protein", "sequence_type_confidence": "high"}
    if protein_only_fraction >= 0.02:
        return {"sequence_type": "protein", "sequence_type_confidence": "medium"}
    if nucleotide_fraction >= 0.98:
        if total_letters >= 200:
            confidence = "high"
        elif total_letters >= 30:
            confidence = "medium"
        else:
            confidence = "low"
        return {"sequence_type": "nucleotide", "sequence_type_confidence": confidence}
    if nucleotide_fraction >= 0.90:
        return {"sequence_type": "nucleotide", "sequence_type_confidence": "low"}
    if protein_only_fraction > 0:
        return {"sequence_type": "protein", "sequence_type_confidence": "low"}
    return {"sequence_type": "ambiguous", "sequence_type_confidence": "low"}


def _median_int(values: List[int]) -> Optional[float]:
    if not values:
        return None
    return float(statistics.median(values))


def inspect_fasta(
    path: str,
    expected_type: str = "auto",
    max_file_size_mb: float = DEFAULT_MAX_FILE_SIZE_MB,
    allowed_roots: Optional[List[str]] = None,
    compute_sha256: bool = True,
) -> Dict[str, Any]:
    """Inspect a FASTA file and report structural and composition statistics.

    Args:
        path: User supplied input path.
        expected_type: ``auto`` / ``protein_fasta`` / ``nucleotide_fasta``.
            ``protein`` and ``nucleotide`` are accepted aliases.
        max_file_size_mb: Upper bound for the file size in megabytes
            (``<= 0`` disables the limit).
        allowed_roots: Optional whitelist of directories the path must live in.
        compute_sha256: Compute and return the digest (needed for provenance).

    Returns:
        A JSON-serialisable inspection dictionary.

    Raises:
        MethArCTError: for missing/unreadable/oversized/invalid input.
    """
    resolved = enforce_allowed_roots(normalize(path), allowed_roots, code=ErrorCode.INPUT_FILE_UNREADABLE)
    target = str(resolved)

    warnings: List[str] = []
    errors: List[str] = []

    if not resolved.exists():
        raise MethArCTError(
            ErrorCode.INPUT_FILE_NOT_FOUND,
            "Input file does not exist.",
            details={"path": target},
        )
    if not resolved.is_file():
        raise MethArCTError(
            ErrorCode.INPUT_FORMAT_INVALID,
            "Input path is not a regular file.",
            details={"path": target},
        )
    if not os.access(target, os.R_OK):
        raise MethArCTError(
            ErrorCode.INPUT_FILE_UNREADABLE,
            "Input file is not readable.",
            details={"path": target},
        )

    size_bytes = os.path.getsize(target)
    max_bytes = 0
    if max_file_size_mb and float(max_file_size_mb) > 0:
        max_bytes = int(float(max_file_size_mb) * 1024 * 1024)
    if max_bytes > 0 and size_bytes > max_bytes:
        raise MethArCTError(
            ErrorCode.INPUT_FILE_TOO_LARGE,
            "Input file exceeds the configured size limit.",
            details={"path": target, "size_bytes": size_bytes, "max_file_size_mb": max_file_size_mb},
        )
    if size_bytes == 0:
        raise MethArCTError(
            ErrorCode.INPUT_FORMAT_INVALID,
            "Input file is empty.",
            details={"path": target},
        )

    if _is_probably_compressed(target):
        raise MethArCTError(
            ErrorCode.INPUT_FORMAT_INVALID,
            "Compressed inputs are not supported; decompress the file first.",
            details={"path": target},
        )
    if _looks_binary(target):
        raise MethArCTError(
            ErrorCode.INPUT_FORMAT_INVALID,
            "Input file appears to be binary, not FASTA text.",
            details={"path": target},
        )

    sequence_count = 0
    total_letters = 0
    lengths: List[int] = []
    ids = set()
    duplicate_ids = 0
    invalid_sequence_count = 0
    invalid_char_samples: List[str] = []
    protein_only_count = 0
    nucleotide_letter_count = 0
    other_letter_count = 0
    current_length = 0
    current_has_invalid = False
    header_seen = False

    stream = _open_text(target)
    try:
        for raw_line in stream:
            line = raw_line.rstrip("\r\n")
            if not line.strip():
                continue
            if line.startswith(">"):
                if sequence_count > 0:
                    lengths.append(current_length)
                    if current_has_invalid:
                        invalid_sequence_count += 1
                header_seen = True
                sequence_count += 1
                current_length = 0
                current_has_invalid = False

                header = line[1:].strip()
                if not header:
                    warnings.append("Empty FASTA header at sequence {0}.".format(sequence_count))
                else:
                    seq_id = header.split()[0]
                    if seq_id in ids:
                        duplicate_ids += 1
                    else:
                        ids.add(seq_id)
                continue

            if not header_seen:
                raise MethArCTError(
                    ErrorCode.INPUT_FORMAT_INVALID,
                    "The first non-empty line of a FASTA file must be a header "
                    "starting with '>'.",
                    details={"path": target},
                )

            for char in line:
                if char.isspace():
                    continue
                upper = char.upper()
                if upper in PROTEIN_LETTERS:
                    current_length += 1
                    total_letters += 1
                    if upper in PROTEIN_ONLY_LETTERS:
                        protein_only_count += 1
                    elif upper in NUCLEOTIDE_LETTERS:
                        nucleotide_letter_count += 1
                elif upper in IGNORED_SEQUENCE_CHARS:
                    continue
                else:
                    current_has_invalid = True
                    if len(invalid_char_samples) < 5:
                        name = unicodedata.name(char, "U+{0:04X}".format(ord(char)))
                        if name not in invalid_char_samples:
                            invalid_char_samples.append(name)
                    other_letter_count += 1
        if sequence_count > 0:
            lengths.append(current_length)
            if current_has_invalid:
                invalid_sequence_count += 1
    finally:
        stream.close()

    if sequence_count == 0:
        raise MethArCTError(
            ErrorCode.INPUT_FORMAT_INVALID,
            "No FASTA records were found in the input file.",
            details={"path": target},
        )

    classification = _classify_sequence_type(
        protein_only_count, nucleotide_letter_count, total_letters
    )
    sequence_type = classification["sequence_type"]
    confidence = classification["sequence_type_confidence"]

    if invalid_sequence_count:
        errors.append(
            "{0} sequence(s) contain characters outside the protein/nucleotide alphabets.".format(
                invalid_sequence_count
            )
        )
    if duplicate_ids:
        warnings.append("{0} duplicate sequence identifier(s) detected.".format(duplicate_ids))
    if other_letter_count:
        warnings.append(
            "{0} residue character(s) are not standard IUPAC letters ({1}).".format(
                other_letter_count, ", ".join(invalid_char_samples)
            )
        )
    if sequence_type == "ambiguous":
        warnings.append(
            "Sequence type could not be determined from the alphabet; specify it explicitly."
        )

    expected = (expected_type or "auto").strip().lower()
    expected_map = {
        "auto": "auto",
        "protein": "protein",
        "protein_fasta": "protein",
        "nucleotide": "nucleotide",
        "nucleotide_fasta": "nucleotide",
    }
    if expected not in expected_map:
        raise MethArCTError(
            ErrorCode.INVALID_ARGUMENT,
            "Unsupported expected_type value.",
            details={"expected_type": expected_type, "allowed": sorted(expected_map)},
        )
    expected_resolved = expected_map[expected]

    if expected_resolved != "auto":
        if sequence_type == "ambiguous":
            errors.append(
                "Expected {0} input but the sequence alphabet is ambiguous.".format(expected_resolved)
            )
        elif sequence_type != expected_resolved:
            errors.append(
                "Expected {0} input but detected {1}.".format(expected_resolved, sequence_type)
            )

    if compute_sha256:
        sha256 = sha256_file(target)
    else:  # pragma: no cover - provenance always wanted in practice
        sha256 = None

    valid = not errors

    return {
        "schema_version": SCHEMA_VERSION,
        "path": target,
        "filename": resolved.name,
        "file_exists": True,
        "file_size_bytes": size_bytes,
        "format_detected": "fasta",
        "sequence_type": sequence_type,
        "sequence_type_confidence": confidence,
        "expected_type": expected_resolved,
        "sequence_count": sequence_count,
        "total_residues": total_letters,
        "length": {
            "min": int(min(lengths)) if lengths else None,
            "median": _median_int(lengths),
            "max": int(max(lengths)) if lengths else None,
        },
        "invalid_sequence_count": invalid_sequence_count,
        "duplicate_id_count": duplicate_ids,
        "valid": valid,
        "compatible_modules": [],
        "incompatible_modules": [],
        "warnings": warnings,
        "errors": errors,
        "sha256": sha256,
        "modified_time": _iso_mtime(target),
    }


def _iso_mtime(path: str) -> Optional[str]:
    """Return the file modification time as an ISO-8601 UTC string."""
    from .clock import timestamp_to_iso

    try:
        stamp = os.path.getmtime(path)
    except OSError:  # pragma: no cover - defensive
        return None
    return timestamp_to_iso(stamp)


def assert_sequence_type_supported(inspection: Dict[str, Any], module: str, supported_types) -> None:
    """Raise when the inspected input type is not accepted by *module*."""
    sequence_type = inspection.get("sequence_type")
    supported = set(supported_types)
    if sequence_type == "ambiguous":
        raise MethArCTError(
            ErrorCode.INPUT_TYPE_AMBIGUOUS,
            "Input sequence type is ambiguous; specify expected_type explicitly.",
            details={"module": module, "sequence_type": sequence_type},
        )
    if sequence_type not in supported:
        raise MethArCTError(
            ErrorCode.INPUT_MODULE_INCOMPATIBLE,
            "Input sequence type is not compatible with the requested module.",
            details={
                "module": module,
                "sequence_type": sequence_type,
                "supported_types": sorted(supported),
            },
        )
