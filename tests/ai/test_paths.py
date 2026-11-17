# -*- coding: utf-8 -*-
"""Tests for path safety helpers (AI layer)."""

import json
import os

import pytest

from metharct.ai.errors import ErrorCode, MethArCTError
from metharct.ai.paths import (
    atomic_write_json,
    atomic_write_text,
    enforce_allowed_roots,
    is_within,
    normalize,
    relative_posix,
    safe_join,
    sha256_file,
)


def test_normalize_makes_absolute(tmp_path):
    target = tmp_path / "sub" / "file.txt"
    resolved = normalize("sub/file.txt", base=str(tmp_path))
    assert resolved == target
    assert resolved.is_absolute()


def test_is_within_accepts_nested_and_rejects_siblings(tmp_path):
    root = tmp_path / "root"
    child = root / "a" / "b.txt"
    sibling = tmp_path / "other" / "b.txt"
    root.mkdir()
    assert is_within(child, root)
    assert is_within(root, root)
    assert not is_within(sibling, root)


def test_safe_join_rejects_parent_escape(tmp_path):
    root = tmp_path / "job"
    root.mkdir()
    assert safe_join(root, "artifacts", "a.tsv") == root / "artifacts" / "a.tsv"
    with pytest.raises(MethArCTError) as excinfo:
        safe_join(root, "..", "secret.txt")
    assert excinfo.value.code == ErrorCode.ARTIFACT_ACCESS_DENIED


def test_enforce_allowed_roots(tmp_path):
    allowed = tmp_path / "allowed"
    allowed.mkdir()
    inside = allowed / "file.faa"
    inside.write_text(">a\nMKT\n", encoding="utf-8")
    assert enforce_allowed_roots(inside, [str(allowed)]) == inside

    outside = tmp_path / "outside.faa"
    outside.write_text(">a\nMKT\n", encoding="utf-8")
    with pytest.raises(MethArCTError) as excinfo:
        enforce_allowed_roots(outside, [str(allowed)])
    assert excinfo.value.code == ErrorCode.ARTIFACT_ACCESS_DENIED


def test_enforce_allowed_roots_disabled_when_none(tmp_path):
    target = tmp_path / "anything.faa"
    target.write_text(">a\nMKT\n", encoding="utf-8")
    assert enforce_allowed_roots(target, None) == target


def test_atomic_writes_and_hash(tmp_path):
    text_path = tmp_path / "out.txt"
    atomic_write_text(text_path, "hello\n")
    assert text_path.read_text(encoding="utf-8") == "hello\n"

    json_path = tmp_path / "out.json"
    atomic_write_json(json_path, {"b": 1, "a": [2, 3]})
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert payload == {"b": 1, "a": [2, 3]}

    digest = sha256_file(text_path)
    assert len(digest) == 64 and all(char in "0123456789abcdef" for char in digest)


def test_relative_posix_uses_forward_slashes(tmp_path):
    root = tmp_path / "job"
    nested = root / "artifacts" / "x.tsv"
    nested.parent.mkdir(parents=True)
    nested.write_text("a\n", encoding="utf-8")
    assert relative_posix(nested, root) == "artifacts/x.tsv"
    assert os.sep not in relative_posix(nested, root)
