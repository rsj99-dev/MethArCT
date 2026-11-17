# -*- coding: utf-8 -*-
"""Packaging regression tests.

``setup.py`` feeds ``requirements.txt`` into ``install_requires``, so anything
listed there is installed for every end user. These tests keep development and
documentation tooling out of that file.
"""

import os
import re

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REQUIREMENTS = os.path.join(REPO_ROOT, "requirements.txt")
REQUIREMENTS_DEV = os.path.join(REPO_ROOT, "requirements-dev.txt")
SETUP_PY = os.path.join(REPO_ROOT, "setup.py")

#: Packages that must never be installed for end users.
DEV_ONLY = frozenset(
    [
        "pytest",
        "pytest-cov",
        "black",
        "flake8",
        "mypy",
        "sphinx",
        "sphinx-rtd-theme",
        "myst-parser",
        "jsonschema",
        "mcp",
    ]
)

_NAME_PATTERN = re.compile(r"^\s*([A-Za-z0-9._-]+)")


def _read_lines(path):
    with open(path, "r", encoding="utf-8") as stream:
        return [line.strip() for line in stream if line.strip()]


def _requirement_names(path):
    names = []
    for line in _read_lines(path):
        if line.startswith("#") or line.startswith("-"):
            continue
        match = _NAME_PATTERN.match(line.split(";")[0])
        if match:
            names.append(match.group(1).lower())
    return names


def test_requirements_txt_contains_no_dev_tooling():
    names = _requirement_names(REQUIREMENTS)
    assert names, "requirements.txt must list the runtime dependencies"
    offenders = sorted(DEV_ONLY.intersection(names))
    assert offenders == [], "dev-only packages must live in requirements-dev.txt: {0}".format(offenders)


def test_requirements_txt_lists_the_runtime_stack():
    names = set(_requirement_names(REQUIREMENTS))
    for package in ("biopython", "pandas", "numpy", "scikit-learn", "pyyaml", "joblib", "loguru"):
        assert package in names, package


def test_setup_py_keeps_reading_requirements_txt():
    with open(SETUP_PY, "r", encoding="utf-8") as stream:
        source = stream.read()
    assert "install_requires=read_requirements()" in source
    assert "extras_require=" in source
    assert '"ai"' in source and '"mcp"' in source


def test_mcp_extra_is_guarded_by_a_python_marker():
    with open(SETUP_PY, "r", encoding="utf-8") as stream:
        source = stream.read()
    # The MCP Python SDK requires Python >= 3.10; without the marker a 3.8 user
    # would get an unresolvable requirement instead of a clear runtime message.
    assert 'python_version >= "3.10"' in source


def test_requirements_dev_includes_runtime_and_tooling():
    lines = _read_lines(REQUIREMENTS_DEV)
    # Runtime packages are pulled in through the include line, never duplicated.
    assert any(line.replace(" ", "") == "-rrequirements.txt" for line in lines)
    names = set(_requirement_names(REQUIREMENTS_DEV))
    assert {"pytest", "jsonschema", "sphinx"}.issubset(names)
    runtime = set(_requirement_names(REQUIREMENTS))
    assert not runtime.intersection(names), "runtime packages must not be duplicated"
