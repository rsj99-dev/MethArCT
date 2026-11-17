# -*- coding: utf-8 -*-
"""Environment self-check (``metharct doctor``).

Doctor only *reports* the state of the environment.  It never installs
packages, downloads databases or modifies ``PATH``.
"""

import os
import tempfile
from typing import Any, Dict, List, Optional

from .capabilities import get_capabilities, _cultivation_dir, _database_dir, _diamond_config
from .clock import utc_now_iso
from .constants import (
    SCHEMA_VERSION,
    STATUS_AVAILABLE,
    STATUS_MISCONFIGURED,
    STATUS_NOT_REQUIRED,
    STATUS_UNAVAILABLE,
)
from .paths import is_writable_dir
from .probes import (
    probe_executable,
    probe_glob,
    probe_package_version,
    probe_python_module,
    python_environment,
)

#: Packages MethArCT imports unconditionally at analysis time.
_CORE_PACKAGES = ("biopython", "pandas", "numpy", "scikit-learn", "pyyaml", "joblib")

#: Python modules checked individually (import name -> distribution name).
_CORE_MODULES = {
    "Bio": "biopython",
    "pandas": "pandas",
    "numpy": "numpy",
    "yaml": "pyyaml",
    "joblib": "joblib",
}


def _check(
    check_id: str,
    category: str,
    status: str,
    detail: str,
    required: bool = False,
    hint: Optional[str] = None,
    version: Optional[str] = None,
    extra: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    payload = {
        "id": check_id,
        "category": category,
        "status": status,
        "required": required,
        "detail": detail,
        "hint": hint,
        "version": version,
    }
    if extra:
        payload.update(extra)
    return payload


def _python_version_check() -> Dict[str, Any]:
    import sys

    major, minor = sys.version_info[0], sys.version_info[1]
    supported = (major, minor) >= (3, 8)
    return _check(
        "python_version",
        "python",
        STATUS_AVAILABLE if supported else STATUS_MISCONFIGURED,
        "Python {0}.{1}.{2}".format(major, minor, sys.version_info[2]),
        required=True,
        hint=None if supported else "MethArCT requires Python 3.8 or newer.",
    )


def _package_checks() -> List[Dict[str, Any]]:
    checks = []
    for distribution in _CORE_PACKAGES:
        result = probe_package_version(distribution)
        checks.append(
            _check(
                "package:{0}".format(distribution),
                "python_packages",
                result["status"],
                result["detail"],
                required=True,
                hint="pip install {0}".format(distribution) if result["status"] != STATUS_AVAILABLE else None,
                version=result.get("version"),
            )
        )
    for module_name, distribution in _CORE_MODULES.items():
        result = probe_python_module(module_name)
        if result["status"] != STATUS_AVAILABLE:
            checks.append(
                _check(
                    "import:{0}".format(module_name),
                    "python_packages",
                    result["status"],
                    result["detail"],
                    required=True,
                    hint="pip install {0}".format(distribution),
                )
            )
    return checks


def _check_diamond(config) -> List[Dict[str, Any]]:
    diamond = _diamond_config(config)
    executable = probe_executable(diamond["command"], ["--version"], timeout=20)
    checks = [
        _check(
            "executable:diamond",
            "external_tools",
            executable["status"],
            executable["detail"],
            required=True,
            hint=executable.get("hint")
            or "Install DIAMOND and set tools.diamond.path in metharct_config.yaml.",
            version=executable.get("version"),
        )
    ]

    db_dir = diamond["db_dir"]
    dmnd = probe_glob(db_dir, "*.dmnd")
    checks.append(
        _check(
            "database:diamond",
            "databases",
            dmnd["status"],
            "{0} (directory: {1})".format(dmnd["detail"], db_dir),
            required=True,
            hint="Run the pathway module once with a working DIAMOND to build the .dmnd files.",
            extra={"count": dmnd.get("count", 0), "path": db_dir},
        )
    )
    return checks


def _check_pathway_references(config) -> List[Dict[str, Any]]:
    from .capabilities import MODULE_REGISTRY

    pathway = next(module for module in MODULE_REGISTRY if module.name == "pathway")
    base_dir = _database_dir(config, "")
    checks = []
    for subdir in ("methane", "sulfur", "nitrogen"):
        result = probe_glob(os.path.join(base_dir, subdir), "*.fasta")
        checks.append(
            _check(
                "database:{0}".format(subdir),
                "databases",
                result["status"],
                "{0} (directory: {1})".format(result["detail"], os.path.join(base_dir, subdir)),
                required=False,
                hint="Reference sequences are required for {0} pathway detection.".format(subdir),
                extra={"count": result.get("count", 0)},
            )
        )
    # The module entry is only used to keep the check list in sync with the registry.
    checks.append(
        _check(
            "registry:pathway",
            "registry",
            STATUS_AVAILABLE if pathway else STATUS_MISCONFIGURED,
            "Module 'pathway' is registered.",
            required=False,
        )
    )
    return checks


def _check_cultivation(config) -> List[Dict[str, Any]]:
    checks = []
    for pathway_type in ("amino_acid", "vitamin", "nucleotide", "atp"):
        directory = _cultivation_dir(config, pathway_type)
        result = probe_glob(directory, "*.[fF][aA][sS][tT][aA]")
        checks.append(
            _check(
                "database:cultivation:{0}".format(pathway_type),
                "databases",
                result["status"],
                "{0} (directory: {1})".format(result["detail"], directory),
                required=False,
                extra={"count": result.get("count", 0), "path": directory},
            )
        )
    return checks


def _check_kangshengsu(config) -> List[Dict[str, Any]]:
    directory = _database_dir(config, "kangshengsu")
    result = probe_glob(directory, "*.faa", minimum=5)
    return [
        _check(
            "database:kangshengsu",
            "databases",
            result["status"],
            "{0} (directory: {1})".format(result["detail"], directory),
            required=False,
            extra={"count": result.get("count", 0), "path": directory},
        )
    ]


def _check_checkm2(config) -> List[Dict[str, Any]]:
    command = "checkm2"
    if config is not None:
        command = config.get("tools.checkm2.path", "checkm2")
        if config.get("tools.checkm2.use_wsl", False):
            command = config.get("tools.checkm2.wsl_path", "wsl checkm2")
    result = probe_executable(command, ["--version"], timeout=25)
    checks = [
        _check(
            "executable:checkm2",
            "external_tools",
            result["status"] if result["status"] != STATUS_UNAVAILABLE else STATUS_NOT_REQUIRED,
            result["detail"],
            required=False,
            hint="Optional: install CheckM2 to enable the genome_quality module.",
            version=result.get("version"),
        )
    ]

    database_path = os.environ.get("CHECKM2DB")
    if database_path:
        exists = os.path.isdir(database_path)
        checks.append(
            _check(
                "database:checkm2",
                "databases",
                STATUS_AVAILABLE if exists else STATUS_MISCONFIGURED,
                "CHECKM2DB={0}".format(database_path),
                required=False,
                extra={"path": database_path},
            )
        )
    else:
        checks.append(
            _check(
                "database:checkm2",
                "databases",
                STATUS_UNAVAILABLE,
                "CHECKM2DB is not set and no CheckM2 database path was detected.",
                required=False,
                hint="Set CHECKM2DB to the directory of the downloaded CheckM2 database.",
            )
        )
    return checks


def _check_models() -> List[Dict[str, Any]]:
    package_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    core_dir = os.path.join(package_dir, "core")
    checks = []
    model_specs = (
        ("mlhuainanzi", os.path.join(core_dir, "mlhuainanzi", "models"), "*.pkl", 3),
        ("susha", os.path.join(core_dir, "susha", "models"), "*.pkl", 1),
        ("ph_models", os.path.join(core_dir, "ph_predictor", "models"), "*.joblib", 3),
        ("ph_hmm", os.path.join(core_dir, "ph_predictor", "hmm"), "*.joblib", 1),
    )
    for name, directory, pattern, minimum in model_specs:
        result = probe_glob(directory, pattern, minimum=minimum)
        checks.append(
            _check(
                "model:{0}".format(name),
                "models",
                result["status"],
                "{0} (directory: {1})".format(result["detail"], directory),
                required=False,
                hint="Reinstall MethArCT; embedded model files are shipped with the package.",
                extra={"count": result.get("count", 0), "path": directory},
            )
        )
    return checks


def _check_mcp_sdk() -> List[Dict[str, Any]]:
    """Report whether the optional MCP server can run on this interpreter."""
    import sys

    if sys.version_info < (3, 10):
        return [
            _check(
                "python_module:mcp",
                "optional_features",
                STATUS_NOT_REQUIRED,
                "MCP server unavailable: the MCP SDK requires Python 3.10+ and this "
                "interpreter is {0}.{1}.{2}.".format(*sys.version_info[:3]),
                required=False,
                hint="Use Python 3.10 or newer and run: pip install \"metharct[mcp]\"",
            )
        ]

    result = probe_python_module("mcp")
    available = result["status"] == STATUS_AVAILABLE
    return [
        _check(
            "python_module:mcp",
            "optional_features",
            STATUS_AVAILABLE if available else STATUS_NOT_REQUIRED,
            result["detail"],
            required=False,
            hint=None if available else 'pip install "metharct[mcp]"',
            version=result.get("version"),
        )
    ]


def _check_filesystem(workspace: Optional[str]) -> List[Dict[str, Any]]:
    checks = []
    temp_dir = tempfile.gettempdir()
    checks.append(
        _check(
            "filesystem:temp",
            "filesystem",
            STATUS_AVAILABLE if is_writable_dir(temp_dir) else STATUS_MISCONFIGURED,
            "Temporary directory: {0}".format(temp_dir),
            required=True,
            hint="Set TMPDIR/TEMP to a writable directory.",
        )
    )
    if workspace:
        writable = is_writable_dir(workspace)
        checks.append(
            _check(
                "filesystem:workspace",
                "filesystem",
                STATUS_AVAILABLE if writable else STATUS_MISCONFIGURED,
                "Workspace: {0}".format(workspace),
                required=True,
                hint="Choose a workspace directory the current user can write to.",
            )
        )
    return checks


def run_doctor(
    config=None,
    workspace: Optional[str] = None,
    check_databases: bool = True,
) -> Dict[str, Any]:
    """Run all environment checks and return a structured report."""
    checks: List[Dict[str, Any]] = []
    checks.append(_python_version_check())
    checks.extend(_package_checks())
    checks.extend(_check_diamond(config))
    if check_databases:
        checks.extend(_check_pathway_references(config))
        checks.extend(_check_cultivation(config))
        checks.extend(_check_kangshengsu(config))
        checks.extend(_check_checkm2(config))
        checks.extend(_check_models())
    checks.extend(_check_filesystem(workspace))
    checks.extend(_check_mcp_sdk())

    summary = {
        STATUS_AVAILABLE: 0,
        STATUS_UNAVAILABLE: 0,
        STATUS_MISCONFIGURED: 0,
        STATUS_NOT_REQUIRED: 0,
    }
    for check in checks:
        summary[check["status"]] = summary.get(check["status"], 0) + 1

    blocking = [
        check
        for check in checks
        if check["required"] and check["status"] in (STATUS_UNAVAILABLE, STATUS_MISCONFIGURED)
    ]

    capabilities = get_capabilities(config, include_details=False)

    from .. import __version__

    return {
        "schema_version": SCHEMA_VERSION,
        "generated_at": utc_now_iso(),
        "software": {
            "name": "MethArCT",
            "version": __version__,
            "commit": os.environ.get("METHARCT_COMMIT"),
        },
        "environment": python_environment(),
        "ok": not blocking,
        "blocking_issues": [
            {"id": check["id"], "status": check["status"], "detail": check["detail"]}
            for check in blocking
        ],
        "summary": summary,
        "checks": checks,
        "modules": [
            {
                "name": module["name"],
                "enabled": module["enabled"],
                "status": module["status"],
                "reasons": module.get("reasons", []),
            }
            for module in capabilities["modules"]
        ],
    }
