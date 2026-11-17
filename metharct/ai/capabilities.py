# -*- coding: utf-8 -*-
"""Dynamic capability discovery for the MethArCT AI layer.

This module is the **single source of truth** for the module registry.  The CLI,
the MCP server, the serializers and the report builder all read from
:data:`MODULE_REGISTRY`; nothing else may hard-code a module list.

Availability is always derived from real dependency, database and model-file
checks -- never from the mere existence of a source directory.
"""

import os
import sys
import sysconfig
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .clock import utc_now_iso
from .constants import (
    AI_API_VERSION,
    SCHEMA_VERSION,
    STATUS_AVAILABLE,
    STATUS_MISCONFIGURED,
    STATUS_UNAVAILABLE,
)
from .errors import ErrorCode, MethArCTError
from .probes import (
    aggregate_status,
    probe_executable,
    probe_glob,
    probe_python_module,
    python_environment,
)


def _package_dir() -> str:
    """Return the directory containing the ``metharct`` package."""
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@dataclass(frozen=True)
class OptionSpec(object):
    """A module-specific option accepted through ``options``."""

    name: str
    type: str
    description: str
    default: Any = None
    choices: Tuple[str, ...] = ()
    minimum: Optional[float] = None
    maximum: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        payload = {
            "name": self.name,
            "type": self.type,
            "description": self.description,
            "default": self.default,
        }
        if self.choices:
            payload["choices"] = list(self.choices)
        if self.minimum is not None:
            payload["minimum"] = self.minimum
        if self.maximum is not None:
            payload["maximum"] = self.maximum
        return payload


@dataclass(frozen=True)
class ModuleCapability(object):
    """Static definition of one analysable MethArCT module."""

    name: str
    display_name: str
    description: str
    input_types: Tuple[str, ...]
    required_dependencies: Tuple[str, ...]
    optional_dependencies: Tuple[str, ...]
    output_artifacts: Tuple[str, ...]
    options: Tuple[OptionSpec, ...] = field(default_factory=tuple)
    #: Model identifiers reported in ``provenance.model_versions``.
    model_name: Optional[str] = None
    model_version: Optional[str] = None
    #: Scientific limitations copied verbatim into module results.
    limitations: Tuple[str, ...] = field(default_factory=tuple)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "display_name": self.display_name,
            "description": self.description,
            "input_types": list(self.input_types),
            "required_dependencies": list(self.required_dependencies),
            "optional_dependencies": list(self.optional_dependencies),
            "output_artifacts": list(self.output_artifacts),
            "options": [option.to_dict() for option in self.options],
            "model": {"name": self.model_name, "version": self.model_version},
            "limitations": list(self.limitations),
        }


_PREDICTION_LIMITATION = (
    "Values are model predictions derived from genome sequence composition, "
    "not experimental measurements."
)

MODULE_REGISTRY = (
    ModuleCapability(
        name="pathway",
        display_name="Metabolic pathway prediction",
        description=(
            "DIAMOND homology search of the query proteome against the methane, "
            "sulfur and nitrogen reference pathway databases, reporting pathway "
            "hits and completeness."
        ),
        input_types=("protein_fasta",),
        required_dependencies=("executable:diamond", "database:diamond"),
        optional_dependencies=(),
        output_artifacts=(
            "{prefix}_metabolic_summary.csv",
            "{prefix}_{database}_diamond.tsv",
        ),
        options=(
            OptionSpec("databases", "array[string]", "Restrict the search to specific databases."),
            OptionSpec("threads", "integer", "DIAMOND threads.", default=4, minimum=1, maximum=512),
            OptionSpec("evalue", "number", "DIAMOND E-value threshold.", default=1e-5, minimum=0.0),
            OptionSpec("max_target_seqs", "integer", "Maximum target sequences per query.", default=1, minimum=1),
        ),
        limitations=(
            "Homology detection is not proof of gene function or pathway activity.",
            "A pathway reported as not detected is not evidence that it is absent.",
        ),
    ),
    ModuleCapability(
        name="cultivation",
        display_name="Cultivability and amino-acid auxotrophy assessment",
        description=(
            "Assesses biosynthetic pathway completeness for amino acids, "
            "vitamins/coenzymes, nucleotides and ATPases to estimate cultivation "
            "requirements."
        ),
        input_types=("protein_fasta",),
        required_dependencies=("executable:diamond", "database:cultivation"),
        optional_dependencies=(),
        output_artifacts=("{prefix}_cultivation.csv", "{prefix}_amino_acid_biosynthesis.csv"),
        options=(
            OptionSpec("threads", "integer", "DIAMOND threads.", default=4, minimum=1, maximum=512),
            OptionSpec("evalue", "number", "DIAMOND E-value threshold.", default=1e-5, minimum=0.0),
            OptionSpec(
                "identity_threshold",
                "number",
                "Minimum percent identity for a reference hit.",
                default=30.0,
                minimum=0.0,
                maximum=100.0,
            ),
        ),
        limitations=(
            "Auxotrophy is inferred from missing biosynthetic pathway completeness "
            "and requires experimental confirmation.",
        ),
    ),
    ModuleCapability(
        name="temperature",
        display_name="Growth temperature range prediction",
        description=(
            "Predicts T_min, T_opt and T_max in degrees Celsius from amino acid "
            "and dipeptide composition using the embedded MLHuaiNanzi PLS models."
        ),
        input_types=("protein_fasta",),
        required_dependencies=("python_module:sklearn", "model:mlhuainanzi"),
        optional_dependencies=(),
        output_artifacts=("{prefix}_Huainanzi_Summary.tsv",),
        options=(),
        model_name="MLHuaiNanzi",
        limitations=(_PREDICTION_LIMITATION,),
    ),
    ModuleCapability(
        name="salinity",
        display_name="Salinity adaptation prediction",
        description=(
            "Predicts a salinity adaptation class from amino acid composition "
            "using the embedded SuSha ensemble model."
        ),
        input_types=("protein_fasta",),
        required_dependencies=("python_module:sklearn", "model:susha"),
        optional_dependencies=("python_module:shap", "python_module:openpyxl"),
        output_artifacts=("{prefix}_SuSha_Summary.tsv",),
        options=(),
        model_name="SuSha",
        limitations=(
            _PREDICTION_LIMITATION,
            "Class probabilities are model outputs and are not calibrated "
            "experimental confidence values.",
        ),
    ),
    ModuleCapability(
        name="ph",
        display_name="Growth pH preference prediction",
        description=(
            "Predicts optimum, maximum and minimum growth pH from genome features "
            "using the embedded GenomeSpot-style models."
        ),
        input_types=("protein_fasta",),
        required_dependencies=("python_module:hmmlearn", "model:ph_predictor"),
        optional_dependencies=(),
        output_artifacts=("{prefix}_pH_Summary.tsv",),
        options=(),
        model_name="pH predictor",
        limitations=(_PREDICTION_LIMITATION,),
    ),
    ModuleCapability(
        name="antibiotic",
        display_name="Antibiotic resistance prediction",
        description=(
            "Compares the query proteome with five methanogen reference proteomes "
            "by AAI and applies the published antibiotic selection rules."
        ),
        input_types=("protein_fasta",),
        required_dependencies=("executable:diamond", "database:kangshengsu"),
        optional_dependencies=(),
        output_artifacts=("{prefix}_antibiotic_selection.tsv",),
        options=(
            OptionSpec("threads", "integer", "DIAMOND threads.", default=4, minimum=1, maximum=512),
            OptionSpec("evalue", "number", "DIAMOND E-value threshold.", default=1e-5, minimum=0.0),
        ),
        limitations=(
            "Recommendations follow fixed AAI rules against a reference set of five "
            "methanogen genomes and are not a substitute for susceptibility testing.",
        ),
    ),
    ModuleCapability(
        name="genome_quality",
        display_name="Genome quality assessment (CheckM2)",
        description=(
            "Runs CheckM2 to estimate genome completeness and contamination, which "
            "conditions the reliability of every other prediction."
        ),
        input_types=("protein_fasta",),
        required_dependencies=("executable:checkm2",),
        optional_dependencies=(),
        output_artifacts=("{prefix}_checkm2_summary.csv",),
        options=(
            OptionSpec(
                "input_type",
                "string",
                "CheckM2 input type.",
                default="fasta",
                choices=("fasta", "directory"),
            ),
        ),
        model_name="CheckM2",
        limitations=(
            "CheckM2 requires its own reference database and is not bundled with MethArCT.",
        ),
    ),
)

MODULE_NAMES = tuple(module.name for module in MODULE_REGISTRY)
_MODULE_BY_NAME = {module.name: module for module in MODULE_REGISTRY}

#: Mapping of AI module identifiers onto the unified inspection ``sequence_type``.
_MODULE_INPUT_SEQUENCE_TYPES = {"protein_fasta": ("protein",)}


def get_module(name: str) -> ModuleCapability:
    """Return the definition of *name* or raise ``MODULE_UNKNOWN``."""
    module = _MODULE_BY_NAME.get(name)
    if module is None:
        raise MethArCTError(
            ErrorCode.MODULE_UNKNOWN,
            "Unknown module requested.",
            details={"module": name, "available_modules": list(MODULE_NAMES)},
        )
    return module


def normalize_module_list(modules: Optional[Sequence[str]]) -> List[str]:
    """Validate and de-duplicate a requested module list.

    ``None`` or an empty sequence means "all registered modules".
    """
    if not modules:
        return list(MODULE_NAMES)
    normalized: List[str] = []
    for raw in modules:
        if raw is None:
            continue
        name = str(raw).strip().lower()
        if not name:
            continue
        get_module(name)  # raises MODULE_UNKNOWN for typos
        if name not in normalized:
            normalized.append(name)
    if not normalized:
        return list(MODULE_NAMES)
    return normalized


def _default_base_dir() -> str:
    """Locate the reference database root without a Config object.

    Mirrors ``Config.resolve_database_base_dir`` so the AI layer and the
    analyzers always agree: environment variable, bundled package data, the
    repository checkout and ``share/metharct`` are all considered.
    """
    package_dir = _package_dir()
    candidates = [
        os.environ.get("METHARCT_DB_DIR"),
        os.path.join(package_dir, "data", "databases"),
        os.path.join(os.path.dirname(package_dir), "data", "databases"),
        os.path.join(sysconfig.get_paths().get("data", sys.prefix), "share", "metharct", "databases"),
    ]
    for candidate in candidates:
        if candidate and os.path.isdir(candidate):
            return candidate
    return candidates[1]


def _base_dir(config) -> str:
    """Return the database root for *config* (or for a default Config)."""
    if config is not None:
        configured = config.get("databases.base_dir")
        if configured:
            return configured
    try:
        from ..utils.config import Config

        resolved = Config().get("databases.base_dir")
        if resolved:
            return resolved
    except Exception:  # pragma: no cover - configuration is integral to MethArCT
        pass
    return _default_base_dir()


def _diamond_config(config) -> Dict[str, Any]:
    """Return the resolved DIAMOND command and database directory."""
    if config is None:
        return {
            "command": "diamond",
            "db_dir": os.path.join(_base_dir(None), "diamond"),
        }
    command = config.get("tools.diamond.path", "diamond")
    if config.get("tools.diamond.use_wsl", False):
        command = config.get("tools.diamond.wsl_path", "wsl diamond")
    return {"command": command, "db_dir": os.path.join(_base_dir(config), "diamond")}


def _database_dir(config, name: str) -> str:
    return os.path.join(_base_dir(config), name)


def _cultivation_dir(config, pathway_type: str) -> str:
    """Resolve one cultivation pathway directory (config-aware)."""
    if config is not None:
        try:
            resolved = config.get_cultivation_pathway_dir(pathway_type)
            if resolved:
                return resolved
        except Exception:  # pragma: no cover - fall back to defaults
            pass
    default_names = {
        "amino_acid": "path.aa",
        "vitamin": "path.vc",
        "nucleotide": "path_hesuan",
        "atp": "path_atp",
    }
    return os.path.join(
        _database_dir(config, "cultivation"), default_names.get(pathway_type, pathway_type)
    )


def _public_check(check: Dict[str, Any]) -> Dict[str, Any]:
    """Strip host paths from a probe result before exposing it to agents.

    ``doctor`` reports the full probe (including ``path`` / ``command``) because
    diagnosing an environment requires it; the capability manifest does not.
    """
    return {key: value for key, value in check.items() if key not in ("path", "command")}


def probe_module(module: ModuleCapability, config=None) -> Dict[str, Any]:
    """Probe one module and return its dynamic capability entry."""
    checks: List[Dict[str, Any]] = []
    package_dir = _package_dir()
    core_dir = os.path.join(package_dir, "core")

    if module.name in ("pathway", "cultivation", "antibiotic"):
        diamond = _diamond_config(config)
        checks.append(probe_executable(diamond["command"], ["--version"], timeout=15))
        if module.name == "pathway":
            checks.append(probe_glob(diamond["db_dir"], "*.dmnd"))
        elif module.name == "cultivation":
            for subdir in ("amino_acid", "vitamin", "nucleotide", "atp"):
                checks.append(
                    probe_glob(
                        _cultivation_dir(config, subdir),
                        "*.[fF][aA][sS][tT][aA]",
                    )
                )
        else:
            checks.append(probe_glob(_database_dir(config, "kangshengsu"), "*.faa", minimum=5))

    if module.name == "temperature":
        checks.append(probe_python_module("sklearn"))
        checks.append(probe_python_module("joblib"))
        checks.append(probe_glob(os.path.join(core_dir, "mlhuainanzi", "models"), "*.pkl", minimum=3))

    if module.name == "salinity":
        checks.append(probe_python_module("sklearn"))
        checks.append(probe_python_module("joblib"))
        checks.append(probe_python_module("shap"))
        checks.append(probe_glob(os.path.join(core_dir, "susha", "models"), "*.pkl"))

    if module.name == "ph":
        checks.append(probe_python_module("hmmlearn"))
        checks.append(probe_python_module("sklearn"))
        checks.append(probe_glob(os.path.join(core_dir, "ph_predictor", "models"), "*.joblib"))
        checks.append(probe_glob(os.path.join(core_dir, "ph_predictor", "hmm"), "*.joblib"))

    if module.name == "genome_quality":
        command = "checkm2"
        if config is not None:
            command = config.get("tools.checkm2.path", "checkm2")
            if config.get("tools.checkm2.use_wsl", False):
                command = config.get("tools.checkm2.wsl_path", "wsl checkm2")
        checks.append(probe_executable(command, ["--version"], timeout=20))

    # Optional dependencies never make a module unavailable.
    optional_names = set()
    for dependency in module.optional_dependencies:
        if dependency.startswith("python_module:"):
            optional_name = dependency.split(":", 1)[1]
            optional_names.add("python_module:{0}".format(optional_name))
            checks.append(probe_python_module(optional_name))

    required_checks = [check for check in checks if check["name"] not in optional_names]
    status = aggregate_status(required_checks)
    enabled = status == STATUS_AVAILABLE

    missing = [
        check
        for check in required_checks
        if check.get("status") in (STATUS_UNAVAILABLE, STATUS_MISCONFIGURED)
    ]

    payload = module.to_dict()
    payload.update(
        {
            "enabled": enabled,
            "status": status,
            "checks": [_public_check(check) for check in checks],
            "missing": [_public_check(check) for check in missing],
            "reasons": [check.get("detail") for check in missing if check.get("detail")],
        }
    )
    return payload


def compatible_modules(sequence_type: str) -> Tuple[List[str], List[str]]:
    """Return ``(compatible, incompatible)`` module names for a sequence type."""
    compatible: List[str] = []
    incompatible: List[str] = []
    for module in MODULE_REGISTRY:
        accepted = set()
        for input_type in module.input_types:
            accepted.update(_MODULE_INPUT_SEQUENCE_TYPES.get(input_type, ()))
        if sequence_type in accepted:
            compatible.append(module.name)
        else:
            incompatible.append(module.name)
    return compatible, incompatible


def module_input_types() -> List[str]:
    """Return the supported ``input_type`` values, in registration order."""
    types: List[str] = []
    for module in MODULE_REGISTRY:
        for input_type in module.input_types:
            if input_type not in types:
                types.append(input_type)
    return types


def _software_version() -> str:
    from .. import __version__  # lazy: avoids an import cycle at module load

    return __version__


def _mcp_available() -> bool:
    try:
        import importlib.util

        return importlib.util.find_spec("mcp") is not None
    except (ImportError, ValueError):  # pragma: no cover - defensive
        return False


def get_capabilities(config=None, include_details: bool = True) -> Dict[str, Any]:
    """Build the machine-readable capability description.

    All content is generated at call time from the installed package, the
    current configuration and real dependency probes.
    """
    modules = [probe_module(module, config) for module in MODULE_REGISTRY]
    if not include_details:
        for entry in modules:
            entry.pop("checks", None)

    environment = python_environment()
    # The capability manifest is agent-facing: keep the interpreter name, not
    # its absolute location (doctor keeps the full path for diagnosis).
    if environment.get("python_executable"):
        environment["python_executable"] = os.path.basename(environment["python_executable"])

    return {
        "schema_version": SCHEMA_VERSION,
        "generated_at": utc_now_iso(),
        "software": {
            "name": "MethArCT",
            "version": _software_version(),
            "ai_api_version": AI_API_VERSION,
            "commit": os.environ.get("METHARCT_COMMIT"),
        },
        "environment": environment,
        "modules": modules,
        "input_types": module_input_types(),
        "supports": {
            "async_jobs": True,
            "job_cancel": True,
            "batch_input": False,
            "mcp": _mcp_available(),
            "json_report": True,
        },
        "report_schema_version": SCHEMA_VERSION,
        "error_codes": _error_codes(),
    }


def _error_codes() -> List[str]:
    from .errors import ALL_ERROR_CODES

    return list(ALL_ERROR_CODES)
