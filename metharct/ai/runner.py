# -*- coding: utf-8 -*-
"""Thin adapters that call the existing MethArCT analyzers.

The runner performs **parameter mapping only**: it never re-implements a
comparison, never retrains a model and never recomputes a prediction.  All
science stays in ``metharct.core``.
"""

import os
import time
import traceback
from typing import Any, Dict, List, Optional

from .capabilities import get_module, normalize_module_list
from .errors import ErrorCode, MethArCTError


def validate_options(module_name: str, options: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Validate module options against the registry and apply defaults.

    Unknown options and out-of-range values raise ``INVALID_ARGUMENT``; they are
    never silently ignored.
    """
    capability = get_module(module_name)
    allowed = {option.name: option for option in capability.options}
    supplied = dict(options or {})

    unknown = [key for key in supplied if key not in allowed]
    if unknown:
        raise MethArCTError(
            ErrorCode.INVALID_ARGUMENT,
            "Unknown option for module '{0}'.".format(module_name),
            details={"unknown_options": unknown, "allowed_options": sorted(allowed)},
        )

    resolved: Dict[str, Any] = {}
    for name, spec in allowed.items():
        value = supplied.get(name, spec.default)
        if value is None:
            resolved[name] = None
            continue
        if spec.type == "integer":
            if isinstance(value, bool) or not isinstance(value, int):
                raise MethArCTError(
                    ErrorCode.INVALID_ARGUMENT,
                    "Option '{0}' must be an integer.".format(name),
                    details={"module": module_name, "value": value},
                )
        elif spec.type == "number":
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise MethArCTError(
                    ErrorCode.INVALID_ARGUMENT,
                    "Option '{0}' must be a number.".format(name),
                    details={"module": module_name, "value": value},
                )
        elif spec.type == "array[string]":
            if not isinstance(value, (list, tuple)) or not all(isinstance(item, str) for item in value):
                raise MethArCTError(
                    ErrorCode.INVALID_ARGUMENT,
                    "Option '{0}' must be an array of strings.".format(name),
                    details={"module": module_name, "value": value},
                )
        elif spec.type == "string":
            if not isinstance(value, str):
                raise MethArCTError(
                    ErrorCode.INVALID_ARGUMENT,
                    "Option '{0}' must be a string.".format(name),
                    details={"module": module_name, "value": value},
                )
            if spec.choices and value not in spec.choices:
                raise MethArCTError(
                    ErrorCode.INVALID_ARGUMENT,
                    "Option '{0}' must be one of {1}.".format(name, sorted(spec.choices)),
                    details={"module": module_name, "value": value},
                )
        if spec.minimum is not None and isinstance(value, (int, float)) and value < spec.minimum:
            raise MethArCTError(
                ErrorCode.INVALID_ARGUMENT,
                "Option '{0}' is below the minimum of {1}.".format(name, spec.minimum),
                details={"module": module_name, "value": value},
            )
        if spec.maximum is not None and isinstance(value, (int, float)) and value > spec.maximum:
            raise MethArCTError(
                ErrorCode.INVALID_ARGUMENT,
                "Option '{0}' exceeds the maximum of {1}.".format(name, spec.maximum),
                details={"module": module_name, "value": value},
            )
        resolved[name] = value
    return resolved


def _apply_diamond_options(config, options: Dict[str, Any]) -> None:
    if options.get("threads") is not None:
        config.set("tools.diamond.threads", options["threads"])
    if options.get("evalue") is not None:
        config.set("tools.diamond.evalue", options["evalue"])
    if options.get("max_target_seqs") is not None:
        config.set("tools.diamond.max_target_seqs", options["max_target_seqs"])


def _run_pathway(config, input_file: str, output_prefix: str, options: Dict[str, Any]) -> Dict[str, Any]:
    from ..core.diamond_analyzer import DiamondAnalyzer

    _apply_diamond_options(config, options)
    analyzer = DiamondAnalyzer(config)
    databases = options.get("databases") or None
    native = analyzer.analyze_sequence(
        input_file=input_file,
        output_prefix=output_prefix,
        databases=databases,
    )
    native.setdefault("status", "success")
    return native


def _run_cultivation(config, input_file: str, output_prefix: str, options: Dict[str, Any]) -> Dict[str, Any]:
    from ..core.cultivation_analyzer import CultivationAnalyzer

    if options.get("threads") is not None:
        config.set("cultivation.threads", options["threads"])
    if options.get("evalue") is not None:
        config.set("cultivation.evalue", options["evalue"])
    if options.get("identity_threshold") is not None:
        config.set("cultivation.identity_threshold", options["identity_threshold"])

    analyzer = CultivationAnalyzer(config.config)
    pathway_results = analyzer.analyze_genome(genome_file=input_file)
    assessment = analyzer.generate_cultivability_assessment(pathway_results)
    amino_acids = analyzer.analyze_amino_acid_biosynthesis(pathway_results)

    cultivation_csv = "{0}_cultivation.csv".format(output_prefix)
    amino_acid_csv = "{0}_amino_acid_biosynthesis.csv".format(output_prefix)
    analyzer.export_results(pathway_results, assessment, cultivation_csv, "csv")
    analyzer.export_amino_acid_csv(amino_acids, amino_acid_csv)

    return {
        "status": "success",
        "pathway_results": pathway_results,
        "cultivability_assessment": assessment,
        "amino_acid_biosynthesis": amino_acids,
        "output_files": {"cultivation_csv": cultivation_csv, "amino_acid_csv": amino_acid_csv},
    }


def _run_temperature(config, input_file: str, output_prefix: str, options: Dict[str, Any]) -> Dict[str, Any]:
    from ..core.huainanzi_analyzer import HuainanziAnalyzer

    analyzer = HuainanziAnalyzer(config)
    if not getattr(analyzer, "tool_available", False):
        return {
            "status": "skipped",
            "message": "MLHuaiNanzi module is not importable in this environment.",
        }
    return analyzer.predict_temperature(input_file=input_file, output_prefix=output_prefix)


def _run_salinity(config, input_file: str, output_prefix: str, options: Dict[str, Any]) -> Dict[str, Any]:
    from ..core.susha_analyzer import SuShaAnalyzer

    analyzer = SuShaAnalyzer(config)
    if not getattr(analyzer, "tool_available", False):
        return {
            "status": "skipped",
            "message": "SuSha module is not importable in this environment.",
        }
    return analyzer.predict_salinity(input_file=input_file, output_prefix=output_prefix)


def _run_ph(config, input_file: str, output_prefix: str, options: Dict[str, Any]) -> Dict[str, Any]:
    from ..core.ph_analyzer import PHAnalyzer

    analyzer = PHAnalyzer(config)
    if not getattr(analyzer, "tool_available", False):
        return {
            "status": "skipped",
            "message": "pH predictor module is not importable in this environment.",
        }
    return analyzer.predict_ph(input_file=input_file, output_prefix=output_prefix)


def _run_antibiotic(config, input_file: str, output_prefix: str, options: Dict[str, Any]) -> Dict[str, Any]:
    from ..core.antibiotic_analyzer import AntibioticAnalyzer

    _apply_diamond_options(config, options)
    cpus = options.get("threads") or config.get("tools.diamond.threads", 4)
    evalue = options.get("evalue") or config.get("tools.diamond.evalue", 1e-5)
    analyzer = AntibioticAnalyzer(config=config, cpus=cpus, evalue=evalue)
    return analyzer.predict_antibiotics(input_file=input_file, output_prefix=output_prefix)


def _run_genome_quality(config, input_file: str, output_prefix: str, options: Dict[str, Any]) -> Dict[str, Any]:
    from ..core.checkm2_analyzer import CheckM2Analyzer

    analyzer = CheckM2Analyzer(config)
    if not getattr(analyzer, "tool_available", False):
        return {
            "status": "skipped",
            "message": "CheckM2 executable is not available in this environment.",
        }
    return analyzer.analyze_genome_quality(
        input_path=input_file,
        output_prefix=output_prefix,
        input_type=options.get("input_type") or "fasta",
    )


_RUNNERS = {
    "pathway": _run_pathway,
    "cultivation": _run_cultivation,
    "temperature": _run_temperature,
    "salinity": _run_salinity,
    "ph": _run_ph,
    "antibiotic": _run_antibiotic,
    "genome_quality": _run_genome_quality,
}


def execute_module(
    module_name: str,
    config,
    input_file: str,
    output_prefix: str,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Run one module and return ``(native_result)`` with timing metadata.

    Failures are captured and reported as ``status: failed`` so one failing
    module cannot hide the results of the others.
    """
    runner = _RUNNERS.get(module_name)
    if runner is None:
        raise MethArCTError(
            ErrorCode.MODULE_UNKNOWN,
            "No runner is registered for this module.",
            details={"module": module_name},
        )

    resolved = validate_options(module_name, options)
    started = time.time()
    try:
        native = runner(config, input_file, output_prefix, resolved)
    except Exception as exc:  # noqa: BLE001 - reported, not swallowed
        traceback.print_exc()
        native = {
            "status": "failed",
            "error": "{0}: {1}".format(type(exc).__name__, exc),
            "error_code": ErrorCode.ANALYSIS_FAILED,
        }
    elapsed = round(time.time() - started, 3)

    if not isinstance(native, dict):
        native = {
            "status": "failed",
            "error": "Analyzer returned an unsupported result type: {0}".format(
                type(native).__name__
            ),
            "error_code": ErrorCode.INTERNAL_ERROR,
        }
    native.setdefault("status", "success")
    native["runtime_seconds"] = elapsed
    native["options"] = resolved
    return native


def build_output_prefix(artifact_dir: str, input_file: str) -> str:
    """Build an absolute output prefix inside the job artifact directory."""
    stem = os.path.splitext(os.path.basename(input_file))[0] or "input"
    safe_stem = "".join(char if (char.isalnum() or char in "-_.") else "_" for char in stem)
    return os.path.join(artifact_dir, safe_stem)


def ordered_modules(modules: Optional[List[str]]) -> List[str]:
    """Return the requested modules in a deterministic execution order."""
    return normalize_module_list(modules)
