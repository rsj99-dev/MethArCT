# -*- coding: utf-8 -*-
"""Construction and validation of ``metharct_ai_report.json``."""

import os
from typing import Any, Dict, List, Optional

from .capabilities import MODULE_REGISTRY
from .constants import AI_REPORT_FILENAME, SCHEMA_VERSION
from .errors import ErrorCode, MethArCTError
from .paths import atomic_write_json
from .probes import python_environment
from .serializers import json_safe
from .schema_check import validate_document

REPORT_SCHEMA_FILENAME = "metharct-ai-report.schema.json"
JOB_SCHEMA_FILENAME = "metharct-job.schema.json"
CAPABILITIES_SCHEMA_FILENAME = "metharct-capabilities.schema.json"
INPUT_SCHEMA_FILENAME = "metharct-input-inspection.schema.json"


def schemas_dir() -> str:
    """Directory holding the packaged JSON schemas."""
    return os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "schemas")


def load_schema(filename: str) -> Dict[str, Any]:
    """Load one packaged JSON schema document."""
    from .paths import read_json

    path = os.path.join(schemas_dir(), filename)
    schema = read_json(path)
    if not isinstance(schema, dict):
        raise MethArCTError(
            ErrorCode.SCHEMA_VALIDATION_FAILED,
            "Schema document is missing or unreadable.",
            details={"schema": filename, "path": path},
        )
    return schema


def validate_payload(payload: Any, schema_filename: str) -> List[str]:
    """Validate *payload* against a packaged schema and return error strings."""
    return validate_document(payload, load_schema(schema_filename))


def _software_block() -> Dict[str, Any]:
    from .. import __version__

    return {
        "name": "MethArCT",
        "version": __version__,
        "commit": os.environ.get("METHARCT_COMMIT"),
    }


def _provenance(config, module_results: Dict[str, Any]) -> Dict[str, Any]:
    """Collect the parameters and model identities that produced the results."""
    parameters: Dict[str, Any] = {}
    database_versions: Dict[str, Any] = {}
    model_versions: Dict[str, Any] = {}

    if config is not None:
        parameters["diamond"] = {
            "evalue": config.get("tools.diamond.evalue"),
            "max_target_seqs": config.get("tools.diamond.max_target_seqs"),
            "identity_threshold": config.get("tools.diamond.identity_threshold"),
            "threads": config.get("tools.diamond.threads"),
        }
        parameters["cultivation"] = {
            "evalue": config.get("cultivation.evalue"),
            "identity_threshold": config.get("cultivation.identity_threshold"),
            "threads": config.get("cultivation.threads"),
        }
        counts = config.get("reference_sequence_counts", {}) or {}
        database_versions["reference_sequence_counts"] = {
            str(key): value for key, value in counts.items()
        }
        database_versions["pathway_databases"] = {
            "methane": len(config.get("databases.methane_pathways", {}) or {}),
            "sulfur": len(config.get("databases.sulfur_pathways", {}) or {}),
            "nitrogen": len(config.get("databases.nitrogen_pathways", {}) or {}),
        }

    for module in MODULE_REGISTRY:
        if module.name not in module_results:
            continue
        if module.model_name:
            model_versions[module.name] = {
                "name": module.model_name,
                "version": module.model_version,
            }

    return {
        "parameters": json_safe(parameters),
        "database_versions": json_safe(database_versions),
        "model_versions": json_safe(model_versions),
        "environment": json_safe(python_environment()),
    }


def build_report(
    job: Dict[str, Any],
    input_info: Dict[str, Any],
    module_results: Dict[str, Any],
    artifacts: List[Dict[str, Any]],
    config=None,
    warnings: Optional[List[str]] = None,
    errors: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """Assemble the unified AI report for *job*.

    Missing values stay ``None``; nothing is filled with a placeholder number.
    """
    report = {
        "schema_version": SCHEMA_VERSION,
        "software": _software_block(),
        "run": {
            "job_id": job.get("job_id"),
            "created_at": job.get("created_at"),
            "started_at": job.get("started_at"),
            "finished_at": job.get("finished_at"),
            "status": job.get("status"),
            "duration_seconds": job.get("duration_seconds"),
        },
        "input": {
            "filename": input_info.get("filename") or job.get("input_filename"),
            "sha256": input_info.get("sha256") or job.get("input_sha256"),
            "format": input_info.get("format_detected") or "fasta",
            "type": input_info.get("sequence_type"),
            "sequence_count": input_info.get("sequence_count"),
            "total_residues": input_info.get("total_residues"),
        },
        "requested_modules": list(job.get("requested_modules") or []),
        "modules": module_results,
        "warnings": list(warnings if warnings is not None else job.get("warnings") or []),
        "errors": list(errors if errors is not None else job.get("errors") or []),
        "artifacts": artifacts,
        "provenance": _provenance(config, module_results),
    }
    return json_safe(report)


def write_report(job_dir: str, report: Dict[str, Any]) -> str:
    """Write the report atomically and return its absolute path."""
    target = os.path.join(job_dir, AI_REPORT_FILENAME)
    atomic_write_json(target, report)
    return target


def assert_report_valid(report: Dict[str, Any]) -> None:
    """Validate *report* against the packaged schema, raising on failure."""
    problems = validate_payload(report, REPORT_SCHEMA_FILENAME)
    if problems:
        raise MethArCTError(
            ErrorCode.SCHEMA_VALIDATION_FAILED,
            "The generated AI report does not satisfy its JSON schema.",
            details={"problems": problems[:10], "problem_count": len(problems)},
        )
