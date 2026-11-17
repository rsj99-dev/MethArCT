# -*- coding: utf-8 -*-
"""Normalisation of native MethArCT results into the unified AI report shape.

Serializers only *map* existing results; they never recompute a prediction,
never fabricate a value and never turn "not detected" into "absent".
"""

import math
import os
from typing import Any, Callable, Dict, Iterable, List, Optional

from .constants import (
    MODULE_STATUS_FAILED,
    MODULE_STATUS_NOT_AVAILABLE,
    MODULE_STATUS_PARTIAL,
    MODULE_STATUS_SKIPPED,
    MODULE_STATUS_SUCCESS,
)
from .errors import ErrorCode, MethArCTError
from .capabilities import get_module

#: Maximum number of DIAMOND hits copied into the report per pathway.
MAX_EVIDENCE_HITS_PER_PATHWAY = 5

#: Maximum number of pathways copied into the report.
MAX_PATHWAYS_IN_REPORT = 200


def json_safe(value: Any, _depth: int = 0) -> Any:
    """Convert *value* into plain JSON-serialisable Python objects.

    NumPy scalars/arrays, pandas timestamps and path objects are handled by
    duck typing so this module stays importable without those packages.
    """
    if _depth > 12:  # pragma: no cover - defensive
        return str(value)
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return None
        return value
    if isinstance(value, dict):
        return {str(key): json_safe(item, _depth + 1) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [json_safe(item, _depth + 1) for item in value]

    to_dict = getattr(value, "to_dict", None)
    if callable(to_dict) and value.__class__.__name__ == "DataFrame":
        return json_safe(to_dict(orient="records"), _depth + 1)
    to_list = getattr(value, "tolist", None)
    if callable(to_list):
        return json_safe(to_list(), _depth + 1)
    item = getattr(value, "item", None)
    if callable(item):
        try:
            return json_safe(item(), _depth + 1)
        except (ValueError, TypeError):  # pragma: no cover - non-scalar
            pass
    isoformat = getattr(value, "isoformat", None)
    if callable(isoformat):
        return isoformat()
    return str(value)


def _number(value: Any) -> Optional[float]:
    """Return a finite float or ``None`` (never a silent zero)."""
    if value is None or isinstance(value, bool):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(result) or math.isinf(result):
        return None
    return result


def _integer(value: Any) -> Optional[int]:
    number = _number(value)
    return None if number is None else int(number)


def _text(value: Any) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _detection_status(pathway: Dict[str, Any]) -> str:
    """Map native pathway evidence onto the report vocabulary."""
    search_status = pathway.get("search_status")
    if search_status and search_status != "searched":
        return "not_assessed"
    hits = _integer(pathway.get("low_threshold_hits"))
    cultivability_hits = _integer(pathway.get("cultivability_hits"))
    if (hits or 0) > 0 or (cultivability_hits or 0) > 0:
        return "detected"
    completeness = _number(pathway.get("low_completeness_percentage"))
    if completeness is not None and completeness >= 100.0:
        return "detected"
    return "not_detected"


def _normalise_hit(hit: Dict[str, Any]) -> Dict[str, Any]:
    """Copy one DIAMOND hit into a stable schema fragment."""
    return {
        "query_id": _text(hit.get("qseqid")),
        "reference_id": _text(hit.get("sseqid")),
        "percent_identity": _number(hit.get("pident")),
        "alignment_length": _integer(hit.get("length")),
        "evalue": _number(hit.get("evalue")),
        "bitscore": _number(hit.get("bitscore")),
    }


def _pathway_entries(pathway_results: Dict[str, Any]) -> List[Dict[str, Any]]:
    entries = []
    for database, data in list(pathway_results.items())[:MAX_PATHWAYS_IN_REPORT]:
        if not isinstance(data, dict):
            continue
        hits = data.get("low_threshold_hits")
        entries.append(
            {
                "database": database,
                "pathway_name": data.get("pathway_name"),
                "detection_status": _detection_status(data),
                "search_status": data.get("search_status"),
                "reference_sequences": _integer(data.get("reference_sequences")),
                "high_threshold_hits": _integer(data.get("high_threshold_hits")),
                "low_threshold_hits": _integer(hits),
                "high_completeness_percentage": _number(data.get("high_completeness_percentage")),
                "low_completeness_percentage": _number(data.get("low_completeness_percentage")),
                "average_identity_high": _number(data.get("average_identity_high")),
                "average_identity_low": _number(data.get("average_identity_low")),
                "cultivability_hits": _integer(data.get("cultivability_hits")),
                "cultivability_percentage": _number(data.get("cultivability_percentage")),
                "cultivability_status": data.get("cultivability_status"),
            }
        )
    return entries


def _pathway_evidence(pathway_results: Dict[str, Any]) -> List[Dict[str, Any]]:
    evidence: List[Dict[str, Any]] = []
    for database, data in pathway_results.items():
        if not isinstance(data, dict):
            continue
        for key in ("best_hits_high", "best_hits_low", "best_hits"):
            hits = data.get(key)
            if not isinstance(hits, (list, tuple)):
                continue
            for hit in list(hits)[:MAX_EVIDENCE_HITS_PER_PATHWAY]:
                if not isinstance(hit, dict):
                    continue
                record = _normalise_hit(hit)
                record["database"] = database
                record["threshold"] = "high" if key == "best_hits_high" else "low"
                evidence.append(record)
    return evidence


# ---------------------------------------------------------------------------
# Per-module serializers
# ---------------------------------------------------------------------------

def _serialize_pathway(native: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
    summary = native.get("summary") or {}
    pathway_results = native.get("pathway_results") or {}
    entries = _pathway_entries(pathway_results)
    detected = [entry for entry in entries if entry["detection_status"] == "detected"]

    values = {
        "total_input_sequences": _integer(summary.get("total_input_sequences")),
        "databases_requested": _integer(summary.get("databases_requested")),
        "databases_searched": _integer(summary.get("databases_searched")),
        "raw_hit_pathways": _integer(summary.get("raw_hit_pathways")),
        "threshold_passed_pathways": _integer(summary.get("threshold_passed_pathways")),
        "pathways_detected": _integer(summary.get("pathways_detected")),
        "complete_pathways": _integer(summary.get("complete_pathways")),
        "total_hits": _integer(summary.get("total_hits")),
        "pathways": entries,
    }
    return {
        "result_type": "evidence",
        "summary": "DIAMOND homology search reported {0} pathway(s) with at least one "
        "threshold-passing hit out of {1} searched database(s).".format(
            len(detected), summary.get("databases_searched", 0)
        ),
        "values": values,
        "units": {"completeness_percentage": "percent"},
        "evidence": _pathway_evidence(pathway_results),
        "limitations": [
            "not_detected does not mean absent: absence of homology is not evidence "
            "that the pathway is biologically missing.",
            "Completeness percentages are computed against the bundled reference "
            "sequence counts, not against an experimentally validated gene set.",
        ],
    }


def _serialize_cultivation(native: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
    pathway_results = native.get("pathway_results") or {}
    assessment = native.get("cultivability_assessment") or {}
    amino_acids = native.get("amino_acid_biosynthesis") or {}
    aa_summary = amino_acids.get("summary") or {}

    pathways = []
    for name, data in list(pathway_results.items())[:MAX_PATHWAYS_IN_REPORT]:
        if not isinstance(data, dict):
            continue
        completeness = _number(data.get("completeness"))
        pathways.append(
            {
                "pathway_name": name,
                "completeness": completeness,
                "found_genes": _integer(data.get("found_genes")),
                "total_genes": _integer(data.get("total_genes")),
                "detection_status": "detected"
                if completeness is not None and completeness >= 1.0
                else "not_detected",
            }
        )

    values = {
        "total_pathways": _integer(assessment.get("total_pathways")) or len(pathway_results),
        "overall_cultivability": assessment.get("overall_cultivability"),
        "average_pathway_completeness": _number(assessment.get("average_pathway_completeness")),
        "amino_acid_completeness": _number(assessment.get("amino_acid_completeness")),
        "vitamin_coenzyme_completeness": _number(assessment.get("vitamin_coenzyme_completeness")),
        "nucleotide_completeness": _number(assessment.get("nucleotide_completeness")),
        "atp_completeness": _number(assessment.get("atp_completeness")),
        "amino_acid_summary": {
            "present": _integer(aa_summary.get("present")),
            "missing_or_partial": _integer(aa_summary.get("missing_or_partial")),
            "no_dedicated_pathway": _integer(aa_summary.get("no_dedicated_pathway")),
            "missing_list": json_safe(aa_summary.get("missing_list") or []),
        },
        "pathways": pathways,
    }
    return {
        "result_type": "assessment",
        "summary": "Cultivability assessed as '{0}' from biosynthetic pathway "
        "completeness (average {1}).".format(
            assessment.get("overall_cultivability"),
            assessment.get("average_pathway_completeness"),
        ),
        "values": values,
        "units": {"completeness": "fraction_of_reference_genes"},
        "evidence": [],
        "limitations": [
            "Cultivability categories follow the bundled rule thresholds and are not "
            "a measured growth phenotype.",
            "The rule-based confidence field is a heuristic rule level, not a "
            "calibrated probability.",
        ],
    }


def _serialize_temperature(native: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
    prediction = native.get("prediction") or {}
    summary = native.get("summary") or {}
    values = {
        "T_min": _number(prediction.get("T_min")),
        "T_opt": _number(prediction.get("T_opt")),
        "T_max": _number(prediction.get("T_max")),
        "temperature_category": summary.get("temperature_category"),
    }
    return {
        "result_type": "prediction",
        "summary": "Predicted growth temperature range {0}–{1} °C (optimum {2} °C, {3}).".format(
            values["T_min"], values["T_max"], values["T_opt"], values["temperature_category"]
        ),
        "values": values,
        "units": {"T_min": "degC", "T_opt": "degC", "T_max": "degC"},
        "evidence": [],
        "limitations": [
            "Predicted values are model outputs from sequence composition, not "
            "experimental growth measurements.",
        ],
    }


def _serialize_salinity(native: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
    prediction = native.get("prediction") or {}
    top3 = native.get("top3_predictions") or []
    values = {
        "salinity_label": prediction.get("salinity_label"),
        "salinity_index": _integer(prediction.get("salinity_index")),
        "confidence": _number(prediction.get("confidence")),
    }
    evidence = [
        {
            "rank": _integer(item.get("rank")),
            "label": item.get("label"),
            "probability": _number(item.get("probability")),
        }
        for item in top3
        if isinstance(item, dict)
    ]
    return {
        "result_type": "prediction",
        "summary": "Predicted salinity adaptation class '{0}'.".format(
            values["salinity_label"]
        ),
        "values": values,
        "units": {"confidence": "model_probability"},
        "evidence": evidence,
        "limitations": [
            "The salinity class is a model prediction, not a measured NaCl optimum.",
            "The model probability is not a calibrated experimental confidence value.",
        ],
    }


def _serialize_ph(native: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
    prediction = native.get("prediction") or {}

    def _value(field: str) -> Optional[float]:
        payload = prediction.get(field)
        if isinstance(payload, dict):
            return _number(payload.get("value"))
        return _number(payload)

    values = {
        "ph_optimum": _value("ph_optimum"),
        "ph_max": _value("ph_max"),
        "ph_min": _value("ph_min"),
        "is_novel": bool(native.get("is_novel", False)),
    }

    warnings: List[str] = []
    for payload in prediction.values():
        if isinstance(payload, dict):
            warning = _text(payload.get("warning"))
            if warning and warning not in warnings:
                warnings.append(warning)

    if any(values[field] is not None for field in ("ph_optimum", "ph_max", "ph_min")):
        summary = "Predicted growth pH optimum {0} (range {1}–{2}).".format(
            values["ph_optimum"], values["ph_min"], values["ph_max"]
        )
    else:
        summary = "The pH model returned no value for this input" + (
            " ({0}).".format("; ".join(warnings)) if warnings else "."
        )

    return {
        "result_type": "prediction",
        "summary": summary,
        "values": values,
        "units": {"ph_optimum": "pH", "ph_min": "pH", "ph_max": "pH"},
        "evidence": [{"kind": "warning", "message": warning} for warning in warnings],
        "limitations": [
            "Predicted pH values are model outputs, not experimental measurements.",
            "When is_novel is true the genome lies outside the training set and the "
            "prediction is less reliable.",
        ],
    }


def _serialize_antibiotic(native: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
    aai = native.get("aai_results") or {}
    rules = native.get("all_rules") or []
    values = {
        "recommended_antibiotics": json_safe(native.get("recommended_antibiotics") or []),
        "aai_percent": {str(key): _number(value) for key, value in aai.items()},
        "rules_matched": [
            rule.get("name")
            for rule in rules
            if isinstance(rule, dict) and rule.get("satisfied")
        ],
    }
    evidence = [
        {
            "rule": rule.get("name"),
            "satisfied": bool(rule.get("satisfied")),
            "details": json_safe(rule.get("details") or []),
        }
        for rule in rules
        if isinstance(rule, dict)
    ]
    return {
        "result_type": "prediction",
        "summary": "AAI-based antibiotic selection returned {0} recommendation(s).".format(
            len(values["recommended_antibiotics"])
        ),
        "values": values,
        "units": {"aai_percent": "percent_identity"},
        "evidence": evidence,
        "limitations": [
            "Recommendations follow fixed AAI thresholds against five reference "
            "genomes and are not a susceptibility test.",
        ],
    }


def _serialize_genome_quality(native: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
    summary = native.get("summary") or {}
    values = {
        "completeness": _number(summary.get("average_completeness")),
        "contamination": _number(summary.get("average_contamination")),
        "quality_score": _number(summary.get("average_quality_score")),
        "total_genomes": _integer(summary.get("total_genomes")),
        "high_quality_genomes": _integer(summary.get("high_quality_genomes")),
        "medium_quality_genomes": _integer(summary.get("medium_quality_genomes")),
        "low_quality_genomes": _integer(summary.get("low_quality_genomes")),
    }
    return {
        "result_type": "assessment",
        "summary": "CheckM2 reported completeness {0}% and contamination {1}%.".format(
            values["completeness"], values["contamination"]
        ),
        "values": values,
        "units": {"completeness": "percent", "contamination": "percent"},
        "evidence": [],
        "limitations": [
            "Genome quality estimates depend on the external CheckM2 reference database.",
        ],
    }


_SERIALIZERS = {
    "pathway": _serialize_pathway,
    "cultivation": _serialize_cultivation,
    "temperature": _serialize_temperature,
    "salinity": _serialize_salinity,
    "ph": _serialize_ph,
    "antibiotic": _serialize_antibiotic,
    "genome_quality": _serialize_genome_quality,
}

#: Native status strings produced by the core analyzers.
_STATUS_MAP = {
    "success": MODULE_STATUS_SUCCESS,
    "partial": MODULE_STATUS_PARTIAL,
    "failed": MODULE_STATUS_FAILED,
    "skipped": MODULE_STATUS_NOT_AVAILABLE,
    "not_available": MODULE_STATUS_NOT_AVAILABLE,
    "not_applicable": MODULE_STATUS_SKIPPED,
}


def _has_any_prediction(values: Dict[str, Any]) -> bool:
    """Return True when at least one scalar prediction value is present."""
    for value in values.values():
        if isinstance(value, bool):
            continue
        if isinstance(value, (int, float, str)):
            return True
    return False


def normalize_status(native_status: Optional[str]) -> str:
    """Map a native status string onto the report vocabulary."""
    if not native_status:
        return MODULE_STATUS_SUCCESS
    return _STATUS_MAP.get(str(native_status).lower(), MODULE_STATUS_PARTIAL)


def serialize_module(
    module_name: str,
    native_result: Any,
    runtime_seconds: Optional[float] = None,
    artifact_resolver: Optional[Callable[[str], Optional[str]]] = None,
    extra_limitations: Optional[Iterable[str]] = None,
) -> Dict[str, Any]:
    """Build the unified module-result envelope for one module.

    Args:
        module_name: Registered module identifier.
        native_result: Native analyzer result (dictionary) or ``None``.
        runtime_seconds: Wall-clock runtime of the module.
        artifact_resolver: Callback mapping an absolute output path to an
            artifact id declared in the job manifest.
        extra_limitations: Additional limitations appended to the registry ones.

    Returns:
        A JSON-serialisable module result following the report schema.
    """
    capability = get_module(module_name)
    limitations = list(capability.limitations)
    if extra_limitations:
        limitations.extend(extra_limitations)

    model = {"name": capability.model_name, "version": capability.model_version}

    if native_result is None:
        return {
            "module": module_name,
            "status": MODULE_STATUS_NOT_AVAILABLE,
            "result_type": "none",
            "summary": "The module produced no result.",
            "values": {},
            "units": {},
            "evidence": [],
            "limitations": limitations,
            "source_artifacts": [],
            "model": model,
            "error": None,
            "runtime_seconds": json_safe(runtime_seconds),
        }

    if not isinstance(native_result, dict):
        raise MethArCTError(
            ErrorCode.INTERNAL_ERROR,
            "Native module result is not a dictionary.",
            details={"module": module_name, "type": type(native_result).__name__},
        )

    status = normalize_status(native_result.get("status"))
    error_message = native_result.get("error") or native_result.get("message")
    serializer = _SERIALIZERS.get(module_name)

    payload: Dict[str, Any] = {
        "module": module_name,
        "status": status,
        "result_type": "none",
        "summary": None,
        "values": {},
        "units": {},
        "evidence": [],
        "limitations": limitations,
        "source_artifacts": [],
        "model": model,
        "error": _text(error_message) if status != MODULE_STATUS_SUCCESS else None,
        "runtime_seconds": json_safe(runtime_seconds),
    }

    if serializer is not None and status in (MODULE_STATUS_SUCCESS, MODULE_STATUS_PARTIAL):
        payload.update(serializer(native_result, {"module": capability, "status": status}))
        # A prediction module that produced no usable value is a partial result,
        # never a silent success.
        if (
            payload.get("status") == MODULE_STATUS_SUCCESS
            and payload.get("result_type") == "prediction"
            and not _has_any_prediction(payload.get("values") or {})
        ):
            payload["status"] = MODULE_STATUS_PARTIAL
    elif status == MODULE_STATUS_SUCCESS:
        payload["summary"] = "Module completed without a registered serializer."

    if status == MODULE_STATUS_NOT_AVAILABLE and not payload.get("summary"):
        payload["summary"] = _text(error_message) or "Module is not available in this environment."

    # Map native output files onto artifact ids declared in the manifest.
    source_artifacts = []
    output_files = native_result.get("output_files") or {}
    if isinstance(output_files, dict):
        for key, path in output_files.items():
            if not path:
                continue
            artifact_id = artifact_resolver(os.fspath(path)) if artifact_resolver else None
            if artifact_id:
                source_artifacts.append({"artifact_id": artifact_id, "kind": key})
    payload["source_artifacts"] = source_artifacts

    return json_safe(payload)
