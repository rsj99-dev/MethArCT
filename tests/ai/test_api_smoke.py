# -*- coding: utf-8 -*-
"""End-to-end tests: submit a job, run the worker, read the report.

These tests need at least one analysis module to be available in the current
environment; they are skipped (never faked) when nothing can run here.
"""

import json
import os

import pytest

from metharct.ai.capabilities import get_capabilities
from metharct.ai.errors import ErrorCode, MethArCTError
from metharct.ai.report import REPORT_SCHEMA_FILENAME, validate_payload


@pytest.fixture(scope="module")
def available_modules():
    capabilities = get_capabilities()
    return {module["name"] for module in capabilities["modules"] if module["enabled"]}


@pytest.fixture(scope="module")
def runnable_modules(available_modules):
    """Modules that run without DIAMOND and finish quickly."""
    candidates = [name for name in ("temperature", "salinity", "ph") if name in available_modules]
    if not candidates:
        pytest.skip("No embedded-model module is available in this environment.")
    return candidates


def test_run_analysis_end_to_end(client, protein_fasta, runnable_modules):
    module = runnable_modules[0]
    job = client.run_analysis(
        input_file=protein_fasta,
        modules=[module],
        poll_interval=0.2,
        wait=True,
    )
    assert job["status"] in ("succeeded", "partially_succeeded"), job.get("errors")
    assert job["job_id"].startswith("job_")
    assert "input_file" not in job

    status = client.get_job_status(job["job_id"])
    assert status["status"] == job["status"]
    assert status["started_at"] and status["finished_at"]

    report = client.get_result(job["job_id"], validate=True)
    problems = validate_payload(report, REPORT_SCHEMA_FILENAME)
    assert problems == [], problems[:5]
    assert report["run"]["job_id"] == job["job_id"]
    assert report["input"]["sequence_count"] == 3
    assert len(report["input"]["sha256"]) == 64

    module_result = client.get_module_result(job["job_id"], module)
    assert module_result["module"] == module
    assert module_result["model"]["name"]

    summary = client.get_result_summary(job["job_id"], modules=[module])
    assert summary["report_available"] is True
    assert module in summary["modules"]


def test_artifact_pages_do_not_leak_host_paths(client, protein_fasta, runnable_modules):
    """A job artifact must not expose the host input path to an agent."""
    job = client.run_analysis(protein_fasta, modules=runnable_modules, poll_interval=0.2)
    artifacts = client.list_artifacts(job["job_id"])
    input_dir = os.path.dirname(protein_fasta)
    for entry in artifacts:
        if entry["media_type"] not in (
            "text/csv",
            "text/tab-separated-values",
            "application/json",
            "text/plain",
        ):
            continue
        page = client.read_artifact(job["job_id"], entry["artifact_id"], offset=0, limit=50)
        payload = json.dumps(page, ensure_ascii=False)
        assert input_dir not in payload, entry["path"]
        assert ":/" not in payload.replace("metharct://", ""), entry["path"]


def test_worker_logs_are_captured(client, protein_fasta, runnable_modules):
    job = client.run_analysis(protein_fasta, modules=[runnable_modules[0]], poll_interval=0.2)
    job_dir = str(client.manager.job_dir(job["job_id"]))
    stdout_log = os.path.join(job_dir, "logs", "stdout.log")
    stderr_log = os.path.join(job_dir, "logs", "stderr.log")
    assert os.path.isfile(stdout_log)
    assert os.path.isfile(stderr_log)
    with open(stdout_log, "r", encoding="utf-8", errors="replace") as stream:
        content = stream.read()
    assert "metharct-worker" in content
    assert job["status"] in ("succeeded", "partially_succeeded")


def test_report_and_job_record_are_on_disk(client, protein_fasta, runnable_modules):
    job = client.run_analysis(protein_fasta, modules=[runnable_modules[0]], poll_interval=0.2)
    job_dir = str(client.manager.job_dir(job["job_id"]))
    assert os.path.isfile(os.path.join(job_dir, "job.json"))
    assert os.path.isfile(os.path.join(job_dir, "metharct_ai_report.json"))
    assert os.path.isfile(os.path.join(job_dir, "manifest.json"))
    with open(os.path.join(job_dir, "job.json"), "r", encoding="utf-8") as stream:
        record = json.load(stream)
    assert record["input_sha256"]
    assert record["report_path"] == "metharct_ai_report.json"


def test_artifacts_are_registered_and_paged(client, protein_fasta, runnable_modules):
    job = client.run_analysis(protein_fasta, modules=[runnable_modules[0]], poll_interval=0.2)
    artifacts = client.list_artifacts(job["job_id"])
    assert artifacts, "expected at least one artifact"
    for entry in artifacts:
        assert entry["path"].startswith("artifacts/")
        assert ".." not in entry["path"]
        assert ":" not in entry["path"], "artifact paths must stay relative"

    tabular = next(
        (entry for entry in artifacts if entry["media_type"] in ("text/tab-separated-values", "text/csv")),
        None,
    )
    if tabular is None:
        pytest.skip("The available module produced no tabular artifact.")
    page = client.read_artifact(job["job_id"], tabular["artifact_id"], offset=0, limit=1)
    assert page["returned_rows"] <= 1
    assert page["columns"]


def test_idempotency_key_reuses_the_job(client, protein_fasta, runnable_modules):
    module = runnable_modules[0]
    first = client.run_analysis(
        protein_fasta, modules=[module], idempotency_key="same-key", poll_interval=0.2
    )
    second = client.run_analysis(
        protein_fasta, modules=[module], idempotency_key="same-key", poll_interval=0.2
    )
    assert second["deduplicated"] is True
    assert second["job_id"] == first["job_id"]


def test_nucleotide_input_is_rejected_for_protein_modules(client, nucleotide_fasta):
    with pytest.raises(MethArCTError) as excinfo:
        client.run_analysis(nucleotide_fasta, modules=["temperature"])
    assert excinfo.value.code == ErrorCode.INPUT_MODULE_INCOMPATIBLE


def test_unknown_module_is_rejected(client, protein_fasta):
    with pytest.raises(MethArCTError) as excinfo:
        client.run_analysis(protein_fasta, modules=["metabolomics"])
    assert excinfo.value.code == ErrorCode.MODULE_UNKNOWN


def test_all_unavailable_modules_are_rejected(client, protein_fasta):
    capabilities = get_capabilities()
    disabled = [module["name"] for module in capabilities["modules"] if not module["enabled"]]
    if not disabled:
        pytest.skip("Every module is available here; nothing to assert.")
    if len(disabled) != len(capabilities["modules"]):
        pytest.skip("At least one module is available; the all-unavailable path is unreachable.")
    with pytest.raises(MethArCTError) as excinfo:
        client.run_analysis(protein_fasta, modules=disabled)
    assert excinfo.value.code == ErrorCode.MODULE_UNAVAILABLE


def test_parallel_job_limit_is_enforced(client, protein_fasta, runnable_modules):
    """An already-active job must block a second submission (deterministic)."""
    active = client.manager.create_job(
        input_file=protein_fasta,
        input_type="protein_fasta",
        modules=[runnable_modules[0]],
    )
    assert active["status"] == "queued"

    with pytest.raises(MethArCTError) as excinfo:
        client.run_analysis(protein_fasta, modules=[runnable_modules[0]])
    assert excinfo.value.retryable is True
    assert excinfo.value.details["reason"] == "parallel_job_limit"


def test_cancel_job_returns_an_honest_contract(client, protein_fasta, runnable_modules):
    job = client.submit_analysis(protein_fasta, modules=[runnable_modules[0]])
    result = client.cancel_job(job["job_id"])
    assert result["cancellation_supported"] is True
    assert isinstance(result["accepted"], bool)
    final = client.get_job_status(job["job_id"])
    # Whatever happened, the reported job status must be a real terminal or
    # in-flight state -- never a fabricated "cancelled" without a status change.
    assert final["status"] in (
        "cancelled",
        "succeeded",
        "partially_succeeded",
        "failed",
        "running",
        "validating",
        "queued",
    )
    if result["accepted"]:
        assert result["status"] == final["status"]
