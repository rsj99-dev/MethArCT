# -*- coding: utf-8 -*-
"""Tests for the job lifecycle store."""

import json
import os

import pytest

from metharct.ai.constants import (
    JOB_STATUS_FAILED,
    JOB_STATUS_QUEUED,
    JOB_STATUS_RUNNING,
    JOB_STATUS_SUCCEEDED,
)
from metharct.ai.errors import ErrorCode, MethArCTError
from metharct.ai.jobs import JobManager, summarize_status


@pytest.fixture()
def manager(tmp_path):
    return JobManager(str(tmp_path / "workspace"))


def _create(manager, **overrides):
    payload = {
        "input_file": "input.faa",
        "input_type": "protein_fasta",
        "modules": ["temperature"],
        "input_sha256": "a" * 64,
    }
    payload.update(overrides)
    return manager.create_job(**payload)


def test_create_job_writes_layout_and_record(manager):
    job = _create(manager)
    job_dir = manager.job_dir(job["job_id"])
    assert os.path.isdir(str(job_dir))
    assert os.path.isdir(os.path.join(str(job_dir), "artifacts"))
    assert os.path.isdir(os.path.join(str(job_dir), "logs"))

    record = os.path.join(str(job_dir), "job.json")
    assert os.path.isfile(record)
    with open(record, "r", encoding="utf-8") as stream:
        assert json.load(stream)["job_id"] == job["job_id"]


def test_job_id_is_unique_and_well_formed(manager):
    first = _create(manager)["job_id"]
    second = _create(manager)["job_id"]
    assert first != second
    assert first.startswith("job_")


def test_transition_lifecycle(manager):
    job = _create(manager)
    job = manager.transition(job["job_id"], JOB_STATUS_RUNNING)
    assert job["status"] == JOB_STATUS_RUNNING
    assert job["started_at"] is not None

    job = manager.transition(job["job_id"], JOB_STATUS_SUCCEEDED)
    assert job["status"] == JOB_STATUS_SUCCEEDED
    assert job["finished_at"] is not None
    assert job["duration_seconds"] is not None


def test_queued_job_cannot_jump_to_succeeded(manager):
    job = _create(manager)
    with pytest.raises(MethArCTError) as excinfo:
        manager.transition(job["job_id"], JOB_STATUS_SUCCEEDED)
    assert excinfo.value.code == ErrorCode.JOB_ALREADY_FINISHED
    assert manager.load(job["job_id"])["status"] == JOB_STATUS_QUEUED


def test_illegal_transition_raises(manager):
    job = _create(manager)
    manager.transition(job["job_id"], JOB_STATUS_RUNNING)
    manager.transition(job["job_id"], JOB_STATUS_FAILED)
    with pytest.raises(MethArCTError) as excinfo:
        manager.transition(job["job_id"], JOB_STATUS_RUNNING)
    assert excinfo.value.code == ErrorCode.JOB_ALREADY_FINISHED


def test_load_unknown_job_raises(manager):
    with pytest.raises(MethArCTError) as excinfo:
        manager.load("job_does_not_exist")
    assert excinfo.value.code == ErrorCode.JOB_NOT_FOUND


def test_load_rejects_path_traversal(manager):
    with pytest.raises(MethArCTError) as excinfo:
        manager.load("../escape")
    assert excinfo.value.code == ErrorCode.JOB_NOT_FOUND


def test_list_and_active_jobs(manager):
    first = _create(manager)
    second = _create(manager)
    manager.transition(second["job_id"], JOB_STATUS_RUNNING)
    manager.transition(second["job_id"], JOB_STATUS_SUCCEEDED)

    assert {job["job_id"] for job in manager.list_jobs()} == {first["job_id"], second["job_id"]}
    assert [job["job_id"] for job in manager.active_jobs()] == [first["job_id"]]

    limited = manager.list_jobs(limit=1)
    assert len(limited) == 1
    assert limited[0]["job_id"] in {first["job_id"], second["job_id"]}


def test_list_jobs_is_ordered_newest_first(manager):
    older = _create(manager)
    newer = _create(manager)
    manager.update(older["job_id"], created_at="2020-01-01T00:00:00Z")
    manager.update(newer["job_id"], created_at="2030-01-01T00:00:00Z")
    assert [job["job_id"] for job in manager.list_jobs()] == [newer["job_id"], older["job_id"]]


def test_idempotency_key_lookup(manager):
    job = _create(manager, idempotency_key="abc")
    manager.transition(job["job_id"], JOB_STATUS_RUNNING)
    manager.transition(job["job_id"], JOB_STATUS_SUCCEEDED)

    found = manager.find_by_idempotency_key("abc")
    assert found is not None and found["job_id"] == job["job_id"]
    assert manager.find_by_idempotency_key("other") is None


def test_failed_jobs_are_not_reused_for_idempotency(manager):
    job = _create(manager, idempotency_key="abc")
    manager.transition(job["job_id"], JOB_STATUS_FAILED)
    assert manager.find_by_idempotency_key("abc") is None


def test_public_view_strips_absolute_paths(manager):
    job = _create(manager, input_file=str(manager.workspace / "genomes" / "in.faa"))
    view = JobManager.public_view(job)
    assert "input_file" not in view
    assert "workspace" not in view
    assert view["input_filename"] == "in.faa"
    assert ":" not in view["report_path"]


def test_cancel_request_on_finished_job_raises(manager):
    job = _create(manager)
    manager.transition(job["job_id"], JOB_STATUS_FAILED)
    with pytest.raises(MethArCTError) as excinfo:
        manager.request_cancel(job["job_id"])
    assert excinfo.value.code == ErrorCode.JOB_ALREADY_FINISHED


def test_summarize_status_has_no_host_paths(manager):
    job = _create(manager)
    summary = summarize_status(job)
    assert summary["job_id"] == job["job_id"]
    assert summary["status"] == JOB_STATUS_QUEUED
    assert "input_file" not in summary


def test_errors_and_warnings_are_appended(manager):
    job = _create(manager)
    manager.append_warning(job["job_id"], "careful")
    updated = manager.append_error(job["job_id"], "ANALYSIS_FAILED", "boom")
    assert updated["warnings"] == ["careful"]
    assert updated["errors"][0]["message"] == "boom"
    assert updated["error_code"] == "ANALYSIS_FAILED"


def test_job_record_matches_schema(manager):
    from metharct.ai.report import JOB_SCHEMA_FILENAME, validate_payload

    job = _create(manager)
    manager.transition(job["job_id"], JOB_STATUS_RUNNING)
    record = manager.load(job["job_id"])
    problems = validate_payload(record, JOB_SCHEMA_FILENAME)
    assert problems == [], problems[:5]
