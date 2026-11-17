# -*- coding: utf-8 -*-
"""Tests for the machine-readable CLI mode.

Every JSON invocation is executed in a subprocess so stdout purity (no banner,
no progress text, no BOM) is verified exactly as an agent would see it.
"""

import json
import os
import subprocess
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CLI = [sys.executable, "-m", "metharct.cli.main"]


def _run(args, cwd=REPO_ROOT):
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONPATH"] = REPO_ROOT + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        CLI + args,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
        timeout=300,
    )


def _json(result):
    """Decode stdout as UTF-8 JSON and assert it is a single object."""
    text = result.stdout.decode("utf-8")
    assert not text.startswith("\ufeff"), "stdout must not start with a BOM"
    assert text.lstrip().startswith("{"), text[:80]
    return json.loads(text)


def test_root_help_lists_legacy_and_ai_commands():
    result = _run(["--help"])
    assert result.returncode == 0
    output = result.stdout.decode("utf-8", "replace")
    for command in (
        "comprehensive",
        "diamond",
        "checkm2",
        "cultivation",
        "susha",
        "mlhuainanzi",
        "ph",
        "antibiotic",
        "batch",
        "capabilities",
        "doctor",
        "inspect",
        "run",
        "result",
        "artifact",
        "cancel",
        "jobs",
    ):
        assert command in output, command


def test_capabilities_json_is_pure_json():
    result = _run(["capabilities", "--format", "json", "--no-details"])
    payload = _json(result)
    assert result.returncode == 0
    assert payload["ok"] is True
    assert payload["error"] is None
    assert set(payload) == {"ok", "request_id", "data", "warnings", "error"}
    assert payload["data"]["software"]["name"] == "MethArCT"
    assert payload["data"]["modules"]


def test_doctor_json_is_pure_json():
    result = _run(["doctor", "--format", "json", "--skip-databases"])
    payload = _json(result)
    assert result.returncode == 0
    assert payload["ok"] is True
    assert "checks" in payload["data"]
    assert payload["data"]["summary"]["available"] > 0


def test_inspect_json_reports_protein(protein_fasta):
    result = _run(["inspect", protein_fasta, "--format", "json"])
    payload = _json(result)
    assert result.returncode == 0
    assert payload["data"]["sequence_type"] == "protein"
    assert payload["data"]["sequence_count"] == 3
    assert payload["data"]["valid"] is True


def test_inspect_missing_file_returns_stable_error(tmp_path):
    result = _run(["inspect", str(tmp_path / "nope.faa"), "--format", "json"])
    payload = _json(result)
    assert result.returncode == 1
    assert payload["ok"] is False
    assert payload["error"]["code"] == "INPUT_FILE_NOT_FOUND"
    assert "message" in payload["error"]


def test_json_out_writes_file_atomically(tmp_path, protein_fasta):
    target = tmp_path / "inspection.json"
    result = _run(["inspect", protein_fasta, "--format", "json", "--json-out", str(target)])
    assert result.returncode == 0
    assert os.path.isfile(str(target))
    with open(str(target), "r", encoding="utf-8") as stream:
        payload = json.load(stream)
    assert payload["data"]["sequence_type"] == "protein"
    leftovers = [name for name in os.listdir(str(tmp_path)) if name.startswith(".tmp_")]
    assert leftovers == []


def test_result_unknown_job_returns_job_not_found(tmp_path):
    result = _run(
        [
            "result",
            "job_20260101T000000_deadbeef",
            "--format",
            "json",
            "--workspace",
            str(tmp_path / "ws"),
        ]
    )
    payload = _json(result)
    assert result.returncode == 1
    assert payload["error"]["code"] == "JOB_NOT_FOUND"


def test_artifact_unknown_job_returns_error(tmp_path):
    result = _run(
        [
            "artifact",
            "job_20260101T000000_deadbeef",
            "--format",
            "json",
            "--workspace",
            str(tmp_path / "ws"),
        ]
    )
    payload = _json(result)
    assert result.returncode == 1
    assert payload["error"]["code"] == "JOB_NOT_FOUND"


def test_cancel_unknown_job_returns_error(tmp_path):
    result = _run(
        [
            "cancel",
            "job_20260101T000000_deadbeef",
            "--format",
            "json",
            "--workspace",
            str(tmp_path / "ws"),
        ]
    )
    payload = _json(result)
    assert result.returncode == 1
    assert payload["error"]["code"] == "JOB_NOT_FOUND"


def test_text_mode_capabilities_still_human_readable():
    result = _run(["capabilities"])
    assert result.returncode == 0
    output = result.stdout.decode("utf-8", "replace")
    assert "MethArCT" in output
    assert "Modules:" in output


@pytest.mark.parametrize(
    "command",
    ["comprehensive", "diamond", "checkm2", "cultivation", "susha", "mlhuainanzi", "ph", "antibiotic", "batch"],
)
def test_legacy_subcommand_help_is_unchanged(command):
    """Every legacy subcommand must still parse its own --help."""
    result = _run([command, "--help"])
    assert result.returncode == 0, result.stderr
    assert "usage:" in result.stdout.decode("utf-8", "replace")
