# MethArCT AI API (version 1.0)

The AI layer is an **additive** adapter on top of the existing MethArCT
analyzers. It adds a stable Python API, a machine-readable CLI mode, a job
manager and an optional MCP server. It does not change any scientific
algorithm, threshold, database, model file or legacy output.

Three independent usage paths exist:

| Path | Entry point | Audience |
|---|---|---|
| Legacy | `metharct comprehensive ...`, `metharct.core.*` | existing users, unchanged |
| Python / CLI JSON | `metharct.ai.MethArCTClient`, `metharct capabilities --format json` | scripts, workflows |
| MCP | `metharct-mcp` (optional) | agent hosts (Gemma, LFM, Qwen, DeepSeek, ...) |

---

## 1. Installation

```bash
pip install .                # legacy MethArCT only (no AI extras required)
pip install ".[ai]"          # + jsonschema, for third-party schema validation
pip install ".[mcp]"         # + MCP Python SDK, for the MCP server (Python >= 3.10)
```

The AI layer itself only uses the Python standard library, so nothing changes
for users who only run the classic commands.

`requirements.txt` contains **runtime** dependencies only (it is also
`install_requires`). Development, test and documentation tooling lives in
`requirements-dev.txt`:

```bash
pip install -r requirements-dev.txt   # pytest, flake8, black, mypy, sphinx, jsonschema
```

### Python version notes

| Feature | Minimum Python |
|---|---|
| Legacy CLI / Python API / AI layer | 3.8 |
| MCP server (`metharct-mcp`) | **3.10** (MCP Python SDK constraint) |

On 3.8/3.9 the `mcp` extra is skipped by a `python_version` marker, `doctor`
reports the reason, and `metharct-mcp` exits with `DEPENDENCY_MISSING` whose
`details` names the required interpreter. The supported SDK range is
`mcp>=1.2` and both SDK generations are handled: 1.x (`FastMCP`) and 2.x
(`MCPServer`, where FastMCP was renamed).

---

## 2. Module registry

`metharct.ai.capabilities.MODULE_REGISTRY` is the single source of truth used
by the CLI, the MCP tools, the serializers and the report.

| Module id | What it does | Options |
|---|---|---|
| `pathway` | DIAMOND search against the methane / sulfur / nitrogen reference pathways | `databases` (array), `threads`, `evalue`, `max_target_seqs` |
| `cultivation` | Biosynthetic pathway completeness, amino-acid auxotrophy | `threads`, `evalue`, `identity_threshold` |
| `temperature` | MLHuaiNanzi T_min / T_opt / T_max | – |
| `salinity` | SuSha salinity adaptation class | – |
| `ph` | pH optimum / min / max | – |
| `antibiotic` | AAI-based antibiotic selection | `threads`, `evalue` |
| `genome_quality` | CheckM2 completeness / contamination | `input_type` (`fasta` \| `directory`) |

`enabled` is derived from real probes (executable present, `*.dmnd` present,
model files present, Python module importable). A module directory on disk is
never treated as proof of availability.

Unknown module ids and unknown option names raise `MODULE_UNKNOWN` /
`INVALID_ARGUMENT`; nothing is silently ignored.

---

## 3. Python API

```python
from metharct.ai import MethArCTClient

client = MethArCTClient(
    workspace="./metharct_workspace",   # job store; created if missing
    # allowed_roots=["D:/genomes"],     # optional read whitelist
    config_file=None,                   # optional metharct_config.yaml
    max_parallel_jobs=1,
    max_file_size_mb=500,
    default_timeout_seconds=None,
)
```

### 3.1 Discovery

| Method | Returns |
|---|---|
| `get_capabilities(include_details=True)` | version, environment, modules, options, supports |
| `doctor(check_databases=True)` | environment checks with `available` / `unavailable` / `misconfigured` / `not_required` |
| `inspect_input(path, expected_type="auto", max_file_size_mb=None)` | FASTA validation + compatible modules |
| `get_version_info()` | package version, AI API version, schema version, workspace |

`doctor` only reports. It never installs packages, downloads databases or
modifies `PATH`.

### 3.2 Running a job

```python
job = client.run_analysis(
    input_file="MZD129.faa",
    input_type="auto",                  # auto | protein_fasta | nucleotide_fasta
    modules=["temperature", "ph"],      # None = every registered module
    output_dir=None,                    # recorded in metadata; results always go to the job dir
    options={"pathway": {"threads": 8}},
    timeout_seconds=None,
    metadata={"strain": "MZD129"},
    idempotency_key=None,
    wait=True,
    poll_interval=1.0,
)
print(job["job_id"], job["status"])
```

Returns a job view without host absolute paths:

```json
{
  "job_id": "job_20261117T074626_d3f1ebc6",
  "status": "succeeded",
  "phase": "succeeded",
  "requested_modules": ["temperature", "ph"],
  "completed_modules": ["temperature", "ph"],
  "failed_modules": [],
  "skipped_modules": [],
  "warnings": [],
  "report_path": "metharct_ai_report.json",
  "status_tool": "metharct_get_job_status",
  "result_tool": "metharct_get_result_summary",
  "workspace_job_dir": "jobs/job_20261117T074626_d3f1ebc6",
  "unavailable_modules": [],
  "deduplicated": false
}
```

Other job methods:

| Method | Purpose |
|---|---|
| `submit_analysis(input_file, **kwargs)` | same as `run_analysis(wait=False)` |
| `get_job(job_id)` | job view |
| `get_job_status(job_id)` | compact status (phase, progress, errors) |
| `list_jobs(limit=None)` | newest first |
| `cancel_job(job_id)` | best-effort termination, honest status |
| `get_result(job_id, validate=False)` | full unified report |
| `get_result_summary(job_id, modules=None)` | compact report for an LLM context |
| `get_module_result(job_id, module)` | one module envelope |
| `list_artifacts(job_id)` | manifest entries |
| `read_artifact(job_id, artifact_id, offset=0, limit=50, filter_column=None, filter_value=None, contains=False)` | bounded page |

### 3.3 Job state machine

```text
queued -> validating -> running -> succeeded
                              \-> partially_succeeded
                              \-> failed
                              \-> cancelled
```

Transitions are validated; illegal transitions raise `JOB_ALREADY_FINISHED`.
Timestamps (`created_at`, `started_at`, `finished_at`) and
`duration_seconds` are written by the job store, never inferred from logs.

Each job owns:

```text
<workspace>/jobs/<job_id>/{job.json, metharct_ai_report.json, manifest.json,
                           logs/stdout.log, logs/stderr.log, artifacts/}
```

### 3.4 Idempotency

`idempotency_key` reuses an existing job unless that job failed or was
cancelled. The returned view then has `"deduplicated": true`. This prevents an
agent from re-running an expensive DIAMOND search after a transport timeout.

### 3.5 Parallelism

`max_parallel_jobs` (default `1`) caps concurrent jobs per client. Exceeding it
raises a retryable `INVALID_ARGUMENT` with
`details.reason == "parallel_job_limit"`.

---

## 4. CLI (machine-readable mode)

```bash
metharct capabilities --format json
metharct doctor --format json
metharct inspect MZD129.faa --format json
metharct run MZD129.faa --modules temperature,ph --format json
metharct result <job_id> --format json [--module temperature]
metharct artifact <job_id> [<artifact_id>] --format json [--offset 0 --limit 50]
metharct cancel <job_id> --format json
metharct jobs --format json
```

Rules:

* `--format json` prints **exactly one JSON object** on stdout, UTF-8 encoded
  (written to `sys.stdout.buffer`, so a non-UTF-8 console code page cannot
  corrupt it). All progress and diagnostics go to stderr or the job logs.
* `--json-out <path>` writes the same payload atomically (temp file +
  `os.replace`).
* Default (`--format text`) output stays human readable.
* `--workspace <dir>` selects the job store (default `metharct_workspace`,
  override with `METHARCT_WORKSPACE`).

Exit codes: `0` success, `1` error (JSON still contains the error), `2` usage,
`3` job failed / partially succeeded.

The legacy subcommands and their arguments are unchanged; they do not gain a
`--format` flag.

---

## 5. Input inspection

`inspect` is streaming and never loads the whole file:

* file exists, is a regular readable file, is not compressed or binary;
* size limit (`max_file_size_mb`, `<= 0` disables);
* FASTA header structure, at least one record, duplicate ids, empty headers;
* invalid (non-IUPAC) residue characters;
* sequence length min / median / max, total residues;
* sequence type from the actual alphabet (`protein` / `nucleotide`), with a
  rule-level confidence (`high` / `medium` / `low`) — never a calibrated
  probability;
* SHA-256, size and modification time, so a later run can prove it analysed the
  same bytes (the worker re-hashes and refuses a changed input).

The file extension is **not** used to decide the sequence type.

---

## 6. Error codes

| Code | Meaning |
|---|---|
| `INPUT_FILE_NOT_FOUND` | input path does not exist |
| `INPUT_FILE_UNREADABLE` | not readable, or outside the allowed roots |
| `INPUT_FILE_TOO_LARGE` | exceeds the configured size limit |
| `INPUT_FORMAT_INVALID` | not FASTA / empty / wrong header structure |
| `INPUT_TYPE_AMBIGUOUS` | sequence type could not be established |
| `INPUT_MODULE_INCOMPATIBLE` | input type not supported by the module |
| `MODULE_UNKNOWN` | unknown module id |
| `MODULE_UNAVAILABLE` | no requested module can run here |
| `DEPENDENCY_MISSING` | optional dependency missing (e.g. MCP SDK) |
| `DATABASE_MISSING` / `DATABASE_INVALID` | database missing / unusable |
| `ANALYSIS_TIMEOUT` | job exceeded its timeout |
| `ANALYSIS_FAILED` | the analyzer raised |
| `JOB_NOT_FOUND` | unknown job id in this workspace |
| `JOB_ALREADY_FINISHED` | illegal state transition / cancel after finish |
| `ARTIFACT_NOT_FOUND` | unknown artifact id or file missing |
| `ARTIFACT_ACCESS_DENIED` | path escapes the job dir / media type not readable |
| `SCHEMA_VALIDATION_FAILED` | a payload does not satisfy its JSON schema |
| `INVALID_ARGUMENT` | invalid argument or option |
| `INTERNAL_ERROR` | unclassified failure |

Codes are stable. Extra context lives in `details`; the original exception type
and message are preserved there.

---

## 7. Interpretation contract

The API reports what the models produced and nothing more:

* A missing value is `null`, never `0` and never "not detected".
* A prediction is not a measurement; every prediction module carries that
  limitation string.
* Pathway `not_detected` is not biological absence.
* No calibrated confidence is invented. Where the model emits a probability it
  is labelled `model_probability`; the core cultivability confidence is
  described as a rule heuristic.
* A module that produced no usable value is reported `partial`, and a job with
  no usable result is `failed`, never `succeeded`.

---

## 8. Compatibility

* Legacy CLI text output, arguments, exit behaviour, Python classes, method
  signatures and output filenames are unchanged (frozen by
  `tests/regression/test_legacy_compat.py`).
* `metharct.core` is untouched.
* Only additive files were created; `metharct/cli/__init__.py` became lazy in
  order to remove a `runpy` double-import warning (`from metharct.cli import
  main` still works).
* `schema_version` of the report is versioned independently of the package.
