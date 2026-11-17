# Changelog

All notable changes to MethArCT are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and the project adheres to [Semantic Versioning](https://semver.org/).

## [0.8.0] - 2026-11-17

### Added

**AI-ready adaptation layer (`metharct/ai/`) — additive, optional, no science
changes.**

* `MethArCTClient` — stable synchronous Python API:
  `get_capabilities`, `doctor`, `inspect_input`, `run_analysis`,
  `submit_analysis`, `get_job`, `get_job_status`, `list_jobs`, `cancel_job`,
  `get_result`, `get_result_summary`, `get_module_result`, `list_artifacts`,
  `read_artifact`.
* Single module registry (`MODULE_REGISTRY`) with dynamic availability probes
  for `pathway`, `cultivation`, `temperature`, `salinity`, `ph`, `antibiotic`,
  `genome_quality`.
* Streaming FASTA inspection: encoding, header structure, duplicate ids,
  invalid characters, length statistics, protein/nucleotide classification with
  a rule-level confidence, SHA-256.
* Job manager with an explicit state machine
  (`queued → validating → running → succeeded/partially_succeeded/failed/cancelled`),
  atomic `job.json`, per-job `logs/`, `artifacts/` and `manifest.json`.
* Isolated worker process (`python -m metharct.ai.worker`, `shell=False`),
  best-effort cancellation, idempotency keys, parallel-job limit.
* Unified report `metharct_ai_report.json` (`schema_version` **1.0**) with
  per-module status, values + units, evidence, limitations, artifacts (SHA-256)
  and provenance.
* JSON schemas for the report, job record, capabilities and input inspection
  (`metharct/schemas/`), plus a dependency-free validator.
* Machine-readable CLI: `capabilities`, `doctor`, `inspect`, `run`, `result`,
  `artifact`, `cancel`, `jobs` with `--format json` / `--json-out`.
* Optional MCP server (`metharct-mcp`) with eight tools and four resources;
  the MCP SDK is only required when the server is used. Both SDK generations
  are supported (`mcp>=1.2`: 1.x `FastMCP` and 2.x `MCPServer`).
* `metharct/ai/clock.py` — timezone-aware UTC timestamps (no deprecated
  `datetime.utcnow()`; clean on Python 3.12).
* Curated knowledge base (`knowledge/`) and documentation
  (`docs/AI_ADAPTATION_AUDIT.md`, `AI_API.md`, `MCP_TOOLS.md`,
  `AI_REPORT_SCHEMA.md`, `SECURITY.md`, `docs/AI_ADAPTATION_REPORT.md`).
* Test suite: unit, integration, end-to-end and legacy-compatibility tests
  (`tests/ai/`, `tests/regression/`, `tests/fixtures/`), including a real MCP
  `ClientSession` over stdio (`tests/ai/test_mcp_integration.py`, Python 3.10+)
  and packaging guards for the requirements split.

### Changed

* `metharct/cli/__init__.py` now resolves its exports lazily. This removes the
  `runpy` double-import warning from `python -m metharct.cli.main` and makes
  `import metharct.cli` cheap. `from metharct.cli import main` still works.
* `setup.py` gained the `ai` and `mcp` extra requirements (the MCP extra is
  guarded by `python_version >= "3.10"`, the SDK's own constraint), the
  `metharct-mcp` console script, and ships `metharct/schemas/*.json`.
* **`requirements.txt` now lists runtime dependencies only.** Development, test
  and documentation tooling (`pytest`, `black`, `flake8`, `mypy`, `jsonschema`,
  Sphinx) moved to the new `requirements-dev.txt`. `setup.py` still reads
  `requirements.txt` as `install_requires`, so end users no longer receive
  linters or Sphinx. `tests/regression/test_packaging.py` enforces this.
* `metharct doctor` reports the optional MCP SDK (and explains the Python 3.10+
  requirement when the interpreter is older).
* Package version bumped `0.7.3 → 0.8.0`.

### Unchanged (compatibility statement)

* No scientific algorithm, threshold, database, reference sequence or model file
  was modified. `metharct/core/` is byte-for-byte untouched.
* All legacy CLI subcommands (`comprehensive`, `diamond`, `checkm2`,
  `cultivation`, `susha`, `mlhuainanzi`, `ph`, `antibiotic`, `batch`), their
  arguments, defaults, text output, exit behaviour and output filenames are
  unchanged.
* All legacy Python classes and method signatures are unchanged
  (`tests/regression/test_legacy_compat.py`).
* The classic commands do **not** gain `--format`; only the new subcommands do.

### Known issues

* Cancellation stops a worker between modules; it does not interrupt a running
  DIAMOND process.
* Only Windows was exercised for this release; the CI test path for Linux has
  not been run here.
* The MCP SDK 1.x code path (`mcp.server.fastmcp.FastMCP`) is an import fallback;
  the end-to-end MCP session was verified with SDK 2.3.0.
* The pH model can return `null` when a genome lacks the required features; the
  module is then reported `partial` (never a fabricated value).

## [0.7.3] - previous release

See the git history before `bf77997995796bc346a09c42da3f7a1c44481c9d`.
