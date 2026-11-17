# MethArCT AI-ready upgrade — delivery report

Companion to `AI_ADAPTATION_AUDIT.md` (the pre-upgrade audit) and
`CHANGELOG.md` (user-facing notes).

## 1. Verification (executed)

| Command | Result |
|---|---|
| `python -m pytest tests -q` — Python 3.8.20 | **172 passed, 6 skipped**, 63 s |
| `python -m pytest tests -q` — Python 3.12.7 + `mcp 2.3.0` + `jsonschema 4.26.0` | **177 passed, 1 skipped**, 0 warnings, 50 s |
| `python -m pytest tests/ai/test_mcp_integration.py -q` (3.12 + MCP SDK) | **4 passed** — real `ClientSession` over stdio |
| `python -m flake8 metharct/ai tests/ai tests/regression --select=F,E9` | clean (exit 0) |
| `python -m metharct.cli.main --help` | exit 0, all legacy + 8 AI subcommands listed |
| `python -m metharct.cli.main <legacy> --help` (9 subcommands) | exit 0 each |
| `python -m metharct.cli.main capabilities --format json` | exit 0, single JSON object, UTF-8, no BOM |
| `python -m metharct.cli.main doctor --format json` | exit 0 |
| `python -m metharct.cli.main inspect <fixture> --format json` | exit 0, protein/high, 3 sequences |
| `python -m metharct.cli.main run <fixture> -m temperature --format json` | exit 0, `succeeded` |
| `python -m metharct.cli.main result <job_id> --format json` | exit 0, `T_opt = 80.14 degC` |
| `python -m metharct.cli.main artifact <job_id> a001 --limit 2 --format json` | exit 0, 1 row, 5 columns |
| `python -m metharct.cli.main jobs --format json` | exit 0 |
| `python -m metharct.ai.mcp_server --list-tools` | 8 tools with JSON schemas |
| `metharct-mcp` under MCP SDK 2.3.0 (stdio, `ClientSession`) | initialize, `tools/list` (8), `resources/list` (4), `resources/read` ×2, `run_analysis` → `succeeded`, `get_job_status`, `get_result_summary`, `read_artifact`, `cancel_job`, `inspect_input` outside `--allowed-roots` denied, undeclared argument dropped |
| `python -m metharct.ai.mcp_server` on Python 3.8 | exit 2 with `DEPENDENCY_MISSING` naming the required interpreter |

What the MCP session actually executed: `metharct_capabilities`,
`metharct_inspect_input`, `metharct_run_analysis` (MLHuaiNanzi temperature on the
fixture, `status: succeeded`), `metharct_get_job_status`,
`metharct_get_result_summary`, `metharct_read_artifact` (page of 5 rows),
`metharct_cancel_job` on the finished job (`accepted: false`) and an unknown job
(`JOB_NOT_FOUND`).

Skipped tests (environment-limited, never faked):

1. `tests/ai/test_api_smoke.py::test_all_unavailable_modules_are_rejected` —
   at least one module is always available on both test machines.
2. `tests/ai/test_mcp_integration.py` (4 tests) — skipped on Python 3.8 only,
   because the MCP SDK requires 3.10+. They **pass** on Python 3.12.
3. `tests/ai/test_report.py::test_schemas_are_cross_checked_with_jsonschema_when_available` —
   skipped on 3.8 (no `jsonschema`); it **passes** on 3.12 with
   `jsonschema 4.26.0`, so the bundled schemas are validated by both the
   built-in validator and the reference implementation.

Platforms executed: Windows 24H2 with **Python 3.8.20** (project environment) and
**Python 3.12.7** (fresh venv with `mcp 2.3.0`, `jsonschema 4.26.0`, pandas,
numpy, scikit-learn, biopython, loguru, pyyaml). **Linux was not executed here**;
the code paths avoid POSIX-only APIs and the suite is platform-neutral, but this
must be confirmed in CI before release.

## 2. File inventory

### Added — AI layer (19 modules)

```text
metharct/ai/__init__.py        lazy public surface
metharct/ai/constants.py       schema/API versions, status vocabularies, limits
metharct/ai/clock.py           timezone-aware UTC timestamps (3.12-safe)
metharct/ai/errors.py          stable error codes + mapping
metharct/ai/paths.py           path normalisation, containment, atomic I/O, hashing
metharct/ai/schema_check.py    dependency-free JSON-Schema subset validator
metharct/ai/validators.py      streaming FASTA inspection
metharct/ai/probes.py          dependency/executable/path probes
metharct/ai/capabilities.py    single module registry + dynamic probing
metharct/ai/doctor.py          environment self-check
metharct/ai/serializers.py     native result -> unified module envelope
metharct/ai/report.py          report assembly + schema validation
metharct/ai/artifacts.py       manifest, whitelisting, pagination, path redaction
metharct/ai/jobs.py            job store + state machine
metharct/ai/runner.py          thin adapters over the existing analyzers
metharct/ai/worker.py          isolated job execution process
metharct/ai/api.py             MethArCTClient
metharct/ai/mcp_tools.py       MCP tool contracts + handlers
metharct/ai/mcp_server.py      optional FastMCP server
```

### Added — contracts, docs, knowledge, tests

```text
metharct/schemas/metharct-ai-report.schema.json
metharct/schemas/metharct-job.schema.json
metharct/schemas/metharct-capabilities.schema.json
metharct/schemas/metharct-input-inspection.schema.json
metharct/cli/ai_commands.py
docs/AI_ADAPTATION_AUDIT.md  docs/AI_ADAPTATION_REPORT.md  docs/AI_API.md
docs/MCP_TOOLS.md            docs/AI_REPORT_SCHEMA.md      docs/SECURITY.md
knowledge/README.md
knowledge/modules/{pathway,cultivation,temperature,salinity,ph,antibiotic,genome_quality}.md
knowledge/biology/{methanogenesis,sulfur_metabolism,nitrogen_metabolism,marker_gene_interpretation}.md
knowledge/interpretation/{evidence_levels,not_detected_vs_absent,prediction_vs_measurement,limitations}.md
knowledge/schemas/metharct_ai_report.md
tests/fixtures/{protein_small.faa,nucleotide_small.fna}
tests/ai/{conftest,test_clock,test_paths,test_errors,test_validators,test_capabilities,
          test_serializers,test_report,test_jobs,test_artifacts,test_cli_json,
          test_mcp_tools,test_mcp_integration,test_api_smoke}.py
tests/regression/{test_legacy_compat,test_packaging}.py
requirements-dev.txt
CHANGELOG.md
.gitignore
```

### Modified (5 files, all additive except the CLI init)

| File | Change |
|---|---|
| `metharct/__init__.py` | version `0.8.0` |
| `metharct/cli/main.py` | register AI subparsers + early dispatch for them |
| `metharct/cli/__init__.py` | eager imports -> lazy `__getattr__` (same public names) |
| `setup.py` | `ai` / `mcp` extras (MCP guarded by `python_version >= "3.10"`), `metharct-mcp` script, ship `schemas/*.json`, version fallback |
| `MANIFEST.in` | include `metharct/schemas/*.json`, `docs/*.md`, `knowledge/*.md` |
| `requirements.txt` | runtime dependencies only (was: runtime + dev/doc tooling) |

### Not modified

`metharct/core/**`, `metharct/analyzer.py`, `metharct/utils/**`,
`metharct/data/**`, `data/databases/**`, `metharct/cli/commands.py`,
`tests/test_core_logic.py`. Verified with `git status --short metharct/core data`
(only untracked `__pycache__` entries appear).

## 3. Milestone acceptance

| Milestone | Acceptance criterion | Evidence |
|---|---|---|
| M0 | audit doc + regression baseline | `docs/AI_ADAPTATION_AUDIT.md`; `tests/regression/test_legacy_compat.py` |
| M1 | API unit tests pass, core untouched | `tests/ai/test_errors.py`, `test_serializers.py`, `test_capabilities.py` |
| M2 | input detection + dependency status accurate | `tests/ai/test_validators.py`, `test_capabilities.py`, live `capabilities` run |
| M3 | state tracked, failures preserved | `tests/ai/test_jobs.py`, `test_api_smoke.py` |
| M4 | JSON schema + CLI JSON tests pass | `tests/ai/test_report.py`, `test_cli_json.py` |
| M5 | MCP tool contracts + negative tests | `tests/ai/test_mcp_tools.py` (SDK-free) **and** `tests/ai/test_mcp_integration.py` (real stdio session on Python 3.12 + MCP SDK 2.3.0) |
| M6 | knowledge base, docs, compatibility report | `knowledge/`, `docs/`, this report |

## 4. Open issues

1. **Linux not executed.** Run `pytest` and the CLI walkthrough on Linux before
   tagging the release.
2. **MCP SDK 1.x not exercised locally.** The 2.x generation
   (`mcp.server.mcpserver.MCPServer`, verified with 2.3.0) is covered end to end;
   the 1.x branch (`mcp.server.fastmcp.FastMCP`) is a two-line import fallback
   with the same decorators and is untested here. Pin `mcp<2` in CI or run the
   integration test once against a 1.x release.
3. **Cancellation granularity.** Stops between modules only; a running DIAMOND
   process is not interrupted.
4. **`jsonschema` not installed in the 3.8 project environment.** The schema
   cross-check runs only where the package is available (verified on 3.12).

## 5. Rollback

The upgrade is additive, so rollback is cheap:

```bash
# 1. Revert the five modified files
git checkout -- metharct/__init__.py metharct/cli/main.py \
                metharct/cli/__init__.py setup.py MANIFEST.in

# 2. Remove the new files/directories
rm -rf metharct/ai metharct/schemas metharct/cli/ai_commands.py \
       docs knowledge tests/ai tests/regression tests/fixtures \
       requirements-dev.txt CHANGELOG.md
# (keep .gitignore if you want the cache/workspace ignores)

# 3. Re-run the original tests
python -m pytest tests/test_core_logic.py -q
```

Nothing in `metharct/core/`, `data/`, or the legacy CLI behaviour depends on the
AI layer, so reverting those paths restores the exact 0.7.3 behaviour with no
migration step and no leftover state outside `<workspace>/jobs/`.
