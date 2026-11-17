# MethArCT AI Adaption Audit (Milestone 0)

This document records the *real* state of the repository before the AI-ready
upgrade. Nothing here is inferred; every value was read from the checkout or
printed by the interpreter.

| Item | Value |
|---|---|
| Audit date | 2026-11-17 |
| Branch | `main` |
| Commit | `bf77997995796bc346a09c42da3f7a1c44481c9d` |
| Working tree at audit start | clean (`git status --short` printed nothing) |
| Python | 3.8.20 (Windows 24H2, PowerShell) |
| Package version | `0.7.3` (`metharct/__init__.py`) |
| README version | none pinned in README; the README links to the hosted `v0.6.0` beta (stale, left untouched) |
| Entry point | `metharct=metharct.cli.main:main` (`setup.py`) |

## 1. Directory tree (pre-upgrade)

```text
MethArCT/
├── data/databases/          # reference data, unchanged by this upgrade
│   ├── cultivation/{path.aa,path.vc,path_hesuan,path_atp}
│   ├── diamond/*.dmnd       # built lazily from the reference FASTA files
│   ├── kangshengsu/*.faa    # 5 methanogen reference proteomes
│   ├── methane/ sulfur/ nitrogen/
├── metharct/
│   ├── analyzer.py           # MethArCTAnalyzer
│   ├── cli/{main.py,commands.py,__init__.py}
│   ├── core/
│   │   ├── diamond_analyzer.py  cultivation_analyzer.py
│   │   ├── susha_analyzer.py    huainanzi_analyzer.py
│   │   ├── ph_analyzer.py       checkm2_analyzer.py
│   │   ├── antibiotic_analyzer.py  pathway_predictor.py
│   │   ├── mlhuainanzi/  ph_predictor/  susha/   # embedded models
│   ├── data/default_config.yaml
│   └── utils/{config.py,file_utils.py,logger.py,sequence_utils.py}
├── tests/test_core_logic.py
├── requirements.txt  environment.yml  setup.py  MANIFEST.in
├── metharct_config.yaml  quick_test.py
└── diamond.exe
```

## 2. Existing Python API

| Class | Method | Signature |
|---|---|---|
| `MethArCTAnalyzer` | `analyze_comprehensive` | `(input_file, output_dir=None, run_diamond=True, run_checkm2=True, run_antibiotic=True, **kwargs)` |
| `MethArCTAnalyzer` | `analyze_diamond` / `analyze_checkm2` / `analyze_antibiotic` | `(input_file, output_dir=None)` |
| `PathwayPredictor` | `predict_comprehensive` | `(input_path, output_prefix=None, include_diamond=True, include_checkm2=True, include_cultivation=True, include_susha=True, include_huainanzi=True, include_ph=True, include_antibiotic=True)` |
| `DiamondAnalyzer` | `analyze_sequence` | `(input_file, output_prefix=None, databases=None)` |
| `CultivationAnalyzer` | `analyze_genome` | `(genome_file)` |
| `SuShaAnalyzer` | `predict_salinity` | `(input_file, output_prefix=None)` |
| `HuainanziAnalyzer` | `predict_temperature` | `(input_file, output_prefix=None)` |
| `PHAnalyzer` | `predict_ph` | `(input_file, output_prefix=None)` |
| `AntibioticAnalyzer` | `predict_antibiotics` | `(input_file, output_prefix=None)` |
| `CheckM2Analyzer` | `analyze_genome_quality` | `(input_path, output_prefix=None, input_type='fasta')` |

These signatures are frozen by `tests/regression/test_legacy_compat.py`.

Availability flags already existed on the optional analyzers:
`SuShaAnalyzer.tool_available`, `HuainanziAnalyzer.tool_available`,
`PHAnalyzer.tool_available`, `CheckM2Analyzer.tool_available`. The AI layer
reuses those semantics but probes dependencies without instantiating analyzers
(so `capabilities` has no side effects).

## 3. CLI surface

| Command | Positional | Required options | Optional options |
|---|---|---|---|
| `comprehensive` | `input` | `-o/--output` | `--skip-checkm2`, `--skip-susha`, `--skip-mlhuainanzi`, `--skip-ph`, `--skip-antibiotic`, `--threads`, `--evalue` |
| `diamond` | `input` | `-o/--output` | `--threads`, `--evalue`, `--max-target-seqs` |
| `checkm2` | `input` | – | `-o/--output`, `-t/--input-type`, `--checkm2-path`, `--threads`, `--use-wsl` |
| `cultivation` | `input` | – | `-o/--output`, `--diamond-path`, `--threads`, `--evalue`, `--use-wsl`, `--pathway-weights`, `--cultivability-threshold` |
| `susha` | `input` | `-o/--output` | – |
| `mlhuainanzi` | `input` | `-o/--output` | – |
| `ph` | `input` | `-o/--output` | – |
| `antibiotic` | `input` | `-o/--output` | `--threads`, `--evalue` |
| `batch` | `input` (directory) | `-o/--output` (CSV) | `-mags`, `--threads`, `--evalue` |

Global options: `-v/--verbose`, `-c/--config`, `--version`.

## 4. Output files produced by the existing code

| Module | Files (prefix = `<output-dir>/<input-stem>`) |
|---|---|
| pathway | `<prefix>_metabolic_summary.csv`, `<prefix>_<db>_diamond.tsv` |
| cultivation | `<prefix>_cultivation.csv`, `<prefix>_amino_acid_biosynthesis.csv` |
| temperature | `<prefix>_Huainanzi_Summary.tsv` |
| salinity | `<prefix>_SuSha_Summary.tsv`, `<prefix>_SuSha_Result.xlsx` |
| ph | `<prefix>_pH_Summary.tsv`, `<prefix>_pH_Details.json` |
| antibiotic | `<prefix>_antibiotic_selection.tsv` |
| genome_quality | `<prefix>_checkm2_summary.csv`, `<prefix>_checkm2/` |
| comprehensive | `<prefix>_comprehensive_analysis.json`, `<prefix>_analysis_summary.txt`, `<prefix>_integrated_summary.csv` |

None of these names, columns or semantics were changed by the upgrade.

## 5. Dependency and data discovery

* `Config.resolve_database_base_dir()` searches, in order: `METHARCT_DB_DIR`,
  the configured `databases.base_dir`, `<package>/data/databases`,
  `<repo>/data/databases`, `<prefix>/share/metharct/databases`.
* DIAMOND / CheckM2 executables come from `tools.diamond.path` and
  `tools.checkm2.path` (plus the `use_wsl` / `wsl_path` variants).
* `requirements.txt` mixed runtime and development dependencies
  (`pytest`, `black`, `flake8`, `mypy`, `sphinx`, ...) while `setup.py` reads
  that file as `install_requires`, so `pip install .` also installed the dev
  tools. **Resolved in 0.8.0**: `requirements.txt` now lists runtime packages
  only, the tooling moved to `requirements-dev.txt`, and
  `tests/regression/test_packaging.py` keeps it that way.
* `data/databases/diamond/*.dmnd` are built lazily from the reference FASTA
  files by `DiamondAnalyzer._prepare_databases()`.

## 6. Test baseline (before the upgrade)

| Test | Kind | Needs DIAMOND |
|---|---|---|
| `test_repo_config_uses_unified_schema` | unit | no |
| `test_legacy_cultivation_dict_is_supported` | unit | no |
| `test_decode_encoded_filename` | unit | no |
| `test_discover_pathways_decodes_names` | unit | no |
| `test_extract_missing_amino_acids_uses_analyzer` | unit | no |
| `test_sr_denominator_uses_unique_subjects` | unit | no |
| `test_zero_hit_pathways_are_kept_and_sr_cannot_reach_100` | unit | no |
| `test_antibiotic_rule_labels_match_documentation` | unit | no |
| `test_wsl_path_conversion` | unit | no |

Result before the upgrade: **9 passed**. There was no end-to-end test, no CLI
test and no test that touched the optional model modules.

## 7. Real environment at audit time

`metharct capabilities` (after the upgrade) on the audit machine reports:

| Module | Status | Reason |
|---|---|---|
| pathway | available | `diamond` on PATH, `*.dmnd` present |
| cultivation | available | 4 cultivation pathway directories populated |
| temperature | available | 3 `*.pkl` models + scikit-learn |
| salinity | available | 13 `*.pkl` models + scikit-learn |
| ph | available | 7 `*.joblib` models + `hmmlearn` |
| antibiotic | available | 5 `*.faa` reference proteomes |
| genome_quality | unavailable | CheckM2 executable not installed |

No prediction value was fabricated from this table: it only lists what the
probes found.

## 8. Adaptation mapping (plan → implementation)

| Plan requirement | Implementation |
|---|---|
| Stable Python API (M1) | `metharct/ai/api.py` → `MethArCTClient` |
| Module registry (M1) | `metharct/ai/capabilities.py` → `MODULE_REGISTRY` |
| Input validation (M2) | `metharct/ai/validators.py` → `inspect_fasta` |
| Capabilities / doctor (M2) | `capabilities.get_capabilities`, `doctor.run_doctor` |
| Job manager + worker (M3) | `jobs.py` + `worker.py` (subprocess, `shell=False`) |
| JSON report (M4) | `report.py`, `serializers.py`, `metharct/schemas/*.json` |
| CLI JSON mode (M4) | `metharct/cli/ai_commands.py` |
| MCP server (M5) | `metharct/ai/mcp_tools.py` + `mcp_server.py` (optional) |
| Artifact access (M4) | `artifacts.py` (manifest whitelist + pagination) |
| Knowledge base (M6) | `knowledge/` |
| Docs (M6) | `docs/AI_API.md`, `MCP_TOOLS.md`, `AI_REPORT_SCHEMA.md`, `SECURITY.md` |
| Tests (all) | `tests/ai/`, `tests/regression/`, `tests/fixtures/` |

## 9. Deviations from the suggested layout

* The JSON schemas live in `metharct/schemas/` instead of a repository-root
  `schemas/` directory so that an installed wheel still ships them. The plan
  explicitly allowed a different tree ("目录是建议结构").
* No file was moved, renamed or deleted. `metharct/core/` is untouched.

## 10. Known risks and open issues

| ID | Risk | Status |
|---|---|---|
| R1 | `requirements.txt` doubles as `install_requires`, so dev tools were installed for end users. | **Fixed in 0.8.0**: runtime-only `requirements.txt` + `requirements-dev.txt`, guarded by `tests/regression/test_packaging.py` |
| R2 | The pH model returns `null` when the genome lacks the required features (observed on the 3-sequence fixture). The report marks such a module `partial` instead of inventing a value. | By design |
| R3 | The legacy `comprehensive` command still prints to stdout; the AI layer never calls it (it calls the analyzers directly). | By design |
| R4 | `PathwayPredictor._integrate_results` hard-codes `confidence: 0.85` for the temperature block. The AI report does **not** reuse that number. | Reported, not changed (science code) |
| R5 | Windows only was executed for the upgrade; Linux CI has not been run here. | Open |
| R6 | Cancellation is best effort: a worker is killed between modules (`taskkill` / `SIGTERM`), never mid-DIAMOND. | Documented |
| R7 | The MCP Python SDK requires Python >= 3.10 while MethArCT supports 3.8+. The `mcp` extra is guarded by a `python_version` marker, `doctor` explains the gap, and the MCP integration test skips with that reason on 3.8. | Documented |
