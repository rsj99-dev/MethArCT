# MethArCT AI report schema (version 1.0)

Every job writes `metharct_ai_report.json` next to the original result files.
The machine-readable schema ships with the package:
[`metharct/schemas/metharct-ai-report.schema.json`](../metharct/schemas/metharct-ai-report.schema.json).

`schema_version` is **independent** of the package version: it only changes
when the report structure changes incompatibly.

## Job directory layout

```text
<workspace>/jobs/<job_id>/
├── job.json                  # job record  (metharct-job.schema.json)
├── metharct_ai_report.json   # unified report (this document)
├── manifest.json             # artifact manifest
├── logs/
│   ├── stdout.log            # worker progress
│   └── stderr.log            # tracebacks / tool chatter
└── artifacts/                # the original MethArCT output files
    ├── <stem>_metabolic_summary.csv
    ├── <stem>_<db>_diamond.tsv
    ├── <stem>_cultivation.csv
    ├── <stem>_amino_acid_biosynthesis.csv
    ├── <stem>_Huainanzi_Summary.tsv
    ├── <stem>_SuSha_Summary.tsv
    ├── <stem>_SuSha_Result.xlsx
    ├── <stem>_pH_Summary.tsv
    ├── <stem>_pH_Details.json
    ├── <stem>_antibiotic_selection.tsv
    └── <stem>_checkm2_summary.csv
```

## Top level

| Field | Type | Meaning |
|---|---|---|
| `schema_version` | string (const `1.0`) | report structure version |
| `software` | object | `name`, `version`, `commit` (null unless `METHARCT_COMMIT` is set) |
| `run` | object | `job_id`, `created_at`, `started_at`, `finished_at`, `status`, `duration_seconds` |
| `input` | object | `filename`, `sha256`, `format`, `type`, `sequence_count`, `total_residues` |
| `requested_modules` | array | module ids the caller asked for |
| `modules` | object | one module result per module id |
| `warnings` | array of string | non-fatal notes |
| `errors` | array of `{code, message}` | stable error codes |
| `artifacts` | array of artifact | see below |
| `provenance` | object | `parameters`, `database_versions`, `model_versions`, `environment` |

**Null policy.** A value that could not be computed is `null`. `0` is only ever
written when the program actually measured zero. Nothing is filled with a
placeholder, and no confidence number is invented.

`input.sha256` is the SHA-256 of the analysed file. The worker re-hashes the
input before running and refuses to continue if it changed since inspection.

## Module result

| Field | Type | Meaning |
|---|---|---|
| `module` | string | module id |
| `status` | enum | `success` / `partial` / `failed` / `skipped` / `not_available` / `not_applicable` |
| `result_type` | string | `prediction` / `evidence` / `assessment` / `none` |
| `summary` | string or null | one-sentence human summary |
| `values` | object | module-specific values (may contain `null`) |
| `units` | object | unit per value key (`degC`, `pH`, `percent`, `model_probability`, ...) |
| `evidence` | array | DIAMOND hits, rule evaluations, model warnings |
| `limitations` | array of string | scientific caveats that must accompany the result |
| `source_artifacts` | array of `{artifact_id, kind}` | links into `artifacts` |
| `model` | `{name, version}` | model identity (`version` is null when the code has none) |
| `error` | string or null | native error message when `status` is not `success` |
| `runtime_seconds` | number or null | wall clock time of the module |

### Status semantics

* `success` – the analyzer ran and produced a usable value.
* `partial` – the analyzer ran but produced no usable prediction (for example
  the pH model returned `null` because the genome lacks required features).
  The reason is in `summary` / `evidence`.
* `failed` – the analyzer raised; `error` carries the message.
* `not_available` – the module cannot run here (missing executable, model or
  database).
* `skipped` – the module was deliberately not run.
* `not_applicable` – the input is not suitable for this module.

Job status is derived from module statuses: any `failed` plus any `success`
gives `partially_succeeded`; only failures give `failed`; no usable result at
all gives `failed` (never `succeeded`).

### Pathway results

`values.pathways[]` entries expose the native DIAMOND summary without changing
it:

```json
{
  "database": "CO2-CH4",
  "pathway_name": "CO2 reduction methanogenesis pathway",
  "detection_status": "detected",
  "search_status": "searched",
  "reference_sequences": 15,
  "high_threshold_hits": null,
  "low_threshold_hits": 5,
  "high_completeness_percentage": null,
  "low_completeness_percentage": 33.33,
  "cultivability_hits": null
}
```

`detection_status` is one of `detected` / `not_detected` / `not_assessed`.

> `not_detected` means "no hit passed the threshold". It is **not** evidence of
> biological absence, and the module limitations say so explicitly.

### Temperature / salinity / pH

Values carry units and the model name, and never claim to be measurements:

```json
{
  "values": {"T_min": 66.16, "T_opt": 80.14, "T_max": 92.7, "temperature_category": "Thermophilic"},
  "units": {"T_min": "degC", "T_opt": "degC", "T_max": "degC"},
  "model": {"name": "MLHuaiNanzi", "version": null}
}
```

Salinity probabilities are `model_probability`, not calibrated experimental
confidence. Cultivability `confidence` from the core code is a rule heuristic
and is described as such in the limitations.

## Artifact entry

```json
{
  "artifact_id": "a001",
  "path": "artifacts/protein_small_Huainanzi_Summary.tsv",
  "media_type": "text/tab-separated-values",
  "size_bytes": 168,
  "sha256": "…64 hex chars…",
  "row_count": 1,
  "description": "temperature output file generated by the temperature module.",
  "module": "temperature",
  "filename": "protein_small_Huainanzi_Summary.tsv"
}
```

`path` is always relative to the job directory and may not start with `/` or
contain `..`. Binary artifacts (XLSX) are listed but cannot be read as text by
an agent.

## Reading a report

```python
from metharct.ai import MethArCTClient

client = MethArCTClient(workspace="./metharct_workspace")
summary = client.get_result_summary(job_id, modules=["temperature"])
report = client.get_result(job_id, validate=True)   # full report + schema check
```

Validation uses the bundled schema. Install `jsonschema` (`pip install
"metharct[ai]"`) if you also want third-party validation; the package itself
ships a dependency-free validator for the keyword subset it uses.
