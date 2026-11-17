---
doc_id: KB-SCHEMA-REPORT
title: Field guide for metharct_ai_report.json
metharct_version_range: ">=0.8.0,<1.0"
module: all
source_files:
  - metharct/schemas/metharct-ai-report.schema.json
  - metharct/ai/report.py
  - docs/AI_REPORT_SCHEMA.md
last_verified: 2026-11-17
---

# Report field quick reference

The authoritative description is `docs/AI_REPORT_SCHEMA.md`; the machine
contract is `metharct/schemas/metharct-ai-report.schema.json` (schema version
`1.0`).

## Fast lookup

| Question | Field |
|---|---|
| Which MethArCT version produced this? | `software.version` |
| Which input bytes were analysed? | `input.sha256` |
| Which modules ran? | `requested_modules`, `modules[*].status` |
| Did anything fail? | `errors`, `modules[*].status == "failed"` |
| What are the numbers? | `modules.<id>.values` |
| What are the units? | `modules.<id>.units` |
| What supports the call? | `modules.<id>.evidence` |
| What must I not over-claim? | `modules.<id>.limitations` |
| Where are the original files? | `artifacts[*].path` (relative) |
| Which thresholds were used? | `provenance.parameters` |
| Which model produced this? | `modules.<id>.model`, `provenance.model_versions` |

## Null rules

* `null` = not computed / not available. Never read it as zero.
* `status: partial` means the module ran but produced no usable value.
* `status: not_available` means the module could not run in that environment.

## Reading tables

Do not load a whole artifact. Use the artifact id and paging:

```python
page = client.read_artifact(job_id, "a001", offset=0, limit=50)
while page["next_offset"] is not None:
    page = client.read_artifact(job_id, "a001", offset=page["next_offset"], limit=50)
```

Use `filter_column` / `filter_value` to narrow rows (for example
`Detection_Status = Detected`).
