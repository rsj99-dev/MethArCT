---
doc_id: KB-INT-EVIDENCE
title: Evidence levels in a MethArCT report
metharct_version_range: ">=0.7.3,<1.0"
module: all
source_files:
  - metharct/ai/serializers.py
  - metharct/ai/report.py
last_verified: 2026-11-17
---

# Evidence levels

MethArCT does not emit a calibrated "confidence" for pathway calls. Strength is
judged from the reported numbers and evidence, using these levels.

## Pathway calls

| Level | Observable pattern | How to phrase it |
|---|---|---|
| Strong | completeness near 100 %, several distinct reference subjects, high identity hits | "all reference genes of the route were detected" |
| Moderate | completeness 50–99 %, multiple hits | "the route is partially represented in the reference set" |
| Weak | one or two hits, low identity, low completeness | "a homology signal was observed; the route is not established" |
| None | `not_detected` | "no reference gene passed the threshold" (never "absent") |
| Not assessed | `not_assessed` | "the database was not searched" |

## Prediction modules

| Level | Observable pattern |
|---|---|
| Strong | model value present, `is_novel: false`, genome quality high, genome size typical |
| Moderate | model value present but genome partial, or `is_novel: true` |
| Weak | module `partial`, no value, unusual proteome size |
| None | module `not_available` / `failed` |

## Rules

1. Never upgrade "detected" to "active".
2. Never upgrade "predicted" to "measured".
3. Always state the genome quality context when interpreting any prediction.
4. Quote `evidence[]` identifiers (reference ids, E-values) when a claim matters.
5. If the report contains `null`, say the value is unavailable — do not
   substitute a number from memory or from another genome.
