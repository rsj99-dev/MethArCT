---
doc_id: KB-MOD-ANTIBIOTIC
title: Antibiotic resistance prediction from AAI
metharct_version_range: ">=0.7.3,<1.0"
module: antibiotic
source_files:
  - metharct/core/antibiotic_analyzer.py
  - data/databases/kangshengsu/*.faa
last_verified: 2026-11-17
---

# Antibiotic module

## What it does

Computes the average amino-acid identity (AAI) between the query proteome and
five methanogen reference proteomes, then applies fixed selection rules.

Reference proteomes (`data/databases/kangshengsu/`):
`Methanobacterium.faa`, `Methanomassiliicoccus.faa`, `Methanopyrus.faa`,
`Methanosarcina.faa`, `Methanothrix.faa`.

## Rules (verbatim from `AntibioticAnalyzer.ANTIBIOTIC_RULES`)

| Recommendation | Condition | Logic |
|---|---|---|
| Bacitracin | AAI > 50 % vs `Methanobacterium.faa` **or** `Methanopyrus.faa` | any |
| Bacitracin + Tunicamycin | AAI > 50 % vs `Methanothrix.faa` **or** `Methanomassiliicoccus.faa` | any |
| Vanadate | AAI > 65 % vs `Methanosarcina.faa` | any |
| Tunicamycin | AAI < 50 % vs **all five** references | all |

AAI is computed from DIAMOND blastp hits; the thresholds above are hard-coded
in the analyzer and are not configurable through the AI layer.

## Outputs

```json
{
  "values": {
    "recommended_antibiotics": ["Bacitracin"],
    "aai_percent": {"Methanobacterium.faa": 61.4, "Methanopyrus.faa": 48.9},
    "rules_matched": ["Bacitracin"]
  },
  "units": {"aai_percent": "percent_identity"}
}
```

`evidence[]` contains one entry per evaluated rule with its per-reference
detail strings, so the reasoning is auditable.

## Limits

* Fixed thresholds against five reference genomes: this is a heuristic
  selection aid, not a susceptibility test.
* AAI is averaged over hits; horizontal gene transfer, partial genomes and
  contamination distort it.
* "No recommendation" is not "resistant to everything". It means no rule fired.
* The `Tunicamycin` rule (`AAI < 50 %` with all references) describes a
  reference-set pattern, not a resistance mechanism.
