---
doc_id: KB-MOD-GENOME-QUALITY
title: Genome quality assessment (CheckM2)
metharct_version_range: ">=0.7.3,<1.0"
module: genome_quality
source_files:
  - metharct/core/checkm2_analyzer.py
scientific_references:
  - CheckM2 (external tool and database; not bundled)
last_verified: 2026-11-17
---

# Genome quality module

## What it does

Runs the external CheckM2 `predict` command on the input proteome and reports
completeness, contamination, an aggregate quality score and quality buckets.

Quality buckets used by the analyzer (`checkm2_analyzer.py`):

| Bucket | Completeness | Contamination |
|---|---|---|
| high | `>= 90 %` | `<= 5 %` |
| medium | `>= 50 %` | `<= 10 %` |
| low | below the medium thresholds | – |

## Requirements

* The `checkm2` executable, resolved from `tools.checkm2.path` (or
  `tools.checkm2.wsl_path` when `use_wsl` is true).
* A CheckM2 reference database. MethArCT does not bundle it; `metharct doctor`
  reports `CHECKM2DB` when it is set.

When either is missing the module cannot run, and the report says
`not_available` with the reason — it never returns `0 %` completeness.

## Outputs

```json
{
  "values": {
    "completeness": 98.2,
    "contamination": 1.4,
    "quality_score": 96.8,
    "total_genomes": 1,
    "high_quality_genomes": 1,
    "medium_quality_genomes": 0,
    "low_quality_genomes": 0
  },
  "units": {"completeness": "percent", "contamination": "percent"}
}
```

## Limits

* Every CheckM2 number depends on the external database version; record it if
  the result is used in a publication.
* Genome quality conditions all other predictions: a 60 % complete MAG will give
  unreliable temperature, salinity, pH and pathway results.
* CheckM2 operates on protein FASTA input in MethArCT's integration; nucleotide
  input must be called through a codon-aware workflow first.
