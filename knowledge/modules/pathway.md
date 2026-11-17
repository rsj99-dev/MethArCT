---
doc_id: KB-MOD-PATHWAY
title: Metabolic pathway prediction (DIAMOND homology search)
metharct_version_range: ">=0.7.3,<1.0"
module: pathway
source_files:
  - metharct/core/diamond_analyzer.py
  - metharct/data/default_config.yaml
scientific_references:
  - DIAMOND protein aligner (Buchfink et al.)
  - KEGG orthology identifiers (K numbers) referenced in the FASTA headers
last_verified: 2026-11-17
---

# Pathway module

## What it does

`DiamondAnalyzer.analyze_sequence` runs DIAMOND blastp of the query proteome
against the bundled reference pathway databases and summarises, per database,
how many reference sequences were matched.

## Reference databases

| Group | Databases |
|---|---|
| Methane metabolism (16) | `CO2-CH4`, `CO-CH4`, `Acetate`, `Methanol`, `Methylamine`, `Dimethylamine`, `Trimethylamine`, `Methanethiol`, `Methoxy`, `Formate`, `Alkane`, `Fatty_acid`, `Glycine betaine methanogenesis`, `Methylthiopropionate methanogenesis`, `Tetramethylammonium methanogenesis`, `Methanol dismutation methanogenesis` |
| Sulfur metabolism (6) | `ASR`, `SO`, `SOX`, `S4I`, `SR`, `DSR` |
| Nitrogen metabolism (8) | `ANR`, `DEN`, `DNR`, `NIT`, `AMX`, `CNIT`, `NFX1`, `NFX2` |
| Cultivability | `CULTIVATION` (the cultivation pathway set, see `cultivation.md`) |

The per-database reference sequence counts are declared in
`metharct/data/default_config.yaml` (`reference_sequence_counts`). For `SR` the
denominator used at runtime is the number of *unique FASTA subjects*, not the
smaller configured number — this is covered by a regression test.

## Thresholds

DIAMOND itself is run with the configured `tools.diamond.evalue` (default
`1e-5`) and `tools.diamond.max_target_seqs` (default `1`). Pathway matching
uses the analyzer's internal identity/coverage rules, which are **not**
modified by the AI layer. Do not paraphrase a numeric threshold without reading
`diamond_analyzer.py` at the installed version.

## Outputs

* `values.pathways[]` — one entry per database with
  `detection_status` (`detected` / `not_detected` / `not_assessed`),
  `reference_sequences`, hit counts and completeness percentages.
* `evidence[]` — up to five best DIAMOND hits per threshold class with
  `reference_id`, `percent_identity`, `evalue`, `bitscore` and
  `alignment_length`.
* `values.pathways_detected`, `complete_pathways`, `raw_hit_pathways`,
  `databases_requested`, `databases_searched`.

## Limits

* Homology is not function. A hit proves sequence similarity to a reference
  entry, nothing more.
* `not_detected` means no hit passed the threshold; it is not evidence of
  absence (see `interpretation/not_detected_vs_absent.md`).
* Completeness percentages are relative to the bundled reference set, not to a
  validated minimal gene set for the pathway.
* The databases cover methane, sulfur and nitrogen metabolism. Carbon
  fixation, fermentation and other pathways are out of scope.
