---
doc_id: KB-MOD-PH
title: Growth pH preference prediction
metharct_version_range: ">=0.7.3,<1.0"
module: ph
source_files:
  - metharct/core/ph_analyzer.py
  - metharct/core/ph_predictor/predictor.py
  - metharct/core/ph_predictor/feature_extractor.py
  - metharct/core/ph_predictor/signal_peptide.py
  - metharct/core/ph_predictor/models/*.joblib
  - metharct/core/ph_predictor/hmm/*.joblib
scientific_references:
  - GenomeSpot pre-trained pH models (bundled)
last_verified: 2026-11-17
---

# pH module

## What it does

Predicts three pH values from genome-wide protein features:

| Value | Unit | Meaning |
|---|---|---|
| `ph_optimum` | pH | predicted growth optimum |
| `ph_min` | pH | predicted lower boundary |
| `ph_max` | pH | predicted upper boundary |

Feature extraction combines composition-style features with a signal-peptide
HMM (`hmm/hmm_signal_peptide.joblib`). Three scikit-learn pipelines
(`ph_optimum.joblib`, `ph_min.joblib`, `ph_max.joblib`) produce the values, and
per-target error models (`error_ph_*.joblib`) provide an error estimate. A
one-class SVM (`novelty_ph.joblib`) flags inputs that differ strongly from the
training distribution.

## Novelty flag

`values.is_novel` mirrors the analyzer: when `true`, the genome lies outside the
training set and the prediction is less reliable. The AI report keeps the flag
and adds the corresponding limitation string.

## Outputs

```json
{
  "values": {"ph_optimum": 7.1, "ph_max": 8.4, "ph_min": 5.9, "is_novel": false},
  "units": {"ph_optimum": "pH", "ph_min": "pH", "ph_max": "pH"},
  "model": {"name": "pH predictor", "version": null}
}
```

When the model cannot compute a value (for example a genome that lacks the
required features), the analyzer returns `value: null` with a `warning`.
The AI layer then reports the module as `partial`, echoes the warning in
`summary` / `evidence`, and never substitutes a number.

## Limits

* Predictions, not measurements.
* `is_novel: true` means the genome is outside the training distribution; treat
  the numbers as indicative only.
* Small or partial proteomes may not yield the features the model needs; the
  module then returns `partial` rather than a guess.
* The model was trained on methanogen genomes; other taxa are out of scope.
