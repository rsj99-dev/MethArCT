---
doc_id: KB-MOD-SALINITY
title: Salinity adaptation prediction (SuSha)
metharct_version_range: ">=0.7.3,<1.0"
module: salinity
source_files:
  - metharct/core/susha_analyzer.py
  - metharct/core/susha/predict.py
  - metharct/core/susha/models/*.pkl
scientific_references:
  - SuSha ensemble learning model (bundled)
last_verified: 2026-11-17
---

# Salinity module

## What it does

Predicts a salinity adaptation class from genome-wide amino acid composition
features using the bundled SuSha ensemble model, and returns the model
probabilities for the top three classes.

## Classes

| Index | Label | Conventional NaCl range used in the batch exporter |
|---|---|---|
| 0 | `Salt-sensitive` | 0–1 % |
| 1 | `Halotolerant` | 1–3 % |
| 2 | `Slight halophilic` | 3–5 % |
| 3 | `Moderate halophilic` | 5–15 % |
| 4 | `Extreme halophilic` | 15–30 % |

The index/label mapping is defined in `susha_analyzer.py`
(`SALINITY_LABEL_MAP`). The percentage ranges are the ones used by the
`metharct batch` CSV exporter; they describe classes, not measured optima.

## Bundled estimators

`metharct/core/susha/models/` contains 13 estimators (AdaBoost, Bagging,
DecisionTree, ExtraTrees, GradientBoosting, KNN, LDA, LogisticRegression,
NaiveBayes, RandomForest, SVM and two ensembles:
`Ensemble1_ET_RF_Bag`, `Ensemble2_ET_LR_LDA`).

## Outputs

```json
{
  "values": {
    "salinity_label": "Slight halophilic",
    "salinity_index": 2,
    "confidence": 0.62
  },
  "units": {"confidence": "model_probability"},
  "evidence": [
    {"rank": 1, "label": "Slight halophilic", "probability": 0.62},
    {"rank": 2, "label": "Halotolerant", "probability": 0.28},
    {"rank": 3, "label": "Salt-sensitive", "probability": 0.07}
  ],
  "model": {"name": "SuSha", "version": null}
}
```

`confidence` is the maximum class probability of the model. It is reported as
`model_probability`, never as calibrated experimental confidence.

## Limits

* The class is a prediction from composition; it is not a measured NaCl
  optimum or range.
* Probabilities are not calibrated: a value of 0.62 does not mean "62 % of the
  time this is correct".
* Only methanogen genomes were used for training; applying the model to other
  archaea or bacteria is out of scope.
* SHAP interpretation is optional and non-blocking; when `shap` is missing the
  prediction still runs but no per-feature contributions are available.
