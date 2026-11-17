---
doc_id: KB-MOD-TEMPERATURE
title: Growth temperature range prediction (MLHuaiNanzi)
metharct_version_range: ">=0.7.3,<1.0"
module: temperature
source_files:
  - metharct/core/huainanzi_analyzer.py
  - metharct/core/mlhuainanzi/features.py
  - metharct/core/mlhuainanzi/predict.py
  - metharct/core/mlhuainanzi/models/model_T_{min,opt,max}.pkl
scientific_references:
  - MLHuaiNanzi AAC + dipeptide composition PLS regression models (bundled)
last_verified: 2026-11-17
---

# Temperature module

## What it does

Predicts three growth temperatures from the whole-proteome amino acid
composition:

| Value | Unit | Meaning |
|---|---|---|
| `T_min` | °C | predicted lower growth boundary |
| `T_opt` | °C | predicted optimum |
| `T_max` | °C | predicted upper growth boundary |

Features are amino-acid composition (AAC) plus dipeptide composition (DC);
three PLS regression models (`model_T_min.pkl`, `model_T_opt.pkl`,
`model_T_max.pkl`) are applied to the same feature vector. The models are
ships-with-the-package artefacts; the AI layer does not retrain or recalibrate
them.

## Category labels

The analyzer labels the prediction by `T_opt`:

| `T_opt` | Label |
|---|---|
| `< 20 °C` | `Psychrophilic` |
| `< 40 °C` | `Mesophilic` |
| `< 60 °C` | `Thermophilic` |
| `>= 60 °C` | `Hyperthermophilic` |

## Outputs

```json
{
  "values": {"T_min": 66.16, "T_opt": 80.14, "T_max": 92.7, "temperature_category": "Thermophilic"},
  "units": {"T_min": "degC", "T_opt": "degC", "T_max": "degC"},
  "model": {"name": "MLHuaiNanzi", "version": null}
}
```

## Limits

* These are predictions from sequence composition, not measured growth
  temperatures. Never present them as experimental values.
* The model was trained on a specific reference set. Genomes far outside that
  set (very small proteomes, partial MAGs, heavy contamination) give unreliable
  numbers — check `genome_quality` first.
* No calibrated confidence interval is available; the report therefore contains
  no confidence field for this module.
* Short proteomes (< a few hundred proteins) are dominated by composition noise.
