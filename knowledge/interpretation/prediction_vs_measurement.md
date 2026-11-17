---
doc_id: KB-INT-PREDICTION
title: Prediction versus measurement
metharct_version_range: ">=0.7.3,<1.0"
module: temperature, salinity, ph, cultivation
source_files:
  - metharct/core/huainanzi_analyzer.py
  - metharct/core/susha_analyzer.py
  - metharct/core/ph_analyzer.py
  - metharct/ai/serializers.py
last_verified: 2026-11-17
---

# Prediction vs measurement

MethArCT produces **predictions from sequence**, never measurements.

| Report field | Nature | Correct phrasing |
|---|---|---|
| `temperature.values.T_opt` | model prediction (°C) | "predicted optimum growth temperature" |
| `salinity.values.salinity_label` | model class | "predicted salinity adaptation class" |
| `salinity.values.confidence` | model probability | "model probability", not "confidence" |
| `ph.values.ph_optimum` | model prediction | "predicted pH optimum" |
| `cultivation.values.overall_cultivability` | rule-based category | "assessed cultivability category" |
| `genome_quality.values.completeness` | tool estimate | "CheckM2-estimated completeness" |
| `pathway.values.pathways[].detection_status` | homology evidence | "detected at the configured threshold" |

## Why the distinction matters

* A prediction has a domain of applicability. Outside it, the number can be
  systematically wrong even though it is precise.
* A measurement has units tied to a protocol (medium, temperature, strain
  state); a prediction has no protocol.
* Comparing a predicted `T_opt` with a literature growth optimum requires
  knowing the strain, the medium and the cultivation history.

## Required phrasing in downstream text

* Always name the model: "MLHuaiNanzi predicted…", "SuSha predicted…",
  "the pH model predicted…".
* Always include units from `units`.
* If the module status is `partial`, say no value is available instead of
  quoting a number.
* If `is_novel` is true, state that the genome is outside the training
  distribution.
* Never write "measured", "determined experimentally" or "optimum is" for a
  prediction field.
