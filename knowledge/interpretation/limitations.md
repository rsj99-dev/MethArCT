---
doc_id: KB-INT-LIMITATIONS
title: Consolidated limitations of MethArCT
metharct_version_range: ">=0.7.3,<1.0"
module: all
source_files:
  - metharct/ai/capabilities.py
  - metharct/core/
last_verified: 2026-11-17
---

# Limitations

## Scope

* MethArCT targets **methanogenic archaea**. Predictions on other taxa are out
  of scope.
* The pathway databases cover methane, sulfur and nitrogen metabolism plus
  biosynthetic capacity. Other metabolisms are not assessed.

## Evidence

* Pathway calls are homology-based; homology is not function.
* `not_detected` ≠ absent (see `not_detected_vs_absent.md`).
* Completeness is relative to the bundled reference set and its version.
* Best-hit reporting (`max_target_seqs`) can hide secondary assignments.

## Prediction

* Temperature, salinity and pH values are model outputs, not measurements.
* Model probabilities are uncalibrated.
* Genomes outside the training distribution (flagged `is_novel` for pH) are
  less reliable.
* Partial or contaminated genomes distort every composition-based prediction.

## Cultivability

* Categories come from fixed rule thresholds on pathway completeness.
* Auxotrophy predictions require experimental confirmation.
* "Challenging" does not mean "uncultivable".

## Antibiotics

* Five reference proteomes and four fixed AAI rules; not a susceptibility test.

## Software

* CheckM2 is external; its database version affects genome quality numbers.
* Cancellation is best effort between modules, not mid-DIAMOND.
* Results depend on the reference database contents, which are fixed per
  release. Reproducing a result requires the same MethArCT version and the same
  input SHA-256 (both recorded in the report).

## What must never be claimed

* A quantitative methane, sulfide or nitrogen flux.
* A growth rate, yield, or substrate affinity.
* An experimentally confirmed phenotype.
* That a missing pathway proves the organism cannot perform that metabolism.
