---
doc_id: KB-BIO-SULFUR
title: Sulfur metabolism reference sets
metharct_version_range: ">=0.7.3,<1.0"
module: pathway
source_files:
  - data/databases/sulfur/
last_verified: 2026-11-17
---

# Sulfur metabolism in MethArCT

| Database key | Process |
|---|---|
| `ASR` | assimilatory sulfate reduction (sulfate → sulfide for biosynthesis) |
| `DSR` | dissimilatory sulfate reduction (energy-conserving) |
| `SR` | sulfur reduction (elemental sulfur / polysulfide) |
| `SO` | sulfide oxidation |
| `SOX` | sulfur oxidation via the SOX system |
| `S4I` | sulfur oxidation via the tetrathionate intermediate pathway |

## How to read the result

* `ASR` and `DSR` are biologically distinct even though both reduce sulfate.
  Do not merge them into one "sulfate reduction" statement.
* `SR` uses the number of **unique reference subjects** as the denominator, so a
  genome that matches every configured entry can still be below 100 % if the
  reference file contains more unique subjects than the configured count. This
  behaviour is intentional and regression-tested.
* `SOX` and `S4I` are alternative sulfur-oxidation routes; detecting both is
  possible and does not indicate a contradiction.
* A detected `DSR` set does not prove sulfate is the electron acceptor used in
  situ — that requires geochemical context.

## What MethArCT does not cover

* Sulfur isotope fractionation.
* Thiosulfate/polysulfide speciation.
* Kinetic or thermodynamic modelling of sulfur oxidation.
