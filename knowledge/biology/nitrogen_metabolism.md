---
doc_id: KB-BIO-NITROGEN
title: Nitrogen metabolism reference sets
metharct_version_range: ">=0.7.3,<1.0"
module: pathway
source_files:
  - data/databases/nitrogen/
last_verified: 2026-11-17
---

# Nitrogen metabolism in MethArCT

| Database key | Process |
|---|---|
| `ANR` | assimilatory nitrate reduction (nitrate → ammonium for biosynthesis) |
| `DNR` | dissimilatory nitrate reduction (nitrate → ammonium for energy) |
| `DEN` | denitrification (nitrate/nitrite → NO/N₂O/N₂) |
| `NIT` | nitrification (ammonium → nitrite) |
| `CNIT` | complete nitrification (comammox, ammonium → nitrate) |
| `AMX` | anammox (anaerobic ammonium oxidation) |
| `NFX1` | nitrogen fixation, set 1 (nitrogenase components) |
| `NFX2` | nitrogen fixation, set 2 |

## How to read the result

* `ANR` and `DNR` share the nitrate-reductase chemistry but differ in fate and
  regulation; report them separately.
* `NIT` and `CNIT` overlap: a complete-nitrification genome will often also hit
  the `NIT` set.
* Anammox (`AMX`) requires a very specific gene set (hydrazine synthase and
  related genes); partial hits are common and should not be reported as anammox
  capability.
* Nitrogen fixation (`NFX1`/`NFX2`) requires nitrogenase; partial hits often
  correspond to paralogous proteins, so require high completeness before making
  a statement.

## What MethArCT does not cover

* Nitrogen flux, isotope fractionation, or rate measurements.
* Regulatory responses to nitrogen availability.
* Nitric oxide detoxification pathways outside the reference sets.
