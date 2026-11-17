---
doc_id: KB-BIO-METHANOGENESIS
title: Methanogenesis pathways represented in the reference databases
metharct_version_range: ">=0.7.3,<1.0"
module: pathway
source_files:
  - data/databases/methane/
  - metharct/core/diamond_analyzer.py
last_verified: 2026-11-17
---

# Methanogenesis in MethArCT

## Substrate-specific reference sets

MethArCT does not model a single "methanogenesis" pathway. It ships one
reference set per substrate/route, and each is evaluated independently:

| Database key | Route |
|---|---|
| `CO2-CH4` | hydrogenotrophic CO₂ reduction (including the Wolfe cycle enzymes) |
| `CO-CH4` | carboxidotrophic CO oxidation to methane |
| `Formate` | formate-driven methanogenesis |
| `Acetate` | aceticlastic methanogenesis |
| `Methanol` | methanol reduction / dismutation |
| `Methanol dismutation methanogenesis` | methanol dismutation route |
| `Methylamine` | methylamine utilization |
| `Dimethylamine` | dimethylamine utilization |
| `Trimethylamine` | trimethylamine utilization |
| `Methanethiol` | methanethiol / methyl sulfide utilization |
| `Methoxy` | methoxylated aromatic compound utilization |
| `Methylthiopropionate methanogenesis` | methylthiopropionate route |
| `Glycine betaine methanogenesis` | glycine betaine / betaine-methyl utilization |
| `Tetramethylammonium methanogenesis` | tetramethylammonium utilization |
| `Alkane` | alkane-associated methanogenesis |
| `Fatty_acid` | fatty-acid-associated methanogenesis |

## How to read the result

* A **complete** route (`low_completeness_percentage` near 100) is strong
  evidence that the encoding genes are present.
* A **partial** route can mean a divergent enzyme family, a partial genome, or
  a genuinely incomplete pathway. The report gives the numbers; it does not
  decide which.
* Several routes can be detected in one genome. That is normal: many methanogens
  are metabolically versatile.

## What MethArCT does not cover

* Hydrogenase diversity beyond the reference set.
* Electron bifurcation/confurcation stoichiometry.
* Thermodynamic feasibility, substrate concentrations, or growth yields.
* Non-methanogenic carbon metabolism.

Never convert a detected route into a quantitative methane production rate.
