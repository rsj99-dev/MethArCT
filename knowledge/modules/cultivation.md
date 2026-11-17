---
doc_id: KB-MOD-CULTIVATION
title: Cultivability and amino-acid auxotrophy assessment
metharct_version_range: ">=0.7.3,<1.0"
module: cultivation
source_files:
  - metharct/core/cultivation_analyzer.py
  - metharct/data/default_config.yaml
  - data/databases/cultivation/
last_verified: 2026-11-17
---

# Cultivation module

## What it does

`CultivationAnalyzer.analyze_genome` searches the query proteome against four
reference pathway collections and reports, per pathway, the fraction of
reference genes that were found (`completeness`, 0–1).

| Collection | Config key | Directory |
|---|---|---|
| Amino acids | `databases.cultivation_pathways.amino_acid` | `cultivation/path.aa` |
| Vitamins / coenzymes | `databases.cultivation_pathways.vitamin` | `cultivation/path.vc` |
| Nucleotide synthesis | `databases.cultivation_pathways.nucleotide` | `cultivation/path_hesuan` |
| ATPases | `databases.cultivation_pathways.atp` | `cultivation/path_atp` |

## Aggregation and categories

`generate_cultivability_assessment` averages the four category means and
classifies:

| Mean completeness | Category |
|---|---|
| `>= 0.8` | `high` |
| `>= 0.5` | `medium` |
| `< 0.5` | `low` |

(Thresholds: `cultivation_analyzer.py`.) The AI report copies these categories
verbatim in `values.overall_cultivability`.

`analyze_amino_acid_biosynthesis` reports, for the 20 standard amino acids,
which biosynthetic routes are complete, partial or have no dedicated pathway in
the database. The report exposes `values.amino_acid_summary.missing_list`.

## Outputs

* `values.total_pathways`, `average_pathway_completeness`,
  `amino_acid_completeness`, `vitamin_coenzyme_completeness`,
  `nucleotide_completeness`, `atp_completeness`.
* `values.pathways[]` with `pathway_name`, `completeness`, `found_genes`,
  `total_genes` and a `detection_status` where completeness `1.0` is
  `detected`.
* `values.amino_acid_summary` (`present`, `missing_or_partial`,
  `no_dedicated_pathway`, `missing_list`).

## Limits

* Completeness is relative to the bundled reference gene set. A missing gene can
  be an assembly gap, a divergent orthologue or a genuinely absent gene.
* "Missing biosynthesis" is a hypothesis about auxotrophy; it must be confirmed
  by growth experiments.
* The `confidence` value produced by the core cultivability assessment is a
  fixed heuristic tied to the number of analysed pathways. It is **not** a
  calibrated probability, and the AI report labels it as a rule heuristic.
* ATPase completeness is informative about energy metabolism, not about
  whether an organism grows in a given medium.
