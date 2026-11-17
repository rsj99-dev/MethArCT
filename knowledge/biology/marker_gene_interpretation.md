---
doc_id: KB-BIO-MARKERS
title: Interpreting marker-gene (K number) hits
metharct_version_range: ">=0.7.3,<1.0"
module: pathway, cultivation
source_files:
  - data/databases/
  - metharct/core/diamond_analyzer.py
last_verified: 2026-11-17
---

# Marker-gene interpretation

## What the reference files contain

Each reference FASTA header carries a KEGG orthology identifier (`K#####`)
with a short description, for example:

```text
>K00399 mcrA; methyl-coenzyme M reductase alpha subunit [EC:2.8.1.245]
```

DIAMOND reports the best-matching subject id, so a hit maps onto a K number.
The pathway completeness is derived from how many distinct reference subjects
were matched.

## Rules of thumb

1. **Identity matters.** A 30 % identity hit in a deeply conserved fold is much
   weaker evidence than an 80 % identity hit in a diagnostic enzyme.
2. **One hit is not a pathway.** A single subunit of a multi-subunit complex
   (for example MCR) does not establish the whole route.
3. **Paralogs exist.** Methyl-coenzyme M reductase (MCR) and its paralogs
   (e.g. the divergent MCR-like enzymes) can cross-match; the subtype assignment
   needs more than the best hit.
4. **Best-hit bias.** `max_target_seqs` (default 1) keeps only the top hit per
   query; multi-domain proteins can therefore be assigned to one pathway only.
5. **Contamination and strain variation** shift both identity and completeness.

## What the report gives you

`evidence[]` entries contain `query_id`, `reference_id`, `percent_identity`,
`alignment_length`, `evalue` and `bitscore` for up to five best hits per
threshold class. Use them to judge whether a call is supported, instead of
quoting only the completeness percentage.

## What it does not give you

* No phylogenetic placement, no gene neighbourhood, no synteny.
* No expression, no protein abundance, no activity measurement.
