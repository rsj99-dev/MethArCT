---
doc_id: KB-INT-NOT-DETECTED
title: "Not detected" is not "absent"
metharct_version_range: ">=0.7.3,<1.0"
module: pathway, cultivation, antibiotic
source_files:
  - metharct/ai/serializers.py
  - metharct/core/diamond_analyzer.py
last_verified: 2026-11-17
---

# `not_detected` vs `absent`

In a MethArCT report, `not_detected` means:

> no reference sequence of that pathway passed the configured similarity
> threshold in this query.

It does **not** mean the pathway is biologically absent. Reasons a real pathway
can be `not_detected`:

1. **Divergent sequence** — the enzyme is present but below the identity or
   coverage threshold.
2. **Incomplete genome** — the encoding gene is in an assembly gap. Check
   `genome_quality.completeness`.
3. **Reference gap** — the bundled reference set does not contain that variant.
4. **Wrong reading frame / partial ORF** — the predicted protein is truncated.
5. **Database scope** — the pathway is not represented in the bundled databases
   at all (for example carbon fixation outside the methane/sulfur/nitrogen sets).
6. **Search not performed** — reported separately as `not_assessed`, which is a
   third, distinct state.

## Correct phrasing

| Good | Bad |
|---|---|
| "The CO₂-reduction route did not pass the detection threshold." | "The organism cannot reduce CO₂." |
| "No gene of the pathway was found in this assembly." | "The pathway is absent from the species." |
| "The database was not searched for this pathway (`not_assessed`)." | "The pathway was not detected." |

## How to strengthen a negative call

* Check genome completeness first.
* Look for the specific marker gene by name or by a more sensitive search.
* Compare with a close relative's genome.
* Treat a `not_detected` call as a hypothesis, not a conclusion.
