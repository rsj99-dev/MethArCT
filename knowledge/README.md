# MethArCT knowledge base

Stable, curated knowledge for humans and agents. Dynamic facts (installed
version, module availability, dependency status, job state, real results) are
**never** stored here — query `metharct capabilities`, `metharct doctor` and the
job report instead.

```text
knowledge/
├── README.md                     # this file
├── modules/                      # what each module does, inputs, outputs, limits
│   ├── pathway.md  cultivation.md  temperature.md  salinity.md
│   ├── ph.md       antibiotic.md   genome_quality.md
├── biology/                      # what the reference databases represent
│   ├── methanogenesis.md  sulfur_metabolism.md
│   ├── nitrogen_metabolism.md  marker_gene_interpretation.md
├── interpretation/               # how to read results without over-claiming
│   ├── evidence_levels.md  not_detected_vs_absent.md
│   ├── prediction_vs_measurement.md  limitations.md
└── schemas/
    └── metharct_ai_report.md
```

## Rules for this knowledge base

1. Only claims supported by the source files, the bundled configuration or the
   project authors are written here.
2. Every numeric threshold is quoted with the file that defines it.
3. Nothing here is a prediction. Examples are labelled as fixtures.
4. Documents carry front matter so they can be versioned and re-verified.

## Where to look

| Question | Document |
|---|---|
| Which module answers which biological question? | `modules/*.md` |
| What do the reference databases encode? | `biology/*.md` |
| How strong is a given result? | `interpretation/evidence_levels.md` |
| Is "not detected" the same as "absent"? | `interpretation/not_detected_vs_absent.md` |
| Is this a measurement or a prediction? | `interpretation/prediction_vs_measurement.md` |
| What can MethArCT not tell me? | `interpretation/limitations.md` |
| What does a report field mean? | `docs/AI_REPORT_SCHEMA.md` |
