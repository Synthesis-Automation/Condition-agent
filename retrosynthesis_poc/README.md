# C-X one-step retrosynthesis POC

This standalone package tests whether trusted mapped reaction observations can
generate useful one-step disconnections for products containing C-N, C-O, or
C-S bonds. It does not change `reactive_taxonomy`, `condition_registry`, or
`condition_recommender` behavior.

The POC accepts only observations with:

- a pass- or review-quality, verified or inferred reaction core;
- one reaction event and one reported product;
- verified product completeness;
- exactly one observed single C-N, C-O, or C-S bond formation;
- no bond-order or stereochemical change; and
- a retrosynthetic template that regenerates its own recorded contributing
  reactants when applied to its source product.

For an unmapped source row, the POC recomputes the deterministic internal atom
correspondence and materializes it as temporary atom maps. The mapped reaction
must then produce a pass-quality core whose center transition agrees with the
stored inferred core. Mapping evidence and confidence are retained on every
precedent. This does not promote or rewrite the source observation.

Template extraction and stereochemistry-aware application use RDChiral. Core
quality and edit selection come from the repository's structure-derived
observation contracts. Source reaction names and named families are ignored.

## Build a bounded development library

```powershell
python -m retrosynthesis_poc build-library `
  datasets/literature/shards `
  results/retrosynthesis_poc/cx_templates.json.gz `
  --include "C_N_Coupling*.jsonl.gz" `
  --include "C_O_Coupling*.jsonl.gz" `
  --include "C_S_Coupling*.jsonl.gz" `
  --max-rows 10000
```

Omit the `--include` filters to consider all shards, and omit `--max-rows` for
a full scan. Full scans read the canonical compressed observation shards
rather than the recommendation SQLite index because the shards retain mapped
reaction SMILES.

## Disconnect a target

```powershell
python -m retrosynthesis_poc disconnect `
  results/retrosynthesis_poc/cx_templates.json.gz `
  "CC(=O)c1ccccc1Nc1ccccc1" `
  --concise --top-k 5
```

With `--concise`, the command prints only proposed reaction SMILES, one per line.
Use `--top-k 5` to print the five highest-ranked proposals; omit `--concise` for
the full ranked JSON output. The `--bond` filter remains optional; for example,
add `--bond C-N` to search only C-N disconnections.

Each result contains a proposed `precursors>>target` reaction, supporting
template and precedent IDs, an interpretable score breakdown, and the outcome
of running the proposed forward reaction through the existing reaction
featurizer. The proposed reaction SMILES can subsequently be passed to the
condition recommender; this POC intentionally does not recommend conditions
until a complete precursor proposal exists.

## Current boundaries

- One-step candidate generation only; there is no stock database or route
  search.
- Exact RDChiral radius-one templates are used. Template-radius calibration and
  hierarchical correction are future experiments.
- Ranking is a declared deterministic heuristic, not a calibrated probability.
- Multiple valid retrosyntheses can exist; exact recorded-precursor recovery is
  only one evaluation metric.
