# Unified Recommendation Index Guide

This system builds a single index that merges reaction datasets, protocols, and HTE screens.
DRFP fingerprints are computed for entries with reaction SMILES and stored alongside
feature-tag metadata for feature-only entries.

## Build

```bash
python data-processor/build_unified_recommendation_index.py
```

## Output

```
build/unified_recommendation_index/
  index.jsonl         # metadata + tags for each entry
  fingerprints.npz    # DRFP fingerprints (if available)
  stats.json          # counts + build summary
```

## Notes

- If DRFP is not installed, the builder still creates `index.jsonl` and `stats.json`.
- HTE entries are matched by feature tags (reaction type + reactant types).
