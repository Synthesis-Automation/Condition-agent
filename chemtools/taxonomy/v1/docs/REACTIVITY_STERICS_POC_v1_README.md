# Reactivity feature layer (POC v1): Ar- sterics via ortho substitution counts

Generated: 2025-12-18

## Goal

Provide a **systematic, reusable** steric descriptor for the **Ar- context group** that works for:

- **Ar–Br** (aryl bromides)
- **Ar–CHO** (aromatic aldehydes)
- **Ar–OTf** (aryl triflates)
- any other **Ar-* templated feature**

The key is that templated SMARTS encodes the **Ar- anchor atom** as **atom-map `:1`** (e.g. `[c:1][Br]`, `[c:1][CX3H1](=O)`).
We detect anchor sites from those templated features, then compute sterics around the anchor.

## Files

- `chemtools/taxonomy/v1/specs/reactivity_features.computed.v1.json` — computed-feature specification (steric layer)
- `ar_context_sterics_v1.py` — Ar- analysis module (anchor extraction + ortho count)
- `reactivity_sterics_poc_v1.py` — CLI runner

These require:

- `chemtools/taxonomy/v1/specs/calculable_features.compiled.v1.json`
- `chemtools/taxonomy/v1/specs/organic_groups.v1.json`

## Definition (POC)

For each detected Ar-* anchor site:

1. pick the smallest aromatic ring containing the ipso atom
2. take the two ring neighbors of the ipso atom as **ortho positions**
3. an ortho position is counted as **substituted/blocked** if `GetTotalNumHs() == 0`

This counts:

- true ortho substituents (e.g. o-Me, o-CF3, o-Cl)
- fused ring blocking
- aromatic hetero atoms (pyridine N) as blocked

## Usage

```bash
# From repo root:
python chemtools/taxonomy/v1/reactivity_sterics_poc_v1.py "BrC1=CC=CC=C1C" "O=Cc1ccccc1" "Brc1ccc(cc1)C"

# Or from this folder:
python reactivity_sterics_poc_v1.py "BrC1=CC=CC=C1C" "O=Cc1ccccc1" "Brc1ccc(cc1)C"

# Use the direct implementation (bypass computed-feature JSON spec):
python chemtools/taxonomy/v1/reactivity_sterics_poc_v1.py --mode direct "BrC1=CC=CC=C1C"
```

Output fields:

- `aryl_anchor_site_count`
- `aryl_ortho_sub_count_list`
- `aryl_ortho_sub_count_max`
- `aryl_ortho_sub_count_sum`
- `sites`: per-site label + ipso atom index + ortho count

## Next extensions (later)

- ortho bulk proxies (heavy atom count, MW, etc.)
- meta/para patterns
- electronic descriptors (EWG/EDG flags near ipso)
