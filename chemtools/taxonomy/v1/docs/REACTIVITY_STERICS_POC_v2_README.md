# Reactivity feature layer (v2): Ar- sterics (0–10) via ortho bulk proxy

Generated: 2025-12-18

## Goal

Provide a fast, topological steric descriptor around **Ar-** anchor sites that works uniformly for any **Ar-* templated feature** (Ar–Br, Ar–CHO, Ar–OTf, Ar–CN, …).

This implements `chemtools/taxonomy/v1/docs/CODEX_IMPL_ReactivitySterics_v2.md`.

## Files

- `chemtools/taxonomy/v1/specs/reactivity_features.computed.v2.json` — computed-feature spec (v2)
- `ar_context_sterics_v2.py` — Ar- analysis module (anchor extraction + ortho + bulk + score)
- `reactivity_sterics_poc_v2.py` — CLI runner

These reuse:

- `chemtools/taxonomy/v1/specs/calculable_features.compiled.v1.json`
- `chemtools/taxonomy/v1/specs/organic_groups.v1.json`

## Definition (v2)

Per detected Ar-* anchor site (Ar ipso is atom-map `:1` in templated SMARTS):

- pick the smallest aromatic ring containing the ipso atom
- ortho atoms = ring neighbors of ipso in that ring
- ortho substitution count = number of ortho atoms with `GetTotalNumHs() == 0`
- ortho bulk proxy per ortho:
  - if ortho has a non-ring substituent branch, BFS out to 4 bonds excluding ring atoms and count unique heavy atoms (cap 50)
  - else if ortho is blocked (`GetTotalNumHs()==0`), use `blocked_value=2`, else 0
- bulk score `0..8` = `min(4,o1_heavy)+min(4,o2_heavy)`
- steric score `0..10` = `0` if bulk=0 else `round(1 + 9*bulk/8)` (clamped to 10)

## Usage

```bash
# From repo root:
python chemtools/taxonomy/v1/reactivity_sterics_poc_v2.py "Cc1ccccc1Br" "Brc1ccccc1" "O=Cc1ccccc1"

# Or as a module:
python -m chemtools.taxonomy.v1.reactivity_sterics_poc_v2 "Cc1ccccc1Br"
```

Output fields include:

- `aryl_anchor_site_count`
- `aryl_ortho_sub_count_list`
- `aryl_ortho_bulk_heavy_list`
- `aryl_ortho_bulk_score_list`
- `aryl_steric_score_0_10_list`
- `aryl_steric_score_0_10_max`
- `sites` (per-site explain payload)
