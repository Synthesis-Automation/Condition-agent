# Reactivity feature layer (POC v1): Ar- electronics via Gasteiger delta charge (0–10)

Generated: 2025-12-18

## Goal

Compute a simple, fast “electron-poor” score (0–10) around **Ar-** anchor sites that works uniformly for **any Ar-* templated feature** (Ar–Br, Ar–CHO, Ar–OTf, Ar–CN, …).

Anchor sites are detected the same way as sterics: templated SMARTS encodes the Ar ipso atom as atom-map `:1` (e.g. `[c:1][Br]`, `[c:1][CX3H1](=O)`).

## Files

- `reactivity_features.electronics.computed.v1.json` — computed-feature spec (electronics)
- `ar_context_electronics_v1.py` — Ar- analysis module (anchor extraction + Gasteiger + scoring)
- `reactivity_electronics_poc_v1.py` — CLI runner

These reuse:

- `calculable_features.compiled.v1.json`
- `organic_groups.v1.json`

## Definition (POC v1)

Per Ar-* anchor site:

- compute Gasteiger charges: `rdkit.Chem.rdPartialCharges.ComputeGasteigerCharges(mol)`
- `q_ipso` = charge on the ipso atom
- choose smallest aromatic ring containing ipso (all atoms aromatic), compute `q_ring_mean`
- `delta_q = q_ipso - q_ring_mean` (or `delta_q = q_ipso` if ring not found)
- electron-poor score: `score = clamp(round(5 + k * delta_q), 0, 10)`

Default `k = 25` (tunable).

## Usage

```bash
# From repo root:
python chemtools/taxonomy/v1/reactivity_electronics_poc_v1.py "Brc1ccccc1" "COc1ccc(Br)cc1" "O=[N+]([O-])c1ccc(Br)cc1"

# Tune scaling:
python chemtools/taxonomy/v1/reactivity_electronics_poc_v1.py --k 30 "Brc1ccccc1"
```

Output fields include:

- `aryl_anchor_site_count`
- `aryl_electron_poor_score_0_10_list`, `aryl_electron_poor_score_0_10_max`
- `aryl_ipso_gasteiger_charge_list`, `aryl_ring_mean_gasteiger_charge_list`, `aryl_delta_gasteiger_charge_list`
- `sites` (per-site explain payload)
