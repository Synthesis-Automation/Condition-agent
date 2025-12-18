# Codex Plan — Ar- “Electron Density” Analysis (POC v1)

Generated: 2025-12-18

## Goal

Add a **simple, fast electron-density / electronic character analysis** for the **Ar- context anchor** that works with the same anchor-site workflow as sterics.

It must run on **any Ar-* site** (Ar–Br, Ar–CHO, Ar–OTf, Ar–CN, …) because the anchor is the **Ar- ipso atom** identified via **atom-map `:1`** in templated SMARTS.

Deliver a **0–10 scale** electronic score per Ar-* site plus interpretable raw values.

---

## Inputs

Assume these files exist in the working directory:

- `calculable_features.compiled.v1.json`
  - templated atomic features have `source.kind="templated"` and `source.groups=["Ar", "<core>"]`
  - templated SMARTS includes atom-map `:1` on Ar- anchor (e.g., `[c:1][Br]`, `[c:1][CX3H1](=O)`)
- `organic_groups.v1.json` for labels (`Ar-`, `-Br`, `-CHO`, `-OTf`, …)

Optional:

- reuse the existing anchor extraction code from sterics modules if available.

---

## Deliverables

Create:

1) `reactivity_features.electronics.computed.v1.json` (spec)
2) `ar_context_electronics_v1.py` (module)
3) `reactivity_electronics_poc_v1.py` (CLI)
4) `REACTIVITY_ELECTRONICS_POC_v1_README.md`

All UTF-8.

---

## Core workflow (must match sterics anchor workflow)

### A) Detect Ar-* anchor sites

Same as sterics:

1. Load `compiled.atomic`
2. For each atomic feature where:
   - `source.kind == "templated"`
   - `source.groups[0] == "Ar"`
3. Build RDKit query from `detect.smarts_any[0]`
4. Find the query atom index with `AtomMapNum == 1` (Ar- anchor)
5. For each match, map that query atom to **`ipso_atom_idx`** in the molecule

Site record fields:

- `label` (compose: `"Ar-" + "-Br" => "Ar-Br"`)
- `feature_token`
- `ipso_atom_idx`
Deduplicate by `(ipso_atom_idx, feature_token)`.

### B) Choose aromatic ring for context (optional but recommended)

To normalize “electron density” within the ring:

- pick the **smallest aromatic ring** containing `ipso_atom_idx`:
  - from `mol.GetRingInfo().AtomRings()`
  - candidates where all ring atoms `GetIsAromatic() == True`
  - choose smallest length

Return:

- `ring_atom_indices` (tuple/list) or `None` if not found

---

## Electronic descriptor options (choose ONE for POC)

Use **Option 1** (Gasteiger-based) for simplest + general solution.

### Option 1 (recommended POC): Gasteiger charge-based score

RDKit supports quick partial-charge estimation:

- `rdkit.Chem.rdPartialCharges.ComputeGasteigerCharges(mol)`

Compute per site:

- `q_ipso` = Gasteiger charge on ipso atom (`float`)
- `q_ring_mean` = mean charge across ring atoms (`float`) (if ring found)
- `delta_q = q_ipso - q_ring_mean` (if ring found) else `delta_q = q_ipso`

Interpretation (POC assumption):

- **more positive `delta_q`** → more electron-poor (more electrophilic aryl)
- **more negative `delta_q`** → more electron-rich

#### Map to 0–10 (electron-poor score)

Define an “electron-poor score” `e_poor_0_10`:

- Choose scaling `k` to expand typical `delta_q` values into the 0–10 range.
- Initial mapping (tunable):
  - `raw = 5 + k * delta_q`
  - `e_poor_0_10 = clamp(round(raw), 0, 10)`

Recommended starting constants:

- `k = 25`
  - so `delta_q = +0.20 -> ~10`, `delta_q = -0.20 -> ~0`

Also output raw values for calibration:

- `q_ipso`, `q_ring_mean`, `delta_q`

**Important:** include `k` in the spec/README so it can be tuned later.

---

## Required outputs (per SMILES)

Return JSON:

```json
{
  "smiles": "...",
  "aryl_anchor_site_count": 1,
  "aryl_electron_poor_score_0_10_list": [5],
  "aryl_electron_poor_score_0_10_max": 5,
  "aryl_ipso_gasteiger_charge_list": [-0.02],
  "aryl_ring_mean_gasteiger_charge_list": [-0.01],
  "aryl_delta_gasteiger_charge_list": [-0.01],
  "sites": [
    {
      "label": "Ar-Br",
      "feature_token": "aryl_bromide_present",
      "ipso_atom_idx": 3,
      "q_ipso": -0.02,
      "q_ring_mean": -0.01,
      "delta_q": -0.01,
      "electron_poor_score_0_10": 5
    }
  ]
}
```

For molecules with multiple Ar-* sites, lists must align with `sites` order.

---

## CLI

`python reactivity_electronics_poc_v1.py "<smiles1>" "<smiles2>" ...`

Args:

- `--compiled` default `calculable_features.compiled.v1.json`
- `--groups` default `organic_groups.v1.json`
- `--k` default `25` (scaling factor)

Print JSON to stdout.

---

## Tests (must include)

Use 3–4 sanity checks (assert ordering, not exact values, because charges can vary):

1) **bromobenzene** `Brc1ccccc1`
   - should detect Ar-Br site
   - score should be near midrange (e.g., 4–6)

2) **p-anisyl bromide** `COc1ccc(Br)cc1`
   - EDG → **more electron-rich** than bromobenzene
   - expect **lower** electron-poor score than bromobenzene

3) **p-nitrobromobenzene** `O=[N+]([O-])c1ccc(Br)cc1`
   - strong EWG → **more electron-poor**
   - expect **higher** electron-poor score than bromobenzene

4) **benzaldehyde** `O=Cc1ccccc1`
   - should detect Ar-CHO site
   - score likely somewhat electron-poor vs benzene baseline (directional check only)

Write tests as:

- simple assertions with inequalities rather than absolute thresholds if needed.

---

## Notes / future extensions (not required now)

- Replace/augment Gasteiger with substituent-class heuristic (EWG/EDG weighting) if you want more interpretability.
- Add meta/para position mapping within the ring for substituent-aware scoring.
- Provide optional “electron-rich score” as `10 - electron_poor_score`.

---

## Constraints

- RDKit only (no QM, no conformers required).
- Keep it fast; avoid expensive per-atom operations beyond Gasteiger + ring lookup.
- UTF-8 files.
