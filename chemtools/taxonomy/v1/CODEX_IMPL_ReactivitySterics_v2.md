# Codex Implementation Brief — Reactivity Sterics v2 (Ar- anchor, 0–10 score)

## Goal
Implement a **reactivity feature layer** that computes a fast, topological steric metric around **Ar-** sites for any compound SMILES.

It must work **uniformly** for any Ar-* pattern (Ar–Br, Ar–CHO, Ar–OTf, Ar–CN, etc.) by anchoring on the **Ar- ipso atom**.

Outputs should include:
- per-site **ortho substitution count** (0/1/2)
- per-site **ortho bulk proxy** (heavy-atom counts)
- per-site **steric score 0–10**
- aggregate **max steric score** across sites

## Inputs (files)
Use these v1 files (assume they exist in working directory):
- `calculable_features.compiled.v1.json`
  - contains `atomic` list; templated atomic features have `source.kind="templated"` and `source.groups=["Ar", "<core>"]`
  - templated SMARTS contains atom-map label `:1` on Ar- anchor atom (e.g., `[c:1][Br]`, `[c:1][CX3H1](=O)`)
- `organic_groups.v1.json`
  - group labels, `id` -> `name` (chemist style): `Ar-`, `-Br`, `-CHO`, `-OTf`, ...

Optional (not required for sterics):
- `reactant_types.v1.json`

## Deliverables
Create:
1) `reactivity_features.computed.v2.json` (spec file, see below)
2) `ar_context_sterics_v2.py` (module)
3) `reactivity_sterics_poc_v2.py` (CLI)
4) `REACTIVITY_STERICS_POC_v2_README.md`

All files UTF-8.

## Definitions

### A. Anchor site detection (Ar- site)
For each **templated atomic feature** in `compiled.atomic` with:
- `source.kind == "templated"`
- `source.groups[0] == "Ar"`

Do:
- compile SMARTS from `feature.detect.smarts_any[0]` (assume single SMARTS for templated)
- find the query atom index where `AtomMapNum == 1` (this is the Ar- ipso in the query)
- `mol.GetSubstructMatches(query)` gives matches; map `anchor_qidx` → `ipso_atom_idx` in each match

Each match yields an AnchorSite record:
- `label` = compose chemist labels: `"Ar-" + "-Br" => "Ar-Br"` (remove duplicated dash)
  - `left = groups["Ar"].name` e.g. `"Ar-"`
  - `right = groups[core].name` e.g. `"-Br"`, `"-CHO"`, `"-OTf"`
- `feature_token` = atomic feature token (e.g., `aryl_bromide_present`, `aromatic_aldehyde_present`)
- `ipso_atom_idx` = int

Deduplicate by `(ipso_atom_idx, feature_token)`.

### B. Pick ring + ortho atoms
For each `ipso_atom_idx`:
- choose the **smallest aromatic ring** containing `ipso`:
  - use `mol.GetRingInfo().AtomRings()`
  - candidates where `ipso` in ring and **all ring atoms aromatic**
  - pick smallest length ring
- ortho atoms = ring neighbors of ipso (neighbors whose idx is in ring set)
- typical count is 2; handle non-2 gracefully

### C. Ortho substitution count
Per site:
- `ortho_sub_count = number of ortho atoms with GetTotalNumHs() == 0`
Return `0..2` for normal phenyl.

### D. Ortho bulk proxy (fast, topological)
For each ortho atom `o`:
- find branch starts = neighbors of `o` that are **not in ring_set**
- if branch starts exist:
  - for each branch start, BFS outward excluding ring atoms
  - count **unique heavy atoms** visited (atomic num > 1)
  - limit BFS by:
    - `max_depth = 4` bonds from branch start (tunable)
    - `max_atoms = 50` safety cap
  - `ortho_bulk_heavy = max(branch_heavy_counts)`
- else (no external branch):
  - if `o.GetTotalNumHs() == 0` (blocked by fusion/hetero), set `ortho_bulk_heavy = blocked_value = 2`
  - else set 0

Return per site: `ortho_bulk_heavy_pair = [o1_heavy, o2_heavy]` (pad missing with 0).

### E. Bulk score (0..8) and steric score (0..10)
- cap each ortho heavy count: `cap_each = 4`
- `bulk_score = min(4, o1_heavy) + min(4, o2_heavy)` → `0..8`
- steric score mapping:
  - if `bulk_score == 0`: `steric_0_10 = 0`
  - else: `steric_0_10 = round(1 + 9 * bulk_score / 8)`
  - clamp to `<= 10`

### F. Output format
For each SMILES:
```json
{
  "smiles": "...",
  "aryl_anchor_site_count": 1,
  "aryl_ortho_sub_count_list": [0],
  "aryl_ortho_bulk_heavy_list": [[0, 0]],
  "aryl_ortho_bulk_score_list": [0],
  "aryl_steric_score_0_10_list": [0],
  "aryl_steric_score_0_10_max": 0,
  "sites": [
    {
      "label": "Ar-Br",
      "feature_token": "aryl_bromide_present",
      "ipso_atom_idx": 3,
      "ortho_sub_count": 0,
      "ortho_bulk_heavy": [0, 0],
      "ortho_bulk_score": 0,
      "steric_score_0_10": 0
    }
  ]
}
```

## computed-feature spec file (`reactivity_features.computed.v2.json`)
Create a JSON containing:
- version string
- list of computed tokens:
  - `aryl_anchor_site_count`
  - `aryl_ortho_sub_count_list`
  - `aryl_ortho_bulk_heavy_list`
  - `aryl_ortho_bulk_score_list`
  - `aryl_steric_score_0_10_list`
  - `aryl_steric_score_0_10_max`
- include mapping/parameters in `notes` (cap_each=4, max_depth=4, blocked_value=2, mapping formula)

## CLI script
`python reactivity_sterics_poc_v2.py "<smiles1>" "<smiles2>" ...`

Arguments:
- `--compiled` default `calculable_features.compiled.v1.json`
- `--groups` default `organic_groups.v1.json`

Print JSON to stdout.

## Tests (must include)
Add minimal tests (pytest or assertions), at least:
1) `Cc1ccccc1Br` (o-tolyl bromide)
   - should detect at least one site
   - `steric_score_0_10_max > 0`
2) `Brc1ccccc1` (bromobenzene)
   - `steric_score_0_10_max == 0`
3) `O=Cc1ccccc1` (benzaldehyde)
   - should detect Ar-CHO site
   - `steric_score_0_10_max == 0`

## Constraints
- Must use RDKit only (no conformer generation).
- Keep it fast; avoid O(N^2) per molecule.
- Always UTF-8.
