# Motif-based Steric & Electronic Analysis (Groups → Compounds)

## Goal

Implement a **chemist-intuitive** analysis pipeline that:

1. Detects **compound motifs** (e.g., `Ar-Br`, `Ar-CHO`, `Ar-OTf`, `R-Br`, `R-OH`) in a molecule using `organic_compounds.*.json`.
2. For each detected motif hit, computes:
   - **Steric score** (0–10) around the **reactive center**
   - **Electronic score** (0–10) around the same reactive center (simple tag-weighted + optional Gasteiger)

This system **avoids calculable_features**. It uses only:
- `organic_groups.v1.2.json`
- `smarts_templates.v1.json`
- `organic_compounds.v1.2.json`

All outputs must be **explainable** (include per-position contributions).

---

## Critical Requirement: Analysis must be agnostic to X (functional group identity)

For motifs like `Ar-X` / `Arom-X` / `R-X` where `X` may be `Br`, `CHO`, `OH`, etc.:

- **Never branch logic on the identity of `X`** (no special cases for halides vs aldehydes vs alcohols).
- Always anchor analysis on the motif match atom maps:
  - `:1` (A-side) = reactive-center anchor (e.g., ipso aryl carbon for `Ar-*`, alpha carbon for `R-*`)
  - `:2` (B-side) = the attached group branch that defines the motif
- The only allowed “use” of X is structural:
  - when computing **scaffold** descriptors, **exclude the `:2` branch** so results are comparable across different X
  - optionally compute a separate **group bulk** for the `:2` branch, but still without recognizing what X is

Therefore, output **both**:
- `scaffold_*` scores (default, excludes `:2` branch)
- optional `group_*` descriptors (computed generically from the `:2` branch)

---

## Input Files

### 1) `organic_groups.v1.2.json`
- Contains group SMARTS and lightweight `tags`.
- Use tags for simple electronics strength:
  - `ewg_strong` → +2.0
  - `ewg_moderate` → +1.0
  - (future) `edg_moderate` → -1.0, `edg_strong` → -2.0

### 2) `smarts_templates.v1.json`
- String templates connecting group SMARTS, e.g.:
  - `single_bond`: `"{A}{B}"`
  - `via_oxygen`: `"{A}O{B}"`

### 3) `organic_compounds.v1.2.json`
- Each compound has:
  - `id` like `Ar-Br`, `Ar-CHO`, `R-OTs`, `Ar-CCH`…
  - `template`, `A`, `B`
- Generate compound SMARTS from `A/B/template`.

---

## Required Outputs

### Motif hit object

```json
{
  "compound_id": "Ar-Br",
  "match_atom_map": {"1": 12, "2": 18},
  "a_atom_idx": 12,
  "b_atom_idx": 18,
  "bond": [12, 18]
}
```

- `a_atom_idx` is always the **analysis anchor** (`:1`).
- `b_atom_idx` is the **X-branch root** (`:2`) and must be excluded for scaffold descriptors.

### Analysis result object (per hit)

```json
{
  "compound_id": "Ar-Br",
  "center": {"anchor": 12, "bond": [12, 18]},
  "steric": {
    "scaffold_score_0_10": 6.5,
    "group_bulk_0_10": 0.8,
    "method": "scaffold_agnostic_v1",
    "details": {}
  },
  "electronic": {
    "scaffold_score_0_10": 7.2,
    "including_group_score_0_10": 8.2,
    "method": "tag_weighted_scaffold_v1",
    "details": {},
    "optional": {
      "gasteiger": {"q_anchor": 0.08}
    }
  }
}
```

Notes:
- `steric.scaffold_score_0_10` excludes the `:2` branch (agnostic).
- `steric.group_bulk_0_10` measures size of the `:2` branch generically (optional).
- `electronic.scaffold_score_0_10` excludes the ipso/alpha `:2` group contribution.
- `electronic.including_group_score_0_10` includes it (optional; still no special-casing by X).

---

## Implementation Plan

### A) Build a compound SMARTS registry

Compile `compound_id → (queryMol, A, B, B_tags)`.

Steps:
1. Load groups into dict `groups[group_id]`.
2. Load templates into dict `templates[template_id]`.
3. For each compound:
   - `A_smarts = groups[A].smarts`
   - `B_smarts = groups[B].smarts`
   - `compound_smarts = templates[template].format(A=A_smarts, B=B_smarts)`
   - compile query with `Chem.MolFromSmarts`
4. Store also `B_tags = groups[B].tags`.

**Important:** preserve atom-map labels `:1` and `:2`.

### B) Motif detection (RDKit)

`detect_motifs(mol, registry) -> list[MotifHit]`

- For each compiled query, run `GetSubstructMatches(uniquify=True)`.
- Determine which query atoms have `AtomMapNum == 1` and `== 2`.
- Convert match tuple → `a_atom_idx`, `b_atom_idx`.
- De-duplicate by `(compound_id, a_atom_idx, b_atom_idx)`.

### C) Analysis dispatcher (agnostic)

Dispatch only on **context prefix** (A), never on X:

```python
if compound_id.startswith(("Ar-", "Arom-")):
    analyze_aryl(...)
elif compound_id.startswith(("R-", "Bn-", "Allyl-")):
    analyze_alkyl(...)
```

---

## Steric Analyzer (Agnostic): `scaffold_score` + optional `group_bulk`

### Common rule

- Let `anchor = a_atom_idx` and `x_root = b_atom_idx`.
- When exploring substituents to compute **scaffold steric**, never traverse into `x_root`.
- If you compute `group_bulk_0_10`, do it by exploring only the `x_root` branch away from `anchor`.

### Aryl scaffold steric (`Ar-*` / `Arom-*`)

1) Find aromatic ring containing `anchor` (choose smallest aromatic ring containing anchor).
2) Find ortho ring atoms: ring neighbors of anchor (typically 2).
3) For each ortho ring atom:
   - substituents are neighbors not in ring
   - compute bulk of each substituent by BFS outward (do not enter ring)
4) `bulk_total = bulk(o1) + bulk(o2)` (cap 20)
5) `scaffold_score_0_10 = round(10 * bulk_total / 20, 1)`

**Optional** `group_bulk_0_10`:
- BFS from `x_root` outward, do not traverse back through `anchor`.
- Compute `heavy_atoms + ring_penalty + branching` and scale to 0–10.

### Alkyl scaffold steric (`R-*` / `Bn-*` / `Allyl-*`)

Delegate to the alkyl plan (see separate MD), but enforce the same principle:
- alpha substituents = neighbors of anchor excluding `x_root`
- compute a radius-2 bulk sum + baseline by substitution degree
- optional `group_bulk_0_10` for the X branch

---

## Electronic Analyzer (Agnostic): scaffold default

This is primarily useful for `Ar-*` / `Arom-*`.

### Aryl electronics (tag-weighted)

1) Find ring containing `anchor` as in sterics.
2) For each ring atom, find substituent branches (neighbors not in ring).
3) Recognize substituent group ids minimally (v1.0):
   - NO2, CN, CHO, CO2H, CO2R, CONHR
   - Use group tags to map to strengths (`ewg_strong`/`ewg_moderate`).
4) Position weights by ring distance from anchor:
   - ortho 1.0, meta 0.6, para 1.2
5) Compute:
   - `E_scaffold` excluding the `x_root` branch attached at `anchor`
   - `E_including` including it (optional)

Convert to 0–10:
- `score = clamp(round(5 + E, 1), 0, 10)`

Output:
- `scaffold_score_0_10` (default)
- `including_group_score_0_10` (optional, but recommended to return both)

### Optional: Gasteiger

If enabled:
- compute Gasteiger charges
- report `q_anchor`

Do not make it required.

---

## Public API

```python
def analyze_smiles(smiles: str, registry_paths: dict, options: dict) -> dict:
    # {
    #   "smiles": "...",
    #   "motifs": [...],
    #   "analyses": [...]
    # }
```

Options:
- `include_gasteiger: bool`
- `electronics_include_ipso_group: bool | "both"` (default `"both"`; always compute scaffold score)
- `max_hits_per_compound: int` (optional)

---

## Tests (trend-based)

1) `Brc1ccccc1` (bromobenzene)
- detect `Ar-Br`
- steric scaffold ~0
- electronics scaffold ~5

2) `O=CC1=CC=CC=C1` (benzaldehyde)
- detect `Ar-CHO`
- steric scaffold ~0
- electronics scaffold ~5 (since scaffold excludes CHO branch)
- electronics including_group > 5

3) `O=[N+]([O-])c1ccc(Br)cc1` (p-nitro bromobenzene)
- electronics scaffold > 5 due to NO2 substituent (not the Br)

---

## Deliverables

- `motif_registry.py`, `motif_detect.py`
- `aryl_steric.py` (scaffold + optional group_bulk)
- `aryl_electronics.py` (scaffold + optional include)
- `alkyl_steric.py` (from separate plan; scaffold + optional group_bulk)
- `analyze.py` (dispatcher)
- `tests/` (pytest)
