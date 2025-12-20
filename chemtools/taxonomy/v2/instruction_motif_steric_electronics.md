# Motif-based Steric & Electronic Analysis (Groups → Compounds)

## Goal

Implement a **chemist-intuitive** analysis pipeline that:

1. Detects **compound motifs** (e.g., `Ar-Br`, `Ar-CHO`, `Ar-OTf`) in a molecule using `organic_compounds.*.json`.
2. For each detected motif hit, computes:
   - **Steric score** (0–10) around the **reactive center** (e.g., ortho sterics for `Ar-*`)
   - **Electronic score** (0–10) around the same reactive center (simple tag-based + optional Gasteiger)

This system **avoids calculable_features**. It uses only:

- `organic_groups.v1.2.json`
- `smarts_templates.v1.json`
- `organic_compounds.v1.2.json`

All outputs must be **explainable** (include per-position contributions).

---

## Input Files

### 1) `organic_groups.v1.2.json`

- Contains group definitions with SMARTS and **tags** (lightweight categories).
- Some key expectations:
  - Context groups (`Ar`, `Arom`, `Vinyl`, `R`, `Bn`, `Allyl`…) have tag `context`.
  - Functional groups have tags like:
    - `ewg_strong` (`NO2`, `CN`)
    - `ewg_moderate` (`CHO`, `CO2R`, `CO2H`, `CONHR`)
    - `o_h_nucleophile` (`OH`)
    - `n_h_nucleophile` (`NH2`, `NHR`)
    - `sulfonate_leaving_group` (`OTf`, `OTs`, `OMs`)
    - `halide` (`F`, `Cl`, `Br`, `I`)
- Use tags to assign electronics strength without hard-coding specific group IDs everywhere.

### 2) `smarts_templates.v1.json`

- Contains string templates for connecting groups:
  - `single_bond`: `"{A}{B}"`
  - `via_oxygen`: `"{A}O{B}"`

### 3) `organic_compounds.v1.2.json`

- Each compound entry contains at least:
  - `id` like `Ar-Br`, `Ar-CHO`, `R-OTs`, `Ar-CCH`…
  - `template` (e.g., `single_bond`, `via_oxygen`)
  - `A` (context group id) and `B` (functional group id)
- You must generate compound SMARTS by expanding `A/B/template` to final SMARTS.

---

## Required Outputs

### Motif hit object

For each detected motif in a molecule, return:

```json
{
  "compound_id": "Ar-Br",
  "match_atom_map": {"1": 12, "2": 18},
  "a_atom_idx": 12,
  "b_atom_idx": 18,
  "bond": [12, 18]
}
```

- Atom mapping labels come from the SMARTS (e.g., `:1`, `:2`).
- `a_atom_idx` is the **context anchor** (e.g., ipso aryl carbon for `Ar-*`).

### Analysis result object (per hit)

```json
{
  "compound_id": "Ar-Br",
  "center": {"ipso": 12, "bond": [12, 18]},
  "steric": {
    "score_0_10": 6.5,
    "method": "ortho_bulk_v1",
    "ortho": [
      {"ring_atom": 11, "bulk": 9, "heavy_atoms": 8, "has_ring": true, "branching": 1},
      {"ring_atom": 13, "bulk": 4, "heavy_atoms": 4, "has_ring": false, "branching": 0}
    ]
  },
  "electronic": {
    "score_0_10": 7.2,
    "method": "tag_weighted_v1",
    "including_group": true,
    "contributions": [
      {"pos": "para", "group_id": "NO2", "strength": 2.0, "weight": 1.2, "term": 2.4},
      {"pos": "meta", "group_id": "CN", "strength": 2.0, "weight": 0.6, "term": 1.2}
    ],
    "optional": {
      "gasteiger": {"q_ipso": 0.08}
    }
  }
}
```

---

## Implementation Plan

### A) Build a Compound SMARTS Registry

**Goal:** compile `compound_id → (smarts, A, B, tags_of_B)` for fast matching.

Steps:

1. Load `organic_groups` into a dict: `group_id → group_record`.
2. Load `smarts_templates` into a dict: `template_id → template_string`.
3. Load `organic_compounds` and for each entry:
   - Resolve `A_smarts = groups[A].smarts`
   - Resolve `B_smarts = groups[B].smarts`
   - Expand to `compound_smarts = template.format(A=A_smarts, B=B_smarts)`
4. Store:
   - `compound_id`
   - `compound_smarts`
   - `A`, `B`
   - `B_tags = groups[B].tags`

**Important:** SMARTS must preserve atom-map labels `:1` and `:2` from group SMARTS.

### B) Motif Detection (RDKit)

Create `detect_motifs(mol, compiled_compounds) -> list[MotifHit]`

Implementation notes:

- Use `Chem.MolFromSmarts` for each compound SMARTS (precompiled).
- Use `mol.GetSubstructMatches(query, uniquify=True)` to get matches.
- Extract the atom indices that correspond to map numbers:
  - For each atom in query, read `atom.GetAtomMapNum()`.
  - Find which matched atom index corresponds to map `1` and map `2`.
- Return a MotifHit object for each match.

Edge cases:

- Some patterns may match multiple times (multiple functional groups). Return one hit per match.
- Avoid duplicates: consider canonical `(compound_id, a_idx, b_idx)` uniqueness.

### C) Determine Analysis Mode per Motif Hit

Implement a dispatcher:

```python
def analyze_hit(mol, hit, groups_dict):
    if hit.compound_id.startswith(("Ar-", "Arom-")):
        return analyze_aryl_ipso(mol, hit, groups_dict)
    # later:
    # elif hit.compound_id.startswith(("R-", "Bn-", "Allyl-")): analyze_alkyl_alpha(...)
    # elif hit.compound_id.startswith(("Vinyl-",)): analyze_vinyl_beta(...)
```

For this task, **implement aryl only**:

- `Ar-*` and `Arom-*`

---

## Steric Analyzer: `Ar-*` ortho bulk

### 1) Find the ring and ortho atoms

Inputs:

- `ipso = hit.a_atom_idx`

Algorithm:

1. Use RDKit ring info (`mol.GetRingInfo().AtomRings()`) to find a ring containing `ipso`.
2. Choose the **smallest aromatic ring** containing `ipso` (prefer 6-membered, but use smallest length).
3. Let `ring_atoms` be the set.
4. Ortho atoms are the **two ring neighbors** of `ipso` within `ring_atoms`:
   - `ortho = [n for n in mol.GetAtomWithIdx(ipso).GetNeighbors() if n.GetIdx() in ring_atoms]`
   - Expect 2 for aromatic rings.

### 2) Compute substituent bulk for each ortho atom

For each ortho ring atom `o`:

- Find substituent neighbors **not in ring**:
  - `subs = [nbr for nbr in o.GetNeighbors() if nbr.GetIdx() not in ring_atoms]`
- If no substituent neighbor: bulk contribution = 0 (ortho-H).

For each substituent neighbor `s`, compute fragment features:

- BFS from `s` outward, **never re-enter ring atoms**.
- Collect fragment atom indices.
- Compute:
  - `heavy_atoms = count(atom atomicNum > 1 in fragment)`
  - `has_ring = fragment has any ring atom (RDKit ring info) -> bool`
  - `branching = 1 if the first fragment atom s has >2 heavy neighbors (excluding o) else 0`

Define bulk per substituent:

- `bulk = heavy_atoms + (2 if has_ring else 0) + (1 if branching else 0)`
- Cap each substituent bulk at 12.

If multiple substituents (rare at ortho): sum them.

### 3) Convert to a 0–10 score

- `bulk_total = bulk(o1) + bulk(o2)` (cap at 20)
- `steric_score_0_10 = round(10 * bulk_total / 20, 1)`

Return detailed contributions per ortho site.

---

## Electronic Analyzer: Tag-weighted ring electronics at ipso

### 1) Identify substituents on the aryl ring

Using the same `ring_atoms`:

- For each ring atom `r` in ring_atoms:
  - substituents are neighbors not in ring_atoms
  - for each substituent neighbor, map it to one of the **group ids** from `organic_groups` if possible.
- In v1.0 keep it minimal and only recognize:
  - `NO2`, `CN`, `CHO`, `CO2R`, `CO2H`, `CONHR`
- Use group tags to assign strengths.

### 2) Position classification (ortho/meta/para)

- Build ring adjacency graph restricted to ring_atoms.
- BFS shortest path distance from `ipso` to every ring atom.
- For a 6-membered ring:
  - distance 1 -> ortho
  - distance 2 -> meta
  - distance 3 -> para

Weights:

- `w_ortho = 1.0`, `w_meta = 0.6`, `w_para = 1.2`

### 3) Strength from tags

Define:

- `ewg_strong` -> +2.0
- `ewg_moderate` -> +1.0
- unknown -> 0.0

### 4) Include vs exclude the motif group itself

Provide both modes:

- `including_group=True`: include the ipso substituent group `B` (e.g., `CHO` for `Ar-CHO`).
- `including_group=False`: exclude the ipso substituent group and score only other substituents.

### 5) Convert to 0–10 score

Compute:

- `E = Σ (weight(position) * strength(group))`

Then:

- `score = clamp(round(5 + E, 1), 0, 10)`

### Optional: Gasteiger partial charge

If enabled:

- compute Gasteiger charges with RDKit
- store `q_ipso` (charge on ipso atom)
- keep tag score as primary

---

## Public API

Implement:

```python
def analyze_smiles(smiles: str, registry_paths: dict, options: dict) -> dict:
    # Return dict:
    # {
    #   "smiles": "...",
    #   "motifs": [MotifHit...],
    #   "analyses": [AnalysisResult...]
    # }
    ...
```

Options:

- `include_gasteiger: bool`
- `electronics_include_ipso_group: bool | "both"`
- `max_hits_per_compound: int` (optional safety)

---

## Tests / Examples

Create unit tests (pytest):

1) `Brc1ccccc1` (bromobenzene)

- detect `Ar-Br`
- steric score ~0 (both ortho H)
- electronic score baseline ~5

1) `O=CC1=CC=CC=C1` (benzaldehyde)

- detect `Ar-CHO`
- steric ~0
- electronic including group > 5 (CHO is EWG_moderate)

1) `CC(C)c1cc(Br)ccc1` (ortho-alkyl bromobenzene; any ortho-alkyl example)

- steric > 0 with ortho contribution

1) `O=[N+]([O-])c1ccc(Br)cc1` (p-nitro bromobenzene)

- electronic including other substituent strongly > 5

Use tolerance and check monotonic trends rather than exact numbers.

---

## Deliverables

1) `motif_registry.py`
   - loads JSONs, compiles SMARTS, holds `compiled_compounds`

2) `motif_detect.py`
   - `detect_motifs(mol, compiled_compounds)`

3) `aryl_steric.py`
   - `analyze_aryl_steric(mol, hit)`

4) `aryl_electronics.py`
   - `analyze_aryl_electronics(mol, hit, include_ipso_group=True, include_gasteiger=False)`

5) `analyze.py`
   - glue: `analyze_smiles(...)`

6) `tests/`
   - pytest suite for detection + scoring trends

---

## Coding Notes

- Use RDKit throughout.
- Keep everything fast:
  - precompile SMARTS queries
  - ring BFS only on ring atoms
- Keep constants (weights, caps) at top of file for easy tuning.
- Output must be JSON-serializable (plain dict/list/float/bool).
