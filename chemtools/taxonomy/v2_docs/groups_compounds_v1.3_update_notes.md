# Update notes: organic_groups.v1.3.json + organic_compounds.v1.3.json

Purpose: keep the definition JSONs **simple to maintain**, while still making **attachpoints explicit** for motif compilation and downstream typing (e.g., Ar‑Br, Ar‑OTf, Ar‑B(OH)2).

---

## 0) Naming System Overhaul (R vs Alkyl)

As of version 1.3, we have refined the naming of generic groups and alkyl groups:

- **Alkyl**: Now specifically refers to an $sp^3$ carbon scaffold (`[CX4:1]`).
- **R**: Redefined as a generic wildcard for **any non-hydrogen group/atom** (`[!#1:1]` for scaffold, `[!#1:2]` for substituent).
- **RCH2, R2CH, R3C**: Maintained with the "R" prefix because these refer to the degree of substitution on a central alkyl carbon (where R = anything except H).

This ensures that `Alkyl-X` specifically means a saturated carbon linkage, while `R-X` is the most generic non-H classification.

---

## 1) Key design decision

**Attachpoint is encoded only via atom-mapping inside `group.smarts`.**  
No extra fields are needed.

- **Scaffold groups** (e.g., `Ar`, `Vinyl`, `Alkyl`): attach atom is mapped as `:1`
  - Example: `Ar.smarts = "[c:1]"`

- **Substituent groups** (e.g., `Br`, `Cl`, `OTf`, `B(OH)2`, `CHO`, `CO2H`, `NH2`): attach atom is mapped as `:2`
  - Example: `Br.smarts = "[Br:2]`

> Practical implication: code can always find the attach atom by scanning SMARTS for `:1` or `:2`.

---

## 2) Breaking changes (what to update in code)

### 2.1 `organic_groups`

Removed fields (no longer present):

- `attach_points`
- `attach_labels`

ID Changes:

- The old `R` ID (which matched `[CX4]`) has been renamed to `Alkyl`.
- A new `R` ID has been created as the generic non-hydrogen wildcard.

New expectation:

- Every attachable group must have its attach atom mapped in `smarts`:
  - `:1` for scaffold groups
  - `:2` for substituent groups

### 2.2 `organic_compounds`

Renamed IDs:

- Motifs previously using `R-` as a prefix for alkyl groups have been renamed to `Alkyl-` (e.g., `R-Br` -> `Alkyl-Br`).
- Motifs using `RCH2`, `R2CH`, `R3C`, or `R_acidic` are NOT renamed.

Removed field where it was redundant:

- `anchors` was always `{"core":"A","fg":"B"}` and is now omitted.

New expectation:

- If `anchors` is missing, assume:
  - `anchors.scaffold = "A"`
  - `anchors.substituent = "B"`

No other compound fields changed.

---

## 3) Concrete group changes (examples)

### 3.1 Halides now have explicit `:2`

Before (v1.2):

- `Br.smarts = "[Br]"` (attachpoint ambiguous)

After (v1.3):

- `Br.smarts = "[Br:2]"`
- `Cl.smarts = "[Cl:2]"`
- `I.smarts  = "[I:2]"`
- `F.smarts  = "[F:2]"`

### 3.2 Sulfonate leaving groups now anchor on sulfur (`:2`)

Before (v1.2):

- `OTf.smarts = "S(=O)(=O)C(F)(F)F"` (no attach map)

After (v1.3):

- `OTf.smarts = "[S:2](=O)(=O)C(F)(F)F"`
- `OTs.smarts = "[S:2](=O)(=O)c1ccc(C)cc1"`
- `OMs.smarts = "[S:2](=O)(=O)C"`

### 3.3 Many other groups already had a mapped attach atom

Examples (unchanged mapping behavior, but legacy fields removed):

- `Ar.smarts = "[c:1]"`
- `B(OH)2.smarts = "[B:2](O)O"`
- `CHO.smarts = "[CX3H1:2](=O)"`
- `CO2H.smarts = "[C:2](=O)O"`
- `NH2.smarts = "[NH2:2]"`
- `Alkyl.smarts = "[CX4:1]"` (formerly `R`)
- `R.smarts = "[!#1:1]"` (new generic wildcard)
- `R_Subst.smarts = "[!#1:2]"` (new generic substituent wildcard)

---

## 4) How motif SMARTS compilation should work now

### 4.1 Templates used by compounds

The compound entries reference template IDs (currently `single_bond` and `via_oxygen`).

Recommended minimal implementation (hardcode or keep in a template registry):

- `single_bond`: `"{A}{B}"`
- `via_oxygen`: `"{A}O{B}"`

### 4.2 Compile “motif SMARTS” from groups + templates

Given a compound:

- look up group A: `A_smarts = groups[A_id].smarts`
- look up group B: `B_smarts = groups[B_id].smarts`
- format with template

Example: `Ar-Br`

- `Ar.smarts = "[c:1]"`
- `Br.smarts = "[Br:2]"`
- `single_bond` ⇒ `"[c:1][Br:2]"`

Example: `Ar-OTf`

- `Ar.smarts = "[c:1]"`
- `OTf.smarts = "[S:2](=O)(=O)C(F)(F)F"`
- `via_oxygen` ⇒ `"[c:1]O[S:2](=O)(=O)C(F)(F)F"`

### 4.3 Detection SMARTS (recommended)

For molecule typing / vendor catalogs, you typically **do not need atom maps**.
So generate a “detect SMARTS” by stripping atom-map labels:

- `"[c:1][Br:2]"` → `"[c][Br]"`
- `"[c:1]O[S:2](=O)(=O)..."` → `"[c]O[S](=O)(=O)..."`

Implementation note: stripping can be done with a regex like `r":\d+(?=\])"`.

---

## 5) Classification workflow for a vendor SMILES list

For each SMILES:

1) `mol = Chem.MolFromSmiles(smiles)`
2) For each compound motif:
   - compile detect SMARTS (as above)
   - `q = Chem.MolFromSmarts(detect_smarts)`
   - `if mol.HasSubstructMatch(q): record hit`

Return:

- `hits` (list of motif IDs like `["Ar-Br"]`)
- optional `best_hit` chosen by your precedence rules (e.g., prefer heteroaryl motifs like `Pyridine` over generic `Ar-*`).

---

## 6) Codex checklist (code updates)

1) **Group loader**
   - stop reading `attach_points` / `attach_labels` (they are gone)
   - if you need an attach atom, parse it from `group.smarts` (`:1` or `:2`)

2) **Compound loader**
   - allow `anchors` to be absent; default to `{"core":"A","fg":"B"}`

3) **SMARTS compiler**
   - compile motif SMARTS using `template.format(A=A_smarts, B=B_smarts)`
   - provide a helper to generate detect SMARTS by stripping atom-maps

4) **Sanity checks**
   - ensure each group `smarts` contains:
     - `:1` for scaffold groups (kind == "context")
     - `:2` for substituent groups (kind != "context")
   - ensure all template IDs used in compounds exist in your templates registry

---

## 7) Notes on smarts_templates.v1.json

If your code currently reads a separate template registry, you can keep that concept:

- only `templates` is essential (single_bond / via_oxygen formats)
- the `expansions` section is optional and can be considered build-time provenance

Recommended simplification:

- either hardcode the 2 templates in code, **or**
- keep a small templates JSON, but drop/ignore `expansions` if they duplicate your curated `organic_compounds` list.
