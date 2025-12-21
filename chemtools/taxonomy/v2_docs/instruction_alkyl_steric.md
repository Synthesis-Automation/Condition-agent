# Alkyl Steric Analysis (Agnostic to X) — Plan (Groups → Compounds)

## Goal

Add a **fast, chemist-intuitive steric analyzer** for **alkyl-context motifs** such as:

- `R-Cl`, `R-Br`, `R-I`, `R-OTs`, `R-OMs`, `R-OTf`
- `Bn-Br`, `Bn-Cl` (benzylic)
- `Allyl-Br`, `Allyl-OTs` (allylic)
- (optionally later) `R-OH`, `R-SH`, `R-NH2` for nucleophile sterics

This complements **aryl ortho steric** analysis.  
We avoid `calculable_features` and use motif detection from `organic_compounds.*.json`.

---

## Critical Requirement: Alkyl steric analysis must be agnostic to X

For motifs like `R-X`, `Bn-X`, `Allyl-X` where `X` can be `Br`, `CHO`, `OH`, etc.:

- **Do not branch logic on the identity of `X`.**
- Always anchor on the motif match:
  - `alpha = a_atom_idx` (map `:1`)
  - `x_root = b_atom_idx` (map `:2`)
- The default steric score must be a **scaffold score** computed from substituents around `alpha` **excluding** the `x_root` branch.
- Optionally compute `group_bulk_0_10` for the `x_root` branch generically (size of the X branch), without identifying what X is.

---

## Required Output (per motif hit)

```json
{
  "compound_id": "R-Br",
  "center": {"alpha_c": 7, "bond": [7, 15]},
  "steric": {
    "scaffold_score_0_10": 6.0,
    "group_bulk_0_10": 0.5,
    "method": "alkyl_alpha_scaffold_v1",
    "classification": "secondary",
    "alpha": {
      "degree": 3,
      "carbon_subs": 2,
      "benzylic": false,
      "allylic": false,
      "has_beta_branching": true
    },
    "substituents": [
      {"sub_atom": 6, "bulk": 3, "heavy_atoms_r2": 3, "branching": 0, "has_ring": false},
      {"sub_atom": 8, "bulk": 5, "heavy_atoms_r2": 4, "branching": 1, "has_ring": false}
    ]
  }
}
```

---

## Step 1 — Decide which motifs use Alkyl Steric

Dispatcher rule:

- if `compound_id` starts with `R-`, `Bn-`, or `Allyl-` → run alkyl steric analyzer
- `Ar-*` / `Arom-*` handled elsewhere

---

## Step 2 — Alpha-carbon classification (methyl/primary/secondary/tertiary)

Let:

- `alpha = hit.a_atom_idx`
- `x_root = hit.b_atom_idx` (mapped `:2`)

Compute:

- `degree = alpha_atom.GetDegree()` (includes X branch)
- `carbon_subs = number of carbon neighbors of alpha excluding x_root`

Classification:

- methyl if `carbon_subs == 0`
- primary if `carbon_subs == 1`
- secondary if `carbon_subs == 2`
- tertiary if `carbon_subs >= 3`

Flags:

- `benzylic = True` if any neighbor (excluding x_root) is aromatic
- `allylic = True` if alpha is adjacent to a C=C (simple: any neighbor has a DOUBLE bond to carbon)

---

## Step 3 — Scaffold steric score (0–10)

### 3.1 Alpha substituents (scaffold-only)

- `subs = [nbr for nbr in alpha.GetNeighbors() if nbr.GetIdx() != x_root]`

### 3.2 Bulk for a substituent (radius 2, v1)

For each substituent neighbor `s`:

Local set within radius 2 (from alpha via s):

- include `s`
- include neighbors of `s` excluding `alpha`
- do not traverse into the X branch (`x_root`)

Compute:

- `heavy_atoms_r2 = count(atomNum > 1 in local set)`
- `branching = 1 if s has >=2 heavy neighbors excluding alpha else 0`
- `has_ring = True if any atom in local set is in a ring`
- `bulk = heavy_atoms_r2 + (2 if has_ring else 0) + branching` (cap 12)

### 3.3 Combine and scale

Baseline by substitution class:

- `baseline = {methyl:0, primary:2, secondary:4, tertiary:7}`

Then:

- `bulk_total = sum(sub.bulk)`
- `raw = baseline + bulk_total` (cap 20)
- `scaffold_score_0_10 = round(10 * raw / 20, 1)`

### Optional: neopentyl bonus (recommended)

Detect:

- alpha is primary
- the beta carbon (the only carbon substituent) is tertiary

Then `raw += 3` (cap 20).

---

## Step 4 — Optional group bulk for X branch (agnostic)

Compute `group_bulk_0_10` from the `x_root` branch:

- BFS starting at `x_root`, do not traverse back through `alpha`
- compute `heavy_atoms`, `has_ring`, `branching`
- map to 0–10 with a simple cap (same style as other bulk metrics)

This is optional and must not require recognizing what X is.

---

## Step 5 — Tests (trend-based)

1) `CCl` methyl chloride: very low
2) `CCBr` ethyl bromide: > methyl
3) `CC(C)Br` isopropyl: > ethyl
4) `CC(C)(C)Cl` tert-butyl: near top
5) `CC(C)(C)CBr` neopentyl: > normal primary (if bonus)
6) `c1ccccc1CBr` benzyl: ring penalty present
7) `C=CCBr` allyl: low/moderate, allylic flag true

Assertions:

- tert-butyl > isopropyl > ethyl > methyl
- neopentyl > ethyl
- benzyl includes ring penalty

---

## Deliverables

- `alkyl_steric.py`
- update dispatcher in `analyze.py`
- `tests/test_alkyl_steric.py`
