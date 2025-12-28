# Codex implementation notes: classify vendor SMILES into motif/compound types (e.g., `Ar-Br`) using `organic_groups.v1.3.json` + `organic_compounds.v1.3.json`

This document describes how to use the **new** group/compound JSONs to classify a list of SMILES (vendor catalogs) into motif types such as `Ar-Br`, `Vinyl-OTf`, `Ar-B(OH)2`, etc.

**No RDKit reaction SMARTS are needed.**  
We only use **molecule SMARTS** (`MolFromSmarts` + `HasSubstructMatch`).

---

## Files and assumptions

- `organic_groups.v1.3.json`
  - top-level key: `groups` (list)
  - each group has `id`, `kind`, `smarts`, plus `tags/description`
  - **Attachpoints are encoded only by atom-maps in `smarts`**
    - scaffold groups use `:1` (e.g., `Ar = [c:1]`)
    - substituent groups use `:2` (e.g., `Br = [Br:2]`)
  - For **classification**, atom-maps are optional; we typically strip them.

- `organic_compounds.v1.3.json`
  - top-level key: `compounds` (list)
  - compound entries are either:
    1) **template-based motifs** with fields: `template`, `A`, `B`
    2) **direct SMARTS motifs** with field: `smarts` (string or list)

- Templates used in v1.3 compounds:
  - `single_bond`: `A + B`
  - `via_oxygen`: `A + "O" + B`

> Because only these two templates are used, you can hardcode them in code and you do not need `smarts_templates.v1.json` for classification.

---

## Output recommendation

For each input SMILES, return:

```json
{
  "smiles": "...",
  "ok": true,
  "hits": ["Ar-Br", "Pyridine", "..."],
  "best": "Ar-Br"
}
```

- `hits`: all matched compound motif IDs
- `best`: optional single label, chosen by a precedence function (see below)

---

## Step-by-step algorithm

### Step 1 — Load groups and compounds

1) Read JSON with UTF-8.
2) Build a dict: `groups_by_id[group.id] = group`.

### Step 2 — Compile “detect SMARTS” for each compound motif

Create `compound_id -> list[query_mol]` where each query is an RDKit SMARTS query.

For each compound entry:

- If it has `smarts`:
  - compile each SMARTS to a query
- Else if it has `template` + `A` + `B`:
  - retrieve `A_smarts = groups_by_id[A].smarts`
  - retrieve `B_smarts = groups_by_id[B].smarts`
  - compile motif SMARTS by template:
    - `single_bond`: `A_smarts + B_smarts`
    - `via_oxygen`: `A_smarts + "O" + B_smarts`

**Then strip atom-maps for detection:**

- `"[c:1][Br:2]" -> "[c][Br]"`
- `"[c:1]O[S:2](=O)(=O)..." -> "[c]O[S](=O)(=O)..."`

Why strip:

- mapping is only for attachpoint semantics; classification uses substructure matching.

### Step 3 — Classify each SMILES

For each SMILES:

1) `mol = Chem.MolFromSmiles(smiles)`
2) If `mol is None`: mark `ok=false`, `hits=[]`
3) Else:
   - for each compound id:
     - if any query matches (`mol.HasSubstructMatch(query)`):
       - add compound id to `hits`

### Step 4 — Choose an optional “best” label (if needed)

Many motifs overlap. Example:

- bromopyridine may match `Ar-Br` and `Pyridine`.

Two common policies:

**Policy A (recommended): keep `hits` only.**

- downstream logic decides based on use case.

**Policy B: compute `best` using precedence rules.**
Typical precedence heuristics:

1) optionally prefer direct SMARTS motifs (e.g., `Pyridine`) when present
2) prefer longer/more specific motif IDs
3) prefer context prefixes (`Ar-*`, `Vinyl-*`, `R-*`) over `Any-*` when both appear

---

## Reference code skeleton (minimal, production-friendly)

```python
import json, re
from rdkit import Chem

RE_MAP = re.compile(r":\d+(?=\])")  # remove atom maps like ":1" before closing bracket

def strip_atom_maps(smarts: str) -> str:
    return RE_MAP.sub("", smarts)

def load_groups(groups_path: str) -> dict:
    data = json.load(open(groups_path, "r", encoding="utf-8"))
    return {g["id"]: g for g in data["groups"]}

def load_compounds(compounds_path: str) -> list[dict]:
    data = json.load(open(compounds_path, "r", encoding="utf-8"))
    return data["compounds"]

def compile_compound_smarts(comp: dict, groups_by_id: dict) -> list[str]:
    # direct SMARTS compound (string or list)
    if "smarts" in comp:
        s = comp["smarts"]
        return s if isinstance(s, list) else [s]

    tmpl = comp["template"]
    A = groups_by_id[comp["A"]]["smarts"]
    B = groups_by_id[comp["B"]]["smarts"]

    if tmpl == "single_bond":
        return [f"{A}{B}"]
    if tmpl == "via_oxygen":
        return [f"{A}O{B}"]
    raise ValueError(f"Unknown template: {tmpl} for compound {comp.get('id')}")

def build_query_index(groups_path: str, compounds_path: str) -> dict[str, list[Chem.Mol]]:
    groups_by_id = load_groups(groups_path)
    compounds = load_compounds(compounds_path)

    index: dict[str, list[Chem.Mol]] = {}
    for comp in compounds:
        smarts_list = compile_compound_smarts(comp, groups_by_id)

        qlist: list[Chem.Mol] = []
        for s in smarts_list:
            detect = strip_atom_maps(s)
            q = Chem.MolFromSmarts(detect)
            if q is None:
                raise ValueError(f"Bad SMARTS for {comp.get('id')}: {detect}")
            qlist.append(q)

        index[comp["id"]] = qlist

    return index

def classify_smiles(smiles: str, query_index: dict[str, list[Chem.Mol]]) -> tuple[bool, list[str]]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, []
    hits = []
    for cid, qlist in query_index.items():
        if any(mol.HasSubstructMatch(q) for q in qlist):
            hits.append(cid)
    return True, hits

# Optional: best-label precedence (simple heuristic; customize)
def choose_best(hits: list[str]) -> str | None:
    if not hits:
        return None
    for h in hits:
        if h.startswith(("Ar-", "Vinyl-", "R-", "Bn-", "Allyl-")):
            return h
    return hits[0]
```

---

## Performance notes (important for large vendor catalogs)

If you have >100k SMILES:

- Precompile all SMARTS once (already done above).
- Consider an RDKit prefilter to avoid checking every query:
  - use `rdMolDescriptors.PatternFingerprint` on both molecules and queries,
  - skip query if (mol_fp & query_fp) != query_fp (subset test).
- Or use `rdkit.Chem.rdSubstructLibrary.SubstructLibrary` for batch search patterns (depends on RDKit build).

If your list is huge (millions):

- Use multiprocessing to parallelize SMILES parsing + matching.
- Persist compiled query SMARTS and/or fingerprints to disk (pickle) to avoid recompiling.

---

## Common edge cases and recommended handling

1) **Multiple matches**
   - Always keep `hits` list; don’t force a single label unless needed.

2) **Overlapping definitions (Ar vs heteroaryl motifs)**
   - If you need a single label, prefer ring-specific motifs (e.g., `Pyridine`, `Quinoline`) over generic `Ar-*`.

3) **Invalid SMILES**
   - return `ok=false` and store an error message if desired.

4) **Salts / dot-disconnected SMILES**
   - vendor SMILES may include counterions (`.[K+]` etc.)
   - RDKit `MolFromSmiles` handles this; your substructure queries should match within the major fragment.

---

## Quick test set (sanity)

- `Brc1ccccc1` should hit `Ar-Br`
- `Brc1ccncc1` should hit `Ar-Br` and `Pyridine`
- `OB(O)c1ccccc1` should hit `Ar-B(OH)2`
- `Brc1ccc(OC)cc1` should hit `Ar-Br` and also `Ar-OMe` if you have that motif

---

## Minimal “what Codex should change in existing code”

- Stop using `attach_points` / `attach_labels` (removed).
- If you previously generated detect SMARTS from `smarts_templates.v1.json` expansions:
  - switch to compiling motifs from `organic_compounds.v1.3.json` directly.
- Ensure your compiler supports:
  - template-based compounds (`single_bond`, `via_oxygen`)
  - direct SMARTS compounds (field `smarts`)
- Always strip atom-maps for detection matching.
