
# Implementation Plan: SchemeConditionDB Matcher (No Featurization)

**Date:** 2025-10-03  
**Owner:** (to be assigned)  
**Target reaction type (v1):** Buchwald C–N

---

## 0) Goal (what we’re building)

A tiny engine that:

1) Takes a **reaction SMILES** (`Reactants>>Products`) for a known reaction type (start with Buchwald C–N).
2) Runs a small, fixed **Essential-Core Normalization (ECN)** on the **reactant side** (no amine-class featurization).
3) **Substructure-matches** against a SchemeConditionDB (each entry has **reactant-only SMARTS**).
4) Returns the **best matching scheme’s** condition bundle (with sensible fallbacks), plus an **explainable trace**.

---

## 1) Deliverables

- Python package `scdb_matcher/` (RDKit required).
- JSON dataset: `buchwald_scheme_db.match.json` (provided).
- CLI: `scdb match --rxn "<RXN_SMILES>" --db buchwald_scheme_db.match.json`.
- Optional REST microservice (FastAPI) with a `/match` endpoint.
- Unit tests + a quick regression suite with sample reactions.

---

## 2) Data: SchemeConditionDB (kept simple)

**File:** `buchwald_scheme_db.match.json`

```json
{
  "schema_version": "match-1.0",
  "reaction_type": "Buchwald_CN",
  "updated_at": "ISO8601",
  "match_contract": { "how_to_match": [ "...(docs only)..." ] },
  "preprocess": {
    "strip_counterions": true,
    "neutralize": true,
    "normalize_aromaticity": true,
    "masking": [
      {"center": "[c]-[Cl,Br,I]", "keep_radius": 1, "collapse_to": "[*]"},
      {"center": "[N;H2,H1]-[*]", "keep_radius": 1, "collapse_to": "[*]"}
    ]
  },
  "entries": [
    {
      "id": "M-BH-ARBR-ANILINE",
      "type": "scheme",
      "name": "ArBr + aniline (primary)",
      "reactant_smarts": ["[c]-Br","[NH2]-[c]"],
      "priority": 40,
      "conditions": {
        "pd_source": ["Pd2(dba)3","XPhos Pd G3"],
        "ligands": ["SPhos","XPhos","JohnPhos"],
        "base": ["K2CO3","K3PO4","NaOtBu"],
        "solvent": ["toluene","1,4-dioxane"],
        "temperature_C": [80,100],
        "time_h": [6,14],
        "loadings": {"Pd_mol%":[1.0,2.0], "ligand_mol%":[2.0,4.0], "base_eq":[1.5,2.5]}
      }
    },
    {
      "id": "M-BH-DEFAULT-ArBr",
      "type": "default_condition",
      "name": "Default for aryl bromide",
      "reactant_smarts": ["[c]-Br"],
      "priority": 5,
      "default_condition": {
        "pd_source":"Pd2(dba)3","ligand":"SPhos","base":"K2CO3","solvent":"toluene",
        "temperature_C":90,"time_h":12,"concentration_M":0.1,
        "loadings":{"Pd_mol%":1.0,"ligand_mol%":2.0,"base_eq":2.0}
      }
    }
    // ... (other schemes + defaults already in the file)
  ]
}
```

**Key choices**

- **reactant_smarts** is a list; **all** must match among reactants.
- **priority** ranks specificity (e.g., ortho-hindered patterns get higher numbers).
- **defaults** also use reactant_smarts (e.g., `[c]-Cl`) and fire when no scheme matches.

---

## 3) Essential-Core Normalization (ECN) – minimal + fixed

Run this pipeline on the **reactant side** (no product usage for match):

1) **Split reactants** by `.`  
2) **Strip counterions/solvents** (remove inorganic ions and the smallest fragments that are purely inorganic).  
3) **Neutralize common salts** (optional), **normalize aromaticity** (RDKit: Kekulize → SetAromaticity).  
4) **Masking pass (radius = 1)** using a small fixed anchor list (no “classification” required):
   - Anchors:
     - Electrophile: `"[c]-[Cl,Br,I]"` and `"[c]-OS(=O)(=O)[C,F]"` (sulfonates)
     - N-donors: `"[NH2]-[c]"`, `"[NH]([CX4])[CX4]"`, `"[N;H1]C(=O)O"`
   - For each occurrence: **keep** anchor atoms + all atoms within **1 bond**, collapse everything else to `[*]` termini (to avoid false negatives from remote junk).
5) **Canonical SMARTS bag**: from each masked fragment, emit a small SMARTS (or reuse the anchors above).  
6) **Proceed to matching** with the bag (effectively: “does the reactant set contain each required SMARTS?”).

This is enough to make the system robust to big side chains, protecting groups, and random substituents — while still **never** writing a featurizer.

---

## 4) Matching algorithm (deterministic)

**Inputs:** `rxn_smiles`, `db`  
**Output:** `MatchResult` with `scheme_id` or `default_id`, `conditions`, `trace`

**Steps:**

1) Parse `rxn_smiles` → `reactants`: list of RDKit `Mol`.
2) Apply ECN (above) to get a **SMARTS bag** and **masked reactant mols**.
3) For each **scheme** entry:
   - For each `smarts` in `reactant_smarts`:
     - Substructure match across **all** reactant mols (order-independent).
   - If **all** required SMARTS matched → add to candidates with:
     - `priority` (higher is better)
     - `specificity_score` = total atom count across its SMARTS (longer = more specific)
4) If candidates:
   - Sort by `(priority DESC, specificity_score DESC, id ASC)` and pick top-1.
   - Return its `conditions`.
5) Else, evaluate **default** entries in the same way (they may have generic SMARTS like `[c]-Cl`), pick the best.
6) If still none: return **global default** (the “match anything” default in the DB).
7) Include a **trace** JSON: which SMARTS matched where, which entries were considered, and tie-breaks.

---

## 5) Code skeleton (Python/RDKit)

```python
# scdb_matcher/loader.py
import json
from rdkit import Chem

def load_db(path):
    db = json.load(open(path, "r", encoding="utf-8"))
    entries = db["entries"]
    # Precompile SMARTS patterns
    for e in entries:
        pats = []
        for s in e.get("reactant_smarts", []):
            m = Chem.MolFromSmarts(s)
            if not m:
                raise ValueError(f"Bad SMARTS in {{e['id']}}: {{s}}")
            pats.append((s, m))
        e["_pats"] = pats
    return db
```

```python
# scdb_matcher/ecn.py
from rdkit import Chem
from rdkit.Chem import rdMolStandardize

ANCHORS = ["[c]-[Cl,Br,I]","[c]-OS(=O)(=O)[C,F]","[NH2]-[c]","[NH]([CX4])[CX4]","[N;H1]C(=O)O"]
ANCHOR_MOLS = [(s, Chem.MolFromSmarts(s)) for s in ANCHORS]

def essential_core_normalize(reactant_smiles_list):
    mols = [Chem.MolFromSmiles(s) for s in reactant_smiles_list if s]
    # strip inorganics (simple heuristic)
    mols = [m for m in mols if m and any(a.GetAtomicNum() > 6 for a in m.GetAtoms())]
    # normalize aromaticity
    for m in mols: Chem.SanitizeMol(m)
    # collect masked fragments (radius=1 around anchors), else fall back to original
    masked = []
    for m in mols:
        used = False
        for pat_str, pat in ANCHOR_MOLS:
            if m.HasSubstructMatch(pat):
                for match in m.GetSubstructMatches(pat):
                    # expand atoms within radius=1
                    keep = set()
                    for idx in match:
                        keep.add(idx)
                        atom = m.GetAtomWithIdx(idx)
                        for nbr in atom.GetNeighbors():
                            keep.add(nbr.GetIdx())
                    # build a copy with outside collapsed to dummy atoms
                    em = Chem.PathToSubmol(m, list(keep), useQuery=True)
                    masked.append(em)
                    used = True
        if not used:
            masked.append(m)
    # build a simple SMARTS bag from anchors themselves (fallback)
    smarts_bag = [s for s,_ in ANCHOR_MOLS if any(m.HasSubstructMatch(p) for _,p in ANCHOR_MOLS)]
    return masked, smarts_bag
```

```python
# scdb_matcher/matcher.py
from rdkit import Chem
from .ecn import essential_core_normalize

def match(db, rxn_smiles):
    reactants_str = rxn_smiles.split(">>")[0]
    reactant_smiles = [s for s in reactants_str.split(".") if s]
    masked_mols, smarts_bag = essential_core_normalize(reactant_smiles)

    def entry_matches(e):
        for s, pat in e.get("_pats", []):
            if not any(m.HasSubstructMatch(pat) for m in masked_mols):
                return False
        return True

    schemes = [e for e in db["entries"] if e["type"]=="scheme" and entry_matches(e)]
    if schemes:
        ranked = sorted(
            schemes,
            key=lambda e: (
                e.get("priority", 0),
                sum(pat[1].GetNumAtoms() for pat in e.get("_pats", [])),
                e["id"]
            ),
            reverse=True
        )
        top = ranked[0]
        return { "match_type":"scheme", "id":top["id"], "conditions":top["conditions"],
                 "trace":{"considered":[e["id"] for e in ranked]} }

    defaults = [e for e in db["entries"] if e["type"]=="default_condition" and entry_matches(e)]
    if defaults:
        ranked = sorted(defaults, key=lambda e: (e.get("priority", 0), len(e.get("_pats", [])), e["id"]), reverse=True)
        top = ranked[0]
        return { "match_type":"default", "id":top["id"], "conditions":top["default_condition"],
                 "trace":{"considered":[e["id"] for e in ranked]} }

    globals_ = [e for e in db["entries"] if e["type"]=="default_condition" and not e.get("reactant_smarts")]
    if globals_:
        top = max(globals_, key=lambda e: (e.get("priority", 0), e["id"]))
        return { "match_type":"global_default", "id":top["id"], "conditions":top["default_condition"],
                 "trace":{"considered":[e["id"] for e in globals_]} }
    raise RuntimeError("No scheme and no default matched.")
```

**CLI sketch**

```bash
scdb match --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1" --db buchwald_scheme_db.match.json --json
```

Expected output (JSON):

```json
{
  "reaction_type": "Buchwald_CN",
  "matched_entry": "M-BH-ARBR-ANILINE",
  "conditions": { "...": "..." },
  "trace": { "considered": ["M-BH-ARBR-ANILINE"] }
}
```

---

## 6) Tests (must pass)

1) **Happy path** — `Brc1ccccc1.Nc1ccccc1>>...` → `M-BH-ARBR-ANILINE`  
2) **ArCl + aniline** — `Clc1ccccc1.Nc1ccccc1>>...` → `M-BH-ARCL-ANILINE`  
3) **Ortho-hindered** — `Clc1cc(ccc1C)C.Nc1ccccc1>>...` → `M-BH-ARCL-ANILINE-ORTHO`  
4) **Amide N–H** — `Clc1ccccc1.NC(=O)OC>>...` → `M-BH-ARCL-AMIDE`  
5) **No match** — falls back to per-halide default or global.

---

## 7) Performance & robustness

- **Precompile SMARTS** on DB load.
- **SubstructLibrary** prefilter if DB grows.
- **Cache** parsed mols per input reaction.
- **Prefer SMARTS** over SMILES in scheme patterns (aromatic stability).

---

## 8) Edge cases & rules

- Multiple amines/halides → any reactant can satisfy each SMARTS.
- Protecting groups/long chains → removed by ECN masking (radius=1 keep).
- Multiple matches → resolve by **priority** then **specificity** (total SMARTS atoms).
- Defaults: per-halide first; then global fallback.

---

## 9) Acceptance criteria

- Matches the provided test set with correct scheme IDs.
- No featurization logic (no amine class detectors).
- Returns an explainable **trace** for every decision.

---

## 10) Next steps (optional)

- Per-entry `preprocess` overrides when special masking is needed.
- Extend to **Suzuki** and **Ullmann** by the same pattern (reactant-only SMARTS).
- Add a **screening constructor** to expand a bundle into 4–8 orthogonal recipes.
