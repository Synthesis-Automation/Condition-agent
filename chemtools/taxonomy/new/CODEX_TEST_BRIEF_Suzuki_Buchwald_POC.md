# Suzuki + Buchwald POC (SMARTS Features + Reaction Typing) — Codex Test Brief

This repo bundle is a **proof-of-concept** feature system for:

- **Suzuki–Miyaura coupling** (sp² electrophile + organoboron partner)
- **Buchwald–Hartwig amination** (sp² electrophile + N–H nucleophile)

It is designed to be **feature-first**:

1) Compute **atomic SMARTS features** (direct substructure matches)
2) Compute **derived features** (boolean logic over tokens)
3) Use derived features to classify **reactant types** and **reaction types**


## Files in the ZIP

### Features

- `calculable_features.atomic.suzuki_buchwald.poc.json`
  - Atomic features only: `{token, type, description, detect:{smarts_any:[...]}}`
- `calculable_features.derived.suzuki_buchwald.poc.json`
  - Derived features only: `{token, type, description, derive:{any_of/all_of/none_of:[tokens...]}}`
- `calculable_features.v2.suzuki_buchwald.poc.json`
  - Convenience combined file (atomic + derived)

### Taxonomy and reaction typing

- `reactant_types.suzuki_buchwald.poc.json`
  - Tree taxonomy where each node has `feature_token`
- `reaction_types.suzuki_buchwald.poc.json`
  - Reaction-type definitions using feature constraints

### SMARTS templating provenance

- `smarts_templates.suzuki_buchwald.poc.json`
  - Fragments with attachment labels and templates used to generate repetitive SMARTS families
  - Used **at build-time** only; runtime uses generated SMARTS from atomic features

### Python scripts

- `build_suzuki_buchwald_templates_poc.py`
  - Shows how templating expands fragments + templates into atomic SMARTS entries
- `classify_suzuki_buchwald_poc.py`
  - CLI POC to classify reactions from a set of reactant SMILES

## Conceptual Data Model

### Atomic feature

Detects a structural motif directly via SMARTS:

```json
{
  "token": "aryl_bromide_present",
  "type": "bool",
  "description": "Aryl bromide (Ar–Br)",
  "detect": { "smarts_any": ["[c:1][Br]"] }
}
```

### Derived feature

No SMARTS; combines other tokens:

```json
{
  "token": "sp2_electrophile_present",
  "type": "bool",
  "description": "Any sp2 electrophile (aryl/vinyl halide or pseudohalide)",
  "derive": { "any_of": ["aryl_halide_present", "vinyl_halide_present", "aryl_pseudohalide_present", "vinyl_pseudohalide_present"] }
}
```

## Runtime Evaluation Algorithm (what to implement/test)

### Step 1 — Atomic evaluation

Input: SMILES (single molecule) → RDKit Mol  
For each atomic feature:

- compile SMARTS once (cache)
- evaluate: `mol.HasSubstructMatch(smarts)`
- if `smarts_any` list: OR across patterns  
Output: `results[token] = True/False`

### Step 2 — Derived evaluation

Evaluate derived tokens after atomics. Implement one of:

- **Topological sort** based on dependencies, then evaluate in order (preferred)
- or iterative fixpoint evaluation (works if no cycles)

Derived operators:

- `any_of`: OR
- `all_of`: AND
- `none_of`: NOT(OR(...))

### Step 3 — Reactant type matching

Traverse `reactant_types.*.json` tree:

- node matches if `results[node.feature_token] == True`
- optionally return **all matched nodes** or only **leaf matches**

### Step 4 — Reaction type matching (Suzuki/Buchwald)

For a set of reactant SMILES (list):

- compute features per reactant, then aggregate to reaction-level:
  - `reaction_has[token] = any(reactant_has[token])`
  - keep also counts: `reaction_count[token] = sum(...)` for stricter matching
- classify Suzuki if:
  - has sp² electrophile AND has organoboron partner
- classify Buchwald if:
  - has sp² electrophile AND has N–H nucleophile (amine / aniline / amide-with-NH, depending on POC rules)

## What to Ask Codex to Build (Test Harness)

### A) Feature engine (library)

Create a small Python module, e.g. `feature_engine.py`, that exposes:

- `compute_atomic_features(mol, atomic_feature_defs) -> dict[token,bool]`
- `compute_derived_features(base_dict, derived_feature_defs) -> dict[token,bool]`
- `compute_all_features(mol, all_defs) -> dict[token,bool]` (atomic + derived)
- caching for SMARTS patterns

### B) Validation / lint

Implement `validate_features(atomic_defs, derived_defs)`:

- all referenced tokens exist
- derived graph has **no cycles**
- no duplicate token ids
- SMARTS compile cleanly
- optional: detect obviously redundant SMARTS (same pattern repeated)

### C) Unit tests (pytest)

1) **Atomic feature tests**  
Use a small table of SMILES expected to match specific tokens.
2) **Derived feature tests**  
Construct synthetic token dictionaries and verify logic.
3) **Reaction classification tests**

- Suzuki positive:
  - bromobenzene + phenylboronic acid
- Suzuki negative controls:
  - bromobenzene + aniline (should be Buchwald, not Suzuki)
  - phenylboronic acid + aniline (no electrophile)
- Buchwald positive:
  - bromobenzene + aniline
- Buchwald negative controls:
  - bromobenzene + nitrobenzene (no N–H)
  - chlorobenzene + aniline (may be true depending on sp² electrophile rules; keep explicit expectations)

4) **Template provenance sanity**

- Ensure that generated features in `atomic` include the same SMARTS as in the template expansion list (where applicable)

### D) CLI runner

A simple CLI like:

```bash
python -m rxn_classifier --reactant "Brc1ccccc1" --reactant "OB(O)c1ccccc1"
```

Outputs:

- matched reaction type(s)
- matched key tokens that caused the match (explainability)

## Canonical Test SMILES (suggested)

### Electrophiles (sp²)

- Bromobenzene: `Brc1ccccc1`
- Chlorobenzene: `Clc1ccccc1`
- Iodobenzene: `Ic1ccccc1`
- Vinyl bromide (simple): `C=CBr`  (note: depends on your vinyl SMARTS)

### Pseudohalides

- Phenyl triflate: `O=S(=O)(O)c1ccccc1C(F)(F)F` (some representations vary; treat as optional test)

### Boron partners (Suzuki)

- Phenylboronic acid: `OB(O)c1ccccc1`
- Vinylboronic acid: `C=CB(O)O` (may need alternate)

### Amines (Buchwald)

- Aniline: `Nc1ccccc1`
- Morpholine: `O1CCNCC1` (tertiary? actually secondary with N-H? morpholine has N-H if unsubstituted; SMILES is `O1CCNCC1`, yes N-H)
- Triethylamine (no N–H): `CCN(CC)CC` (negative control)

## Notes / Expected POC Limitations

- POC is coarse: it classifies based on **reactant motifs**, not atom-mapped bond changes.
- Tautomer/protonation handling is minimal; consider standardization later.
- Vinyl/allyl patterns can be tricky; allow some tests to be marked xfail until tuned.

## Deliverable

Please implement:

1) a robust feature engine + validator
2) a pytest suite that tests the above scenarios
3) a CLI that prints classification + reasoning tokens

Keep everything UTF-8.

## Reference Implementation (this repo)

- Feature engine + validator: `chemtools/taxonomy/new/feature_engine.py`
- CLI demo: `chemtools/taxonomy/new/rxn_classifier.py`
- Unit tests: `tests/test_suzuki_buchwald_poc.py`

### Run the CLI (from repo root)

```bash
python -m chemtools.taxonomy.new.rxn_classifier --reactant "Brc1ccccc1" --reactant "OB(O)c1ccccc1"
python -m chemtools.taxonomy.new.rxn_classifier --reactant "Brc1ccccc1" --reactant "Nc1ccccc1"
```

### Run the tests

```bash
pytest -q tests/test_suzuki_buchwald_poc.py
```
