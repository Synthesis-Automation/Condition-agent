# Suzuki + Buchwald POC — Feature Engine + CLI

This folder contains a small, self-contained proof-of-concept (POC) system for:

- **Suzuki–Miyaura coupling**: sp² electrophile + organoboron partner
- **Buchwald–Hartwig amination**: sp² electrophile + N–H nucleophile

The POC is **feature-first**:

1. Compute **atomic SMARTS features** (direct RDKit substructure matches)
2. Compute **derived features** (boolean logic over feature tokens)
3. Use features to match **reactant types** (taxonomy tree) and **reaction types** (rule list)

## What’s in here

- Feature engine + validation: `chemtools/taxonomy/new/feature_engine.py`
- CLI demo: `chemtools/taxonomy/new/rxn_classifier.py`
- Attachment-point templating demo: `chemtools/taxonomy/new/demo_attachment_point_templating.py`
- Attachment-point templating evaluator: `chemtools/taxonomy/new/template_eval.py`
- Latest evaluation report: `chemtools/taxonomy/new/TEMPLATE_EVAL_REPORT.md`
- Legacy POC script (now uses the engine): `chemtools/taxonomy/new/classify_suzuki_buchwald_poc.py`
- Unit tests: `tests/test_suzuki_buchwald_poc.py`

### Data files

- Atomic SMARTS features: `chemtools/taxonomy/new/calculable_features.atomic.suzuki_buchwald.poc.json`
- Derived features: `chemtools/taxonomy/new/calculable_features.derived.suzuki_buchwald.poc.json`
- Combined convenience file: `chemtools/taxonomy/new/calculable_features.v2.suzuki_buchwald.poc.json`
- Reactant type taxonomy tree: `chemtools/taxonomy/new/reactant_types.suzuki_buchwald.poc.json`
- Reaction type rules: `chemtools/taxonomy/new/reaction_types.suzuki_buchwald.poc.json`
- Template provenance (build-time only): `chemtools/taxonomy/new/smarts_templates.suzuki_buchwald.poc.json`

## Requirements

- **RDKit** is required for SMARTS/SMILES evaluation.
- If you set `CHEMTOOLS_DISABLE_RDKIT=1`, classification/validation will be unavailable (and tests will skip).

## How the engine works

The engine implements the algorithm described in `chemtools/taxonomy/new/CODEX_TEST_BRIEF_Suzuki_Buchwald_POC.md`.

### 1) Atomic features (SMARTS)

Each atomic feature entry contains `detect.smarts_any`:

- For each SMARTS in `smarts_any`, compile with the **central cache**: `chemtools.util.smarts_cache.compile_smarts()`
- Evaluate `mol.HasSubstructMatch(pattern)`
- Result is an **OR** across `smarts_any`

API: `compute_atomic_features(mol, atomic_feature_defs) -> dict[token, bool]`

### 2) Derived features (token logic)

Derived features contain a `derive` object:

- `any_of`: OR
- `all_of`: AND
- `none_of`: NOT(OR(...))

Evaluation is dependency-aware:

- Build a derived-token dependency graph
- **Topologically sort** derived tokens
- Evaluate in order (cycles are rejected during validation)

API: `compute_derived_features(base_values, derived_feature_defs) -> dict[token, bool]`

### 3) Reactant type matching (taxonomy tree)

`reactant_types.*.json` is a tree. Each node has a `feature_token`.

- A node matches when `features[feature_token] == True`
- The matcher can return only most-specific matches (`leaf_only=True`, default)

API: `match_reactant_types(feature_values, reactant_type_tree, leaf_only=True)`

### 4) Reaction feature aggregation + reaction type matching

For a list of reactant SMILES:

- Compute features per reactant
- Aggregate to reaction-level:
  - `reaction_has[token] = any(reactant_has[token])`
  - `reaction_count[token] = sum(reactant_has[token])`

Reaction type rules (`reaction_types.*.json`) use:

- `when.all_of`
- `when.any_of`
- `when.none_of`

Matches return a `ReactionTypeMatch` with minimal “why” info (`why_all_of`, `why_any_of`, `why_none_of`).

API:

- `aggregate_reaction_features(per_reactant_features)`
- `match_reaction_types(reaction_has, reaction_types)`
- `classify_reaction_smiles(reactant_smiles, ...)` (end-to-end helper)

## Validation

`validate_features()` checks:

- No duplicate feature tokens
- All derived dependencies exist
- Derived graph is acyclic
- SMARTS compile cleanly (requires RDKit)
- Optional: `reactant_types` and `reaction_types` only reference known tokens

## Run the CLI demo

From repo root:

```bash
python -m chemtools.taxonomy.new.rxn_classifier --reactant "Brc1ccccc1" --reactant "OB(O)c1ccccc1"
python -m chemtools.taxonomy.new.rxn_classifier --reactant "Brc1ccccc1" --reactant "Nc1ccccc1"
```

Useful flags:

- `--json`: emit the full classification payload (reactant features, aggregated features, matches)
- `--no-validate`: skip SMARTS/schema validation (faster if you trust the JSON)

## Run the tests

```bash
pytest -q tests/test_suzuki_buchwald_poc.py
```

## Demonstrate attachment-point templating

This prints template expansions (and, if RDKit is available, shows attachment atom maps like `:1`
and runs a few example matches):

```bash
python -m chemtools.taxonomy.new.demo_attachment_point_templating
```

## Notes / scope

- This is a motif-based classifier (no atom mapping / bond-change detection).
- SMARTS and feature coverage are intentionally minimal for a fast POC.
