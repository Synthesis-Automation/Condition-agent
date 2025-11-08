# SMARTS Alignment Audit (Phase 1)

Date: 2025-11-08  
Owner: Codex (GPT-5) via Condition-agent workspace

## Central Specification Status
- `chemtools/featurizers/calculable_features.json` currently exposes **86** functional-group entries and 300+ calculable feature tokens.  
- Every entry bundles SMARTS, optional text fallbacks, labels, and category tags, and is already consumed by `chemtools/util/functional_groups.py` and `chemtools/featurizers/calculable.py`.  
- Validation helpers exist (`scripts/validate_smarts_patterns.py`, `tests/test_smarts_cache.py`) but downstream code still maintains parallel SMARTS dictionaries.

## Legacy SMARTS / Heuristic Inventories

### 1. `chemtools/router.py`
- `_SMARTS_PATTERNS` keys (26): `acid`, `acyl_halide`, `alcohol`, `aldehyde`, `alkene`, `alkoxide`, `alkyl_halide`, `alpha_beta_unsaturated`, `aryl_halide`, `borane`, `boron`, `carbonyl`, `conjugated_diene`, `cyanide`, `ester`, `grignard`, `iodide`, `ketone`, `nucleophile_n`, `nucleophile_o`, `nucleophile_s`, `organolithium`, `organozinc`, `terminal_alkyne`, `triflate`, `vinyl_halide`.  
- Text heuristics in `_rule_hits()` duplicate the same motifs plus SNAr-specific checks (heteroaryl, nitro, cyano, CF₃, sulfone, etc.).  
- Comparison vs. spec: 17 keys (e.g., `aryl_halide`, `grignard`, `organolithium`, `nucleophile_n/o/s`, `boron`) do **not** have literal name matches in the JSON and require mapping to existing tokens (`sp2_halide_*`, `boronic_acid`, `amine_primary`, etc.).

### 2. `chemtools/detection.py`
- Relies on router `_rule_hits()` output and adds inline heuristics such as SNAr EWG detection (`[N+](=O)[O-]`, `C#N`, `C(F)(F)F`, `S(=O)(=O)`), catalyst-derived flags (`catalyst_pd`, `catalyst_cu`, ...), and text scans for heteroaryl fragments.  
- These heuristics should be backed by spec tokens (e.g., `nitro`, `cyano`, `sulfonyl`, `heteroaryl_ring`) to avoid bespoke string checks.

### 3. `chemtools/featurizers/molecular.py`
- `_ELECTROPHILE_PATTERNS`: `aryl_halide`, `triflate`, `vinyl_halide`, `alkyl_halide`.  
- `_EWG_PATTERNS`: `[N+](=O)[O-]`, `C#N`, `C(F)(F)F`, `C(=O)[!O]`.  
- `_NUCLEOPHILE_PATTERNS`: `aniline`, `indole`, `phenol`, `amide`, `sec_amine`, `prim_amine`.  
- Additional textual heuristics exist for leaving-group guessing and nucleophile classification. None of these feed from the JSON spec, leading to divergent SMARTS.

### 4. `chemtools/features/role/smarts.py`
- Defines its own `_SMARTS_PATTERNS` for `amine_prim`, `amine_sec`, `amine_tert`, `aniline`, `alcohol`, and aryl halides per halogen.  
- These patterns overlap with both the router needs and the functional-group spec but are currently siloed.

### 5. `chemtools/selector_payloads.py`
- Uses module-level SMARTS constants:  
  - `_PHENOL_SMARTS = "[cX3][OX2H]"`  
  - `_FREE_ALCOHOL_SMARTS = "[CX4;!$([CX3](=O)[OX2H0-,OX1-])][OX2H]"`  
  - `_CARBOXYLIC_ACID_SMARTS = "C(=O)[OX2H1]"`  
- Alpha-amino-acid detection logic duplicates patterns already present in the JSON (`carboxylic_acid`, `free_alcohol`, `alpha_amino` equivalents).

### 6. `chemtools/util/substrate_classifier.py`
- Maintains `_SPECIAL_POSITION_SMARTS` (`benzylic`, `allylic`, `propargylic`, `alpha_to_carbonyl`) and `_ORTHO_HETEROATOM_SMARTS`.  
- Also compiles per-halogen SMARTS for alkyl halides and reuses `detect_all()` output, mixing centralized and bespoke definitions.

### 7. Miscellaneous Helpers
- `chemtools/selector_payloads`, `chemtools/featurizers/molecular`, and `chemtools/util/substrate_classifier` each call `compile_smarts()` directly on ad-hoc strings.  
- `chemtools/detection_mapper.py` and tests rely on textual heuristics (e.g., `"cross-coupling" in pred`) rather than spec-provided context tags.

## Gaps Identified
1. **Name mismatches:** Router rule names do not correspond 1:1 with spec entries; we need a mapping layer or spec aliases so `detect_all()` can return the legacy booleans expected by rule logic.  
2. **Missing spec entries:** Patterns such as `grignard`, `organozinc`, `organolithium`, and `nucleophile_*` need explicit definitions (or composite derived tokens) inside the JSON to retire custom SMARTS.  
3. **Text fallback duplication:** Router `_rule_hits()` duplicates the `text_patterns` capability already available in the spec.  
4. **Cached detection reuse:** Multiple modules parse/react to the same SMILES independently without sharing the functional-group cache, leading to repeated RDKit work.

## Next Steps (feeding into later phases)
1. Extend the JSON spec with the missing motifs or introduce alias metadata to aggregate existing granular tokens into the router’s coarser categories.  
2. Build a `FunctionalGroupCatalog` facade in `chemtools/util/functional_groups.py` that can emit router-compatible booleans and counts directly from the spec (with memoization).  
3. Refactor router/detection, CN featurizer, role SMARTS, selector payloads, and substrate classifier to consume the shared catalog instead of declaring new SMARTS.  
4. Expand the validation script/tests to ensure every SMARTS used anywhere is declared in the spec and compiled via the centralized cache.

This audit fulfills Phase 1 of the simplification plan; subsequent phases can now focus on implementation and regression coverage.
