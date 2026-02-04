# Reaction Key and Reaction Type Determination Workflow

This document summarizes how the system determines the Reaction_Key and the
reaction type from a reaction SMILES string. It reflects the current code paths
in `chemtools/featurizers/formatters/reaction.py` and `chemtools/detection.py`,
plus taxonomy data in `chemtools/taxonomy/data/`.

## Inputs and Taxonomy Sources

- Input: reaction SMILES with `>>` separator.
- Core taxonomy:
  - `chemtools/taxonomy/data/reaction_types.v4.0.json`
  - `chemtools/taxonomy/data/compound_logic.json` (motif sets)
- Motif detection and aggregation:
  - `chemtools/featurizers/formatters/reaction.py`
  - `chemtools/featurizers/aggregation.py`
  - `chemtools/featurizers/motif_*` helpers (via `build_molecule_bundle`)

## Reaction Key Determination (featurize_reaction)

Entry point: `chemtools/featurizers/formatters/reaction.py:featurize_reaction()`.

1. Normalize reaction SMILES via `chemtools/smiles.normalize_reaction`.
2. Featurize reactants and products with `build_molecule_bundle` to detect motif
   IDs and context motifs.
3. Aggregate reaction features with `aggregate_reaction_features`:
   - Computes motif deltas (reacted, formed, spectator) by comparing reactant
     and product motif counts.
4. Compute bond-change analysis (atom-mapped) when RDKit is available.
5. Infer broad product tags from bond changes, with a SMARTS fallback for
   aryl-sulfur bonds (e.g., `Product_Aryl_S`).
6. Filter reacted motifs to those touching bond changes; if filtering yields
   an empty set, fall back to all reacted motifs.
7. Build the single Condition Recommendation Key (CRK):
   - `|Reactants -> Product | bond_formed: ... | bond_broken: ... | spectators: ...`
   - Product is optional; when disabled it is reported as `[]`.

Outputs stored in the reaction bundle:

- `reaction_key` (CRK)
- `product_broad_tags`
- `aggregates.reacted_motifs`, `aggregates.formed_motifs`,
  `aggregates.spectator_motifs`

## Reaction Type Determination

### Primary detection API (recommended)

Entry point: `chemtools/detection.py:detect_reaction_type()`.

1. Validate reaction SMILES format (`>>` required).
2. Extract the reaction key and motif deltas via
   `extract_reaction_key()` which calls `featurize_reaction()` (roles disabled).
3. Load reaction type definitions and motif sets:
   - `reaction_types.v4.0.json`
   - `compound_logic.json`
4. For each reaction type:
   - Expand motif sets in reactant and product slots.
   - Enforce `exclude_reacted` patterns.
   - Require each reactant slot to match at least one reacted motif.
   - Require at least one product slot to match a formed motif (if defined).
5. Score and rank matches:
   - Confidence = matched slots / total slots.
   - Boost if both electrophile and nucleophile match.
6. Return `DetectionResult`:
   - Ranked matches, top match, reacted/formed motifs, and `reaction_key`.

### Lightweight validation inside featurize_reaction

`featurize_reaction()` does not run the full detection pipeline. It initializes
`reaction_type` as `Unknown` and then applies a taxonomy-driven validation step
using `validate_detection_with_reacted_motifs()`:

- Uses `chemtools/taxonomy/reaction_catalog.py` to load slot definitions and
  motif sets.
- If reacted/formed motifs satisfy a reaction type definition, it can override
  the initial `Unknown` and attach validation metadata to
  `reaction.detection.validation`.

This is a correction/validation step, not the main detection API.

## Outputs Summary

- Reaction key and motif deltas are produced by the featurizer:
  - `reaction_key` (CRK), `product_broad_tags`, `aggregates.*`
- Reaction type detection (primary) is produced by `chemtools/detection.py`:
  - `DetectionResult.matches`, `DetectionResult.top_match`,
    `DetectionResult.reaction_key`
- Featurizer output may include a validated reaction type in the core bundle,
  but the canonical detection workflow is `chemtools/detection.detect_reaction_type`.

## Improvement Plan (Proposed)

### 1) Consolidate reaction type detection paths

Goal: make `chemtools/detection.detect_reaction_type()` the single source of
reaction type truth, and keep featurizer output aligned with it.

- Replace the "Unknown + validation-only" placeholder in
  `featurize_reaction()` with an optional call to
  `chemtools.detection.detect_reaction_type()` when
  `options.get("include_reaction_type")` is true.
- Store the result in `reaction["reaction_type"]` using the
  `format_reaction_type_summary()` shape to keep schema stable.
- Keep the current lightweight validation for cases where detection is skipped,
  but mark it explicitly as `validation_only` in `reaction.detection.validation`.

### 2) Tighten reaction key generation quality

Goal: reduce ambiguous or misleading Reaction_Key strings.

- Introduce a reaction-key quality flag (e.g., `reaction_key_quality`) based on:
  - presence of product SMILES
  - RDKit mapping availability
  - non-empty reacted/formed motifs
- When quality is low, emit a warning in `result.meta.errors`.

### 3) Strengthen taxonomy-driven constraints (general)

Goal: improve chemical specificity without introducing reaction-specific
patches or ad-hoc heuristics.

- Keep constraints driven by taxonomy definitions and motif-slot logic only.
- Prefer constraints that are chemically general (slot requirements,
  reactant/product balance, and allowed/forbidden motif sets).
- Maintain a small, audited motif-set surface in `compound_logic.json` to
  reduce false positives without hardcoding exceptions for particular families.

### 4) Add regression tests for key/type determinism

Goal: prevent drift and regressions across refactors.

- Add tests for:
  - CRK-v1 format and stability across representative reactions.
  - Detection consistency between `detect_reaction_type()` and CRK-derived
    signatures for the same reaction SMILES.
- Prefer fixture-driven tests in `tests/test_reaction_key_generation.py` and
  `tests/test_reaction_type_detection.py`.

### 5) Metrics and coverage tracking

Goal: track unknowns and taxonomy gaps systematically.

- Add a small CLI/diagnostic report:
  - percent reactions with `Unknown` type
  - percent reactions with empty `reaction_key`
  - top motif pairs in reacted/formed sets that yield no match
- Use the report to prioritize taxonomy expansion.
