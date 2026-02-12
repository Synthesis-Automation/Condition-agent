# Synthon Integration Progress

Last updated: 2026-02-12

## Scope
Track phased migration to taxonomy-driven synthon logic for reaction detection and recommendation.

## Phase Status

- Phase 1: Completed
  - Added synthon taxonomy (`chemtools/taxonomy/data/synthons.v1.json`).
  - Added synthon loader support (`chemtools/taxonomy/loader.py`).
  - Implemented synthon classifier/selector (`chemtools/synthon.py`).
  - Replaced heuristic electrophile/nucleophile pickers in:
    - `chemtools/recommend/utils.py`
    - `chemtools/precedent/loader.py`
  - Added tests: `tests/test_synthon_role_assignment.py`.

- Phase 2: Completed
  - Integrated synthon evidence into reaction-type matching (`chemtools/detection.py`):
    - Reactant-level synthon extraction from reaction SMILES.
    - Slot matching fallback to synthon motifs when reacted motifs are missing.
    - Slot source tracking (`motif` vs `synthon`).
    - Detection payload now includes `synthon_evidence`.
  - Propagated synthon slot metadata through compatibility API:
    - `chemtools/featurizers/detection/__init__.py`
  - Exposed synthon slot evidence in primary detection payload formatting:
    - `chemtools/featurizers/formatters/reaction.py`
  - Added tests: `tests/test_synthon_detection_phase2.py`.

- Phase 3: Completed
  - Added synthon schema support to reaction catalog definitions:
    - `ReactionTypeDefinition.synthons` in `chemtools/taxonomy/reaction_catalog.py`
    - `normalize_reaction_synthons(...)` + synthon slot expansion
  - Added synthon set loading for taxonomy-driven slot expansion:
    - `SYNTHON_FILE` + `_load_synthon_sets()` in `chemtools/taxonomy/reaction_catalog.py`
  - Extended reaction taxonomy validator:
    - Validates `synthons` shape and allowed slot keys
    - Flags unknown synthon IDs
    - Normalizes present `synthons` payloads
    - File: `chemtools/taxonomy/validate_reaction_types.py`
  - Added tests:
    - `tests/test_reaction_catalog_synthons.py`
    - `tests/test_validate_reaction_types_synthons.py`

## Validation Log

- Ran: `pytest -q tests/test_synthon_role_assignment.py`
- Ran: `pytest -q tests/test_synthon_detection_phase2.py`
- Ran: `pytest -q tests/test_hte_recommender.py -k "precedent_recommendations_respect_source_group_filter or recommend_returns_precedent_when_structured_match_missing or recommend_returns_precedent_when_structured_match_missing_for_protocols"`
- Ran: `pytest -q tests/test_reaction_catalog_synthons.py tests/test_validate_reaction_types_synthons.py`
- Ran: `pytest -q tests/test_reaction_type_detection_mapping.py tests/test_detection_validation_evidence.py`

## Remaining Roadmap

- Phase 4: Add synthon features to precedent/recommender scoring/ranking.
- Phase 5: Coverage-driven synthon expansion workflow from unresolved clusters.
- Phase 6: Remove remaining legacy heuristic paths and finalize hardening.
