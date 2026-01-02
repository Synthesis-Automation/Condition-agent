# ChemTools Interface Overview

This note captures the key modules, entry points, and CLI utilities exposed by the ChemTools project so that callers can quickly locate the right API surface. Use it as a companion to `docs/README.md` and the inline module docstrings.

## Core Facade (`chemtools.context`)

- `build_context()` / `ChemContext`: high-level factory returning an object with helper methods for reaction normalization, family detection, recommendation, precedent search, and explanation.
- Typical call flow:
  1. `ctx = build_context()`
  2. `ctx.normalize_reaction(smiles)`
  3. `ctx.recommend_conditions(request)` (delegates to rule, ML, or protocol engines depending on `RecommendRequest` flags)
  4. `ctx.find_precedents(query)`
  5. `ctx.explain(result)`

All request/response payloads mirror the dataclasses defined in `chemtools/contracts.py`, which are also reused by the FastAPI endpoints.

## FastAPI Endpoints

Implemented in `app/main.py`; each endpoint returns the same payloads as the Python API:

| Category | Endpoint(s) | Notes |
| --- | --- | --- |
| Health | `GET /health` | Readiness check |
| Reaction utilities | `POST /normalize`, `POST /detect_family`, `POST /detect_type` | SMILES normalization and family/type detection |
| Featurization | `POST /featurize/molecular`, `POST /featurize/role-aware` | Deterministic feature vectors |
| Precedents | `POST /precedent/knn`, `POST /precedent/filters`, `POST /precedent/explain` | DRFP-powered precedent search |
| Recommendations | `POST /api/v1/recommend`, `/recommend/conditions`, `/recommend/fusion` | Unified recommender output (dataset + protocol + HTE) |
| Plate design | `POST /api/v1/design_plate` | Deprecated (unavailable in unified system) |
| Core search | `POST /api/v1/core/search` | Catalyst core lookup |
| Registry | `POST /match`, `POST /registry` | SCDB matcher and reagent registry |
| Protocol | `POST /protocol/recommend` (via CLI/SDK) | DRFP similarity-based protocol suggestion |

## Recommendation Engines

- **Unified**: `chemtools.recommend.UnifiedRecommender` uses a single index spanning datasets, protocols, and HTE data.
- Legacy rule/ML fusion entry points were removed during consolidation.

## Precedent & Routing Utilities

- `chemtools.precedent.search.knn_search()` – DRFP nearest-neighbour lookup.
- `chemtools.precedent.filters.apply_filters()` – domain filters driven by `chemtools.constraints`.
- `chemtools.router.detect_family_from_reaction(smiles, use_rxn_insight=False)` – family detection with optional rxn_insight integration.

## SMARTS and Classification

- **Functional groups**: `chemtools.util.functional_groups.detect()` with RDKit or text fallbacks.
- **Reactant taxonomy** (centralised in `chemtools.reagent.types`):
  - `classify_reactant_smiles(smiles)` – SMARTS-based member identification (`ReactantMatch` dataclass).
  - `build_reactant_lookup()` – alias → canonical ID maps; `reactant_category_for()` resolves categories.
  - `normalize_reaction_type(label)` / `required_reactant_categories(reaction_id)` – reaction taxonomy helpers.

- **Protocol SMARTS tooling** (`chemtools.protocol.smarts_generator_cli`):
  - `SmartsGenerator(reaction).build_applicability()` → `ReactionSmartsApplicability` (core + guard SMARTS).
  - Visualization helpers: `visualize_smarts_pattern`, `visualize_reaction_smarts`, `visualize_pattern_with_examples` (RDKit-backed).
  - Batch updates: `chemtools.protocol.batch_update_protocol_smarts.ProtocolSmartsUpdater` with `process_protocol_file()` / `process_all_protocols()`.

## Reagent & Taxonomy Management (`chemtools.reagent`)

- Runtime lookup: `find_reagent`, `enrich_reagent_info`, `load_reagent_database`, `filter_precedents_by_database_availability`.
- Curation utilities: `TaxonomyStore`, `RoleHeuristics`, `taxonomy_cli` (invoked via `python -m chemtools.reagent.taxonomy_cli`).
- Analytics: `get_database_statistics`, `print_database_summary`, `find_reagents_by_family`, etc.

## CLI Entrypoints

- `python -m chemtools.rule.cli` – Scheme Condition DB search (JSONL/CSV support).
- `python -m chemtools.protocol.cli` – protocol index build/stats/show-family/show-tag.
- `python -m chemtools.protocol.smarts_generator_cli --reaction ...` – generate SMARTS patterns (`--batch`, `--check-rdkit`).
- `python -m chemtools.reagent.taxonomy_cli` – taxonomy maintenance.
- `python app/local_recommendation_cli.py` or `python scripts/local_recommendation_cli.py` – interactive recommendation runner (rule, ML, protocol, LLM).
- `python scripts/ui_gradio.py` – browser demo.
- `python data-processor/process_reactions.py --rdf ... --txt ...` – SciFinder TXT/RDF merge into CSV.
- HTE helpers in `data-processor/HTE_data/`: e.g. `test_recommender.py`, `quick_test.py`, `simple_condition_recommender.py`.

## Key Supporting Modules

- `chemtools/smiles.py` – SMILES parsing & normalization utilities.
- `chemtools/util/rdkit_helpers.py` – RDKit availability checks, molecule parsing, canonicalization.
- `chemtools/util/substrate_classifier.py` – legacy heuristics (prefer `chemtools.reagent.types`).
- `chemtools/explain.py` – explanation payload assembly.
- `chemtools/output_formatter.py` – consistent JSON packing for rule/ML/protocol outputs.
- `tests/` – pytest suite illustrating the public API usage (e.g., `tests/test_step*_*.py`, `tests/test_reagent_types.py`, `tests/test_smarts_generator.py`).

## Usage Notes

- Always prefer the structured contracts in `chemtools/contracts.py` when constructing request objects for rule/ML/protocol recommenders.
- RDKit is optional; set `CHEMTOOLS_DISABLE_RDKIT=1` to force text fallbacks (some SMARTS-intensive features will degrade gracefully).
- For reproducibility, call the CLI wrappers or FastAPI endpoints described above; they mirror the in-process functions and use identical JSON schemas.
- The `data-processor/HTE_data` directory keeps auxiliary scripts and markdown reports that rely on the canonical JSON definitions under `chemtools/reagent/data/`.

Refer to this file whenever you need a quick reminder of where a particular ChemTools interface lives or which command-line entry point to run.
