Legacy note: reactant_types.json has been removed; use organic_compounds.v1.3.json for reactant type definitions.

Taxonomy Centralization â Phases 0 & 1
======================================

Overview
--------
- Objective: consolidate reaction, reactant, and reagent taxonomies into a single, canonical registry that downstream tools (router, recommenders, CLIs, data processors, LLM agents) can trust.
- Scope: finish Phase 0 (recon) and Phase 1 (schema design) deliverables so implementation work can begin with clear requirements, ownership, and compatibility guarantees.

Phase 0 â Current-State Recon
-----------------------------

### Source Assets
- `chemtools/reagent/data/reactant_types.json`
  - JSON keyed by reactant macro-types (`ArX*`, `RNH2`, `BoronicAcid`), each with members supplying `id`, `name`, `smarts`.
  - No explicit category hierarchy beyond implicit naming; contains aliases in free text but no canonical slug format.
- `chemtools/reagent/data/reaction_types.json`
  - Nested dict grouped by coarse category (`C-C_bond_formation`, `C-N_bond_formation`, â¦) with `reactions` list carrying fields (`id`, `name`, `aliases`, `description`, `reactants`, `typical_catalysts`, `typical_conditions`).
  - Reactant references are arrays of reactant tokens (e.g., `["ArX*"]`, `["ArNH2"]`); no enforcement that tokens exist in `reactant_types.json`.
- `data/reagent_db/reagents.csv`
  - Flattened registry for all roles with `role_1`/`role_2` columns; stores concrete reagents with `name`, `abbreviation`, identifiers, and role metadata.
  - Role-specific metadata is stored in the shared CSV rather than per-role JSON files.
- `chemtools/taxonomy/data/reagent_roles.v2.json`
  - Unified roles + families (single source of truth) used for classification and validation.
- `chemtools/reagent/constants.py`
  - Default family fallbacks, role priority order, and registry file mapping for reagent entries.
- `chemtools/reagent/taxonomy_store.py`, `taxonomy_utils.py`
  - Helpers for loading, caching, and searching reagent families (independent schema from reaction/reactant json).
- `data-processor/HTE_data/*`
  - Validation utilities (`validate_reaction_types.py`, `analyze_reactant_mapping.py`, `fix_reactants_auto.py`) that assume the current JSON layout and reactant naming conventions.
- `chemtools/recommend/utils.py`
  - `_FAMILY_ALIASES` dict normalizing legacy family names (e.g., `Buchwald_CN` -> `C_N_Coupling`).
- `chemtools/router.py`
  - SMARTS-based heuristics tie detected reactant patterns and metal hints to family strings; no direct link to taxonomy json.

### Consumer Inventory
- **Core APIs**: `chemtools.reagent.types`, `chemtools.reagent.lookup`, `chemtools.reagent.analytics`, `chemtools.reagent.validator` load taxonomy files directly.
- **Recommendation stack**: `chemtools/recommend/modules/recommender.py`, `structured.py`, `utils.py`, `precedent` package rely on family strings and reactant tokens defined in the JSON.
- **Router/detection**: `chemtools/router.py` and `chemtools/context.ChemTools.detect_family` reference alias strings and manually detect reactant categories (e.g., aryl halide).
- **Protocol tooling**: `chemtools/protocol/*`, `show_protocol_format.py` expect family IDs aligned with precedent datasets.
- **CLI/UI**: `app/local_recommendation_cli.py`, `app/cross_family_recommendation_cli.py`, `app/web_recommendation_cli.py`, plus supporting docs, accept family override names that mirror the current taxonomy.
- **LLM integrations**: `llmtools/agents.py`, `llmtools/prompts.py`, `llmtools/README.md` embed current family tokens (`Buchwald_CN`, `Suzuki`) in examples.
- **Data-processing scripts**: `data-processor/HTE_data/*.py`, `generate_dataset_summary.py`, `show_protocol_format.py`, etc., parse or generate fields referencing the existing JSON structure.
- **Tests**: fixtures and assertions across `tests/test_precedents.py`, `tests/test_cross_family_search.py`, `tests/test_protocol_smarts.py`, `tests/test_batch_update_protocol_smarts.py`, and others lock in current IDs and reactant tokens.

### Compatibility & Pain Points
- Alias sprawl: identical chemistry concepts appear as `Buchwald_CN`, `C_N_Coupling_Pd`, `Pd_Buchwald`, `Buchwald-Hartwig` across modules; alignment is manual.
- Reactant tokens (`ArX*`, `RNH2`, etc.) are not typed; consumers assume naming conventions that differ between JSONs and data sources.
- Reagent role/family taxonomy is siloed from reaction/reactant taxonomies, so we cannot programmatically assert that a `reaction_type` demands particular reagent families.
- Validation scripts are decentralized and cannot guarantee referential integrity (lack of shared registry, versioning, or schema enforcement).
- Migration risk: external users and LLM prompts depend on legacy names; we need alias resolution and telemetry to avoid breaking CLI/API behavior.

Phase 1 â Unified Schema Design
-------------------------------

### Canonical Entities
1. **ReactionCategory**
   - `id` (snake_case, e.g., `cross_coupling`)
   - `name`, `description`
   - Relationship: contains multiple `ReactionType` entries.
2. **ReactionType**
   - `id` (snake_case, e.g., `suzuki_coupling`)
   - `category_id`, `display_name`, `description`
   - `aliases` (list of strings/regex), `source_ids` (external dataset IDs)
   - `required_reactants` (list of `ReactantTypeRequirement`) with stoichiometry hints.
   - `required_roles` (list of reagent-role expectations with optional default families).
   - `metadata` (typical catalysts/conditions, precedence dataset, protocol IDs).
3. **ReactantType**
   - `id` (snake_case, e.g., `aryl_halide`)
   - `display_name`, `description`
   - `smarts` patterns, `members` (optional list with `id`, `smarts`, `aliases`)
   - `category` (`aryl_halide`, `boron_partner`, etc.), optional `functional_group_tags`.
4. **ReagentRole**
   - `id` (`ligand`, `metal_catalyst`, `base`, etc.)
   - `display_name`, `priority`, `default_family_id`
   - Relationship: ties to `ReagentFamily` definitions.
5. **ReagentFamily**
   - `id`, `role_id`, `display_name`, `properties` (role-specific metadata schema)
   - Links to reagent registry entries (`data/reagent_db`) and supports validation.
6. **Alias**
   - Global table mapping legacy tokens (`Buchwald_CN`, `ArX*`, `Pd_Buchwald`) to canonical `entity_type` + `entity_id`.
   - Supports case-insensitive lookup, optional regex patterns for datasets.

### Storage Layout
```
chemtools/
  taxonomy/
    __init__.py           # Registry loader + public API
    models.py             # Dataclasses/Pydantic models
    registry.py           # Cache, lookup, validation logic
    data/
      manifest.json       # Top-level metadata (version, generated_at, schema_version)
      reaction_categories.json
      reaction_types.json
      reactant_types.json
      reagent_roles.json
      reagent_families.json
      aliases.json
```
- `manifest.json` records schema version, build hash, and compatibility flags.
- Each JSON file uses canonical IDs; `reaction_types.json` references `reactant_types` and `reagent_roles` by ID.
- Optionally provide a single `taxonomy.jsonl` for downstream export convenience.

### Registry API (Concept)
```python
from chemtools.taxonomy import registry

reg = registry.load()                      # memoized singleton
rxn = reg.get_reaction_type("suzuki_coupling")
reg.resolve_alias("Buchwald_CN")           # -> ReactionType(id="cn_coupling_palladium")
reg.iter_reactant_types(category="aryl")   # filter helpers
reg.required_reagent_roles("suzuki_coupling")  # structured requirements
```
- Registry enforces referential integrity at load time and raises descriptive errors for missing links.
- Provides typed results (dataclasses) so IDEs and static tools can reason about fields.

### Validation & Versioning
- Schema guardrail: JSON Schema stored alongside models; validation happens during build and CI (`python -m chemtools.taxonomy.validate`).
- Include semantic versioning for taxonomy (`taxonomy_version: "1.0.0"`); bump when breaking changes occur.
- Provide diff tooling to compare old vs new taxonomy (e.g., `compare_taxonomies.py`).

### Migration Strategy Outline
- Author conversion scripts under `scripts/taxonomy/` to map current JSONs into the unified schema (`convert_legacy_taxonomies.py`).
- Maintain legacy exports for one release cycle; shims look up canonical IDs via the alias table.
- Capture telemetry or logging when a legacy alias is used to highlight remaining consumers.

### Documentation Deliverables
- Developer guide covering new structure, migration helpers, and API usage.
- Update CLI/README examples to reference canonical IDs but mention alias support for compatibility.
- Provide data dictionary for each entity to unblock cross-team adoption (protocol, ML, LLM).

Open Questions for Phase 2+
---------------------------
- Where to host the master reagent registry (current `data/reagent_db` is large); consider splitting into normalized tables vs status quo.
- Decide whether to generate taxonomy files via build step (e.g., from spreadsheets) or commit as source-of-truth JSON.
- Determine ownership of alias governance (who approves new aliases, how to track provenance).
- Establish testing fixtures for taxonomy-driven workflows (mock registry vs real data slices).

Next Actions Snapshot
---------------------
- Ran `python scripts/taxonomy/convert_legacy_taxonomies.py --taxonomy-version 0.1.0` to seed the unified dataset.
- Current registry footprint (`python -m chemtools.taxonomy.validate`):
  - 10 reaction categories
  - 48 reaction types
  - 29 reactant types
  - 13 reagent roles
  - 154 reagent families
  - 451 alias entries
- Guardrails in place for further work:
  - `chemtools/taxonomy/registry.py` enforces referential integrity.
  - `chemtools/taxonomy/validate.py` provides a lightweight CLI for sanity checks.
  - `tests/test_taxonomy_registry.py` ensures both registry loading and the conversion dry-run stay green.

With these Phase 0â1 outputs finalized, implementation (Phase 2 onward) can proceed with clear specifications and compatibility constraints.
