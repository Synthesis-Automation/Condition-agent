# Chemtools Refactor Implementation Plan

## Objective

Refactor `chemtools` into a cleaner, simpler, domain-oriented package where each module can be maintained and debugged independently. Backward compatibility is explicitly not a goal. Public callers in `app/`, `llmtools/`, scripts, and tests should be migrated to the new APIs once the clean internal models are ready.

This refactor should prefer deletion and direct migration over compatibility shims. Temporary wrappers are allowed only inside an active phase when they reduce migration risk, and they must be deleted before the phase is marked complete.

## Target Architecture

```text
chemtools/
  core/              # pure primitives: contracts, errors, SMILES, RDKit helpers, SMARTS cache
  taxonomy/          # taxonomy data, validation, catalogs, motif/reaction scope services
  molecule/          # molecule-level detection and feature extraction
  reaction/          # reaction parsing, events, bond changes, CRK, reaction typing
  recommend/         # recommendation engine only
  precedent/         # precedent search only
  reagent/           # reagent taxonomy, lookup, registry
  synthesis/         # forward and retrosynthesis tools
  formatters/        # output formatting only, no chemistry inference
  visualization/     # rendering only
```

## Dependency Rules

Allowed dependency direction:

```text
core
taxonomy -> core
molecule -> core, taxonomy
reaction -> core, taxonomy, molecule
recommend -> core, taxonomy, reaction, molecule, precedent, reagent
synthesis -> core, taxonomy, molecule, reaction
formatters -> core and domain models only
app/llmtools -> explicit chemtools domain APIs
```

Rules:

- `taxonomy`, `molecule`, `reaction`, `precedent`, `reagent`, and `synthesis` must not import `recommend`.
- `formatters` must not perform chemistry inference.
- Non-taxonomy modules must not open taxonomy JSON files directly.
- All SMARTS compilation must go through `compile_smarts()`.
- Recommendation code consumes reaction/molecule features; it must not infer reaction chemistry itself.
- Reaction code owns reaction interpretation: reaction typing, CRK, events, product motifs, and bond changes.
- Do not preserve old import paths unless there is a short-lived migration reason.
- Do not keep duplicate implementations after a domain module is created.
- Do not add facade APIs that hide unclear ownership between modules.

## Current Hotspots

- `chemtools/recommend/recommender.py`: too large; mixes HTE loading, cache, taxonomy scope matching, scoring, filtering, formatting, and orchestration.
- `chemtools/featurizers/formatters/reaction.py`: too large; mixes reaction featurization, CRK construction, bond changes, product motif inference, reaction type enrichment, and formatting.
- `chemtools/featurizers/formatters/`: misleading package name because it contains domain logic.
- `chemtools/__init__.py`: heavy facade with broad imports.
- Direct `Chem.MolFromSmarts()` calls exist outside the centralized SMARTS cache.

## Progress Tracker

Status labels:

- `[ ]` Not started
- `[~]` In progress
- `[x]` Complete

## Phase 0: Architecture Lock-In

Goal: define the refactor rules before moving implementation code.

Tasks:

- [x] Confirm target package responsibilities.
- [x] Confirm allowed dependency direction.
- [x] Define canonical internal data models.
- [x] Mark legacy modules and deletion candidates.
- [x] Mark public compatibility as out of scope for the refactor.
- [x] Define which old imports will be deleted instead of wrapped.
- [x] Add this implementation plan to `docs/`.

Deliverables:

- [x] `docs/chemtools_refactor_implementation_plan.md`
- [x] Initial model/API sketch for `core`, `molecule`, `reaction`, and `recommend`

Legacy candidates:

- `chemtools/featurizers/formatters/reaction.py`
- `chemtools/recommend/recommender.py`
- `chemtools/context.py`
- `chemtools/reaction_inference.py`
- `chemtools/reaction_similarity.py`
- `chemtools/quick_reaction_glance.py`
- `chemtools/router.py`

## Phase 1: Create Clean Core

Goal: isolate stable low-level primitives and remove scattered core utilities.

Target structure:

```text
chemtools/core/
  __init__.py
  errors.py
  contracts.py
  smiles.py
  rdkit.py
  smarts.py
  reactions.py
```

Tasks:

- [x] Create `chemtools/core/`.
- [x] Move `chemtools/exceptions.py` to `chemtools/core/errors.py`.
- [x] Move `chemtools/contracts.py` to `chemtools/core/contracts.py`.
- [x] Move `chemtools/smiles.py` to `chemtools/core/smiles.py`.
- [x] Move `chemtools/util/rdkit_helpers.py` to `chemtools/core/rdkit.py`.
- [x] Move `chemtools/util/smarts_cache.py` to `chemtools/core/smarts.py`.
- [x] Replace direct `Chem.MolFromSmarts()` calls with `compile_smarts()`.
- [x] Update imports in touched modules.
- [x] Add/update tests for SMARTS cache usage.

Completion criteria:

- [x] No direct `Chem.MolFromSmarts()` calls outside `chemtools/core/smarts.py`.
- [x] Core imports no domain modules.
- [x] Tests covering changed imports pass.

## Phase 2: Make Taxonomy the Data Authority

Goal: all taxonomy JSON access and taxonomy-specific reasoning goes through `chemtools/taxonomy`.

Target structure:

```text
chemtools/taxonomy/
  __init__.py
  loader.py
  reaction_catalog.py
  compound_catalog.py
  motif_catalog.py
  scope.py
  compatibility.py
  validation.py
  models.py
  data/
```

Tasks:

- [x] Add taxonomy model objects in `taxonomy/models.py`.
- [x] Add motif catalog APIs in `taxonomy/motif_catalog.py`.
- [x] Add motif/reaction scope APIs in `taxonomy/scope.py`.
- [x] Add compatibility APIs in `taxonomy/compatibility.py`.
- [ ] Move motif parent expansion logic out of recommendation/featurizer modules.
- [ ] Move scaffold/scope map loading out of recommendation/featurizer modules.
- [ ] Move reaction type alias resolution into taxonomy.
- [ ] Move product motif constraint lookup into taxonomy.

Canonical APIs:

```python
get_reaction_type(id_or_label)
list_reaction_types()
get_motif(motif_id)
expand_motif_scope(motif_ids)
motifs_are_compatible(query, candidate)
get_allowed_product_motifs(reaction_type)
get_reagent_roles(reaction_type)
```

Completion criteria:

- [ ] Non-taxonomy modules do not open taxonomy JSON files directly.
- [ ] Recommendation code asks taxonomy APIs for motif compatibility.
- [ ] Reaction code asks taxonomy APIs for product motif constraints.

Phase 2 note: service modules now exist. Existing legacy reaction/recommendation internals still keep their current behavior until they are split and migrated in Phases 4-5.

## Phase 3: Split Molecule Logic Out of Featurizers

Goal: molecule feature extraction becomes independent from reaction feature extraction.

Target structure:

```text
chemtools/molecule/
  __init__.py
  models.py
  parse.py
  functional_groups.py
  motifs.py
  sterics.py
  electronics.py
  featurize.py
```

Tasks:

- [x] Create `chemtools/molecule/`.
- [x] Move `chemtools/util/functional_groups.py` to `molecule/functional_groups.py`.
- [x] Move molecule featurization from `chemtools/featurizers/molecule.py`.
- [x] Move motif detection/classification from `chemtools/featurizers/motifs/`.
- [x] Move aryl steric/electronic logic into `molecule/sterics_aryl.py`, `molecule/sterics_alkyl.py`, and `molecule/electronics.py`.
- [x] Define `MoleculeFeatures` and related models.
- [x] Add explicit public API in `chemtools/molecule/__init__.py`.

Canonical API:

```python
from chemtools.molecule import featurize_molecule, detect_motifs
```

Completion criteria:

- [x] Molecule code depends only on `core` and `taxonomy`.
- [x] Molecule code has no recommendation imports.
- [x] App/tests can call molecule features through `chemtools.molecule`.

## Phase 4: Create Reaction Module

Goal: reaction analysis becomes a self-contained domain, not hidden inside `featurizers/formatters`.

Target structure:

```text
chemtools/reaction/
  __init__.py
  models.py
  parse.py
  atom_mapping.py
  bond_changes.py
  events.py
  crk.py
  roles.py
  product_motifs.py
  typing.py
  feasibility.py
  featurize.py
```

Tasks:

- [x] Create `chemtools/reaction/`.
- [x] Move `chemtools/_atom_mapping.py` to `reaction/atom_mapping.py`.
- [x] Move reaction detection from `chemtools/detection.py` to `reaction/typing.py`.
- [x] Move reaction inference/validation into `reaction/inference.py` and `reaction/feasibility.py`.
- [x] Split CRK construction from `featurizers/formatters/reaction.py` into `reaction/crk.py`.
- [x] Split reaction event logic into `reaction/events.py`.
- [x] Split bond-change logic into `reaction/bond_changes.py`.
- [x] Split product motif inference into `reaction/product_motifs.py`.
- [x] Split role assignment into `reaction/roles.py`.
- [x] Define `ReactionFeatures`, `ReactionEvents`, and `CrkResult`.
- [x] Add explicit public API in `chemtools/reaction/__init__.py`.

Canonical API:

```python
from chemtools.reaction import (
    featurize_reaction,
    build_crk,
    detect_reaction_type,
    analyze_bond_changes,
)
```

Completion criteria:

- [x] Reaction code owns reaction interpretation.
- [ ] Recommendation code no longer contains reaction typing/product motif inference.
- [x] `featurizers/formatters/reaction.py` is no longer a primary implementation module.

Phase 4 note: the old reaction formatter files were moved out of `featurizers/formatters`. `reaction/featurize.py` still contains a large legacy implementation internally; `reaction/crk.py`, `reaction/product_motifs.py`, and `reaction/roles.py` currently expose domain-owned APIs around that implementation. Further internal decomposition should continue during Phases 5-6 while recommendation and formatter dependencies are removed.

## Phase 5: Rebuild Recommendation as an Orchestrator

Goal: reduce `recommend/recommender.py` into a thin coordinator.

Target structure:

```text
chemtools/recommend/
  __init__.py
  models.py
  api.py
  recommender.py
  hte_loader.py
  source_planner.py
  candidate_generation.py
  filtering.py
  scoring.py
  ranking.py
  explanation.py
  cache.py
```

Tasks:

- [ ] Move HTE CSV discovery/loading into `hte_loader.py`.
- [ ] Move HTE cache logic into `cache.py`.
- [ ] Move source group/path planning into `source_planner.py`.
- [ ] Move candidate generation into `candidate_generation.py`.
- [ ] Move condition/catalyst plausibility filters into `filtering.py`.
- [ ] Move score calculation into `scoring.py`.
- [ ] Move final ordering into `ranking.py`.
- [ ] Move text explanation/summary helpers into `explanation.py` or `formatters`.
- [ ] Keep `recommender.py` as orchestration only.

Canonical flow:

```text
RecommendationRequest
  -> reaction.featurize_reaction()
  -> taxonomy compatibility queries
  -> hte_loader.load()
  -> candidate_generation
  -> filtering
  -> scoring
  -> ranking
  -> RecommendationResult
```

Completion criteria:

- [ ] `chemtools/recommend/recommender.py` is small enough to understand as orchestration.
- [ ] Recommendation code does not parse reaction SMILES directly except through `reaction`.
- [ ] Recommendation code does not directly read taxonomy JSON.

## Phase 6: Clean Formatters

Goal: formatting becomes pure output transformation.

Target structure:

```text
chemtools/formatters/
  __init__.py
  normalization.py
  hte_output.py
  protocol.py
  reaction_output.py
  molecule_output.py
```

Tasks:

- [ ] Move all chemistry inference out of formatter modules.
- [ ] Keep only field normalization, display label cleanup, JSON/table/markdown conversion, and protocol text formatting.
- [ ] Remove RDKit usage from formatter modules.
- [ ] Update callers to use domain APIs before formatting.

Completion criteria:

- [ ] `formatters` contains no chemistry inference.
- [ ] `formatters` does not import RDKit directly.
- [ ] Formatters accept domain models or plain dictionaries as input.

## Phase 7: Reorganize Synthesis Tools

Goal: combine forward and retro tools into a clear synthesis domain.

Target structure:

```text
chemtools/synthesis/
  __init__.py
  forward/
    reactor.py
    scoring.py
    models.py
  retro/
    disconnector.py
    templates.py
    retrons.py
    registry.py
```

Tasks:

- [ ] Move `chemtools/forward/*` to `chemtools/synthesis/forward/`.
- [ ] Move `chemtools/retro/*` to `chemtools/synthesis/retro/`.
- [ ] Update imports.
- [ ] Keep synthesis dependent on `core`, `taxonomy`, `molecule`, and `reaction`.

Completion criteria:

- [ ] Synthesis code has no recommendation imports.
- [ ] Forward and retro modules share consistent model/import conventions.

## Phase 8: Simplify Public Entry Points

Goal: make imports obvious and avoid a heavy top-level facade.

Tasks:

- [ ] Slim `chemtools/__init__.py`.
- [ ] Remove or retire the broad `ChemTools` facade if no longer useful.
- [ ] Prefer explicit domain imports:

```python
from chemtools.reaction import featurize_reaction
from chemtools.recommend import recommend_conditions
from chemtools.molecule import featurize_molecule
```

Completion criteria:

- [ ] Importing `chemtools` is lightweight.
- [ ] App and LLM code use explicit domain APIs.
- [ ] Heavy optional modules are not imported at package import time.

## Phase 9: Update App and LLM Integrations

Goal: move callers onto the new clean APIs.

Update targets:

- [ ] `app/Cpd_rxn_featurization_cli.py`
- [ ] `app/Cpd_rxn_featurization_gui.py`
- [ ] `app/A_convert_to_hte_format.py`
- [ ] `app/HTE_recommender_gui.py`
- [ ] `app/services/*`
- [ ] `llmtools/recommendation_llm.py`
- [ ] `llmtools/agents.py`

Migration examples:

```python
# old
from chemtools.featurizers.unified import featurize_reaction

# new
from chemtools.reaction import featurize_reaction
```

```python
# old
from chemtools import chem

# new
from chemtools.recommend import recommend_conditions
from chemtools.precedent import search_precedents
```

Completion criteria:

- [ ] `app/` imports explicit domain APIs.
- [ ] `llmtools/` imports explicit domain APIs.
- [ ] CLI and GUI tools run against the new modules.

## Phase 10: Delete Legacy Code

Goal: remove old confusing paths after migration.

Deletion candidates:

- [ ] `chemtools/featurizers/formatters/`
- [ ] `chemtools/featurizers/unified.py`
- [ ] `chemtools/context.py`
- [ ] `chemtools/reaction_inference.py`
- [ ] `chemtools/quick_reaction_glance.py`
- [ ] `chemtools/router.py`

Completion criteria:

- [ ] No duplicate reaction/recommendation logic remains.
- [ ] No long-lived compatibility wrappers remain unless explicitly justified.
- [ ] Tests and app imports no longer reference deleted modules.

## Phase 11: Add Architecture Tests

Goal: prevent architectural drift.

Test file:

```text
tests/test_chemtools_architecture.py
```

Assertions:

- [ ] `taxonomy` does not import `recommend`.
- [ ] `molecule` does not import `recommend`.
- [ ] `reaction` does not import `recommend`.
- [ ] `synthesis` does not import `recommend`.
- [ ] `formatters` does not import RDKit directly.
- [ ] No module except `core.smarts` calls `Chem.MolFromSmarts`.
- [ ] No module outside `taxonomy` opens taxonomy JSON files directly.
- [ ] `chemtools/__init__.py` remains lightweight.

Completion criteria:

- [ ] Architecture tests pass in CI/local pytest.
- [ ] Future refactors have automated boundary checks.

## Suggested Execution Order

1. Phase 1: `core`
2. Phase 2: `taxonomy` service cleanup
3. Phase 3: `molecule`
4. Phase 4: `reaction`
5. Phase 5: `recommend`
6. Phase 6: `formatters`
7. Phase 7: `synthesis`
8. Phase 8: public entry points
9. Phase 9: app/LLM migration
10. Phase 10: legacy deletion
11. Phase 11: architecture tests

## First Implementation Slice

Recommended first concrete slice:

- [ ] Create `chemtools/core/`.
- [ ] Move SMARTS cache and RDKit helpers into core.
- [ ] Replace direct `Chem.MolFromSmarts()` calls.
- [ ] Add `tests/test_chemtools_architecture.py` with the SMARTS rule first.
- [ ] Run focused tests for RDKit helpers, forward reactor, retro template extraction, and evaluator.

This slice is low risk and establishes the first enforceable cleanup rule.
