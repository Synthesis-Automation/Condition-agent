# Condition Recommendation Refactor Plan

## Goal

Refactor the condition recommendation stack so it behaves like a modular execution engine with declarative tools and reusable atomic tasks, then let different optimizers or LLM policies sit on top of it.

The target is not a full rewrite of `chemtools`. The target is a cleaner public interface with:

- one deterministic reaction-context builder
- one constraints resolver
- a small set of evidence-retrieval functions
- one deterministic candidate composer
- an optional LLM review layer on top

For ChemCoworker specifically, the desired condition workflow is:

1. LLM heuristic first-pass analysis
2. deterministic reaction context / featurization
3. rule-based evidence
4. motif-based evidence
5. literature-based evidence
6. similarity-based evidence
7. optional LLM synthesis, critique, or exploration on top of the retrieved evidence

## Problem Summary

The current recommendation flow is too coarse at the tool boundary.

Observed issues:

- reaction understanding concerns are split conceptually in overlapping ways
- recommendation sources are bundled into heavy calls instead of reusable atomic evidence steps
- the LLM layer mostly receives pre-bundled results and has limited ability to plan or route work
- constraints, chemistry understanding, retrieval, and synthesis are not cleanly separated

In practice, the current "four sources together" interface behaves like a built-in optimizer instead of a reusable execution layer.

This makes the system harder to reason about and harder to extend.

## Design Principles

Keep the design simple.

Rules for the refactor:

- minimize the number of public APIs
- keep internal helpers private unless they have clear standalone reuse
- separate chemistry facts from operating constraints
- ensure all evidence sources return a common shape
- preserve backward compatibility during migration
- prefer adapters over rewrites when existing internals already work

## Proposed Simple Architecture

### Public deterministic APIs

These should become the main stable interface for recommendation orchestration.

1. `build_reaction_context(reaction, optional_metadata=None)`
2. `resolve_constraints(user_request=None, platform_context=None, inventory=None, policy=None)`
3. `get_rule_evidence(context, constraints)`
4. `get_motif_evidence(context, constraints)`
5. `get_similarity_evidence(context, constraints)`
6. `get_literature_evidence(context, constraints)`
7. `compose_condition_candidates(context, constraints, evidence_bundle)`
8. `review_or_refine_candidates(context, constraints, candidates)` optional

### Private/internal helpers

These should usually stay internal to `chemtools`:

- `_normalize_reaction(...)`
- `_extract_motifs(...)`
- `_extract_reactant_roles(...)`
- `_compute_reaction_key(...)`
- `_classify_reaction_family(...)`
- `_compute_quality_flags(...)`
- `_merge_evidence(...)`
- `_score_candidates(...)`

The main simplification is that most overlapping "reaction understanding" functions become internal substeps of `build_reaction_context(...)`.

## Recommended Data Contracts

Use small typed models or dataclasses. They do not need to be complex, but they should be stable.

### `ReactionContext`

Represents deterministic chemistry understanding for one reaction.

Suggested fields:

- `reaction_smiles_raw`
- `reaction_smiles_normalized`
- `mapping_warning`
- `reactants`
- `products`
- `reaction_family_guess`
- `reaction_family_confidence`
- `reaction_key_raw`
- `reaction_key`
- `event_kinds`
- `reacted_motifs`
- `formed_motifs`
- `spectator_motifs`
- `reacted_formed_overlap`
- `product_broad_tags`
- `product_motifs_reactive`
- `stoichiometry_delta`
- `reagent_roles`
- `quality_flags`

### `RecommendationConstraints`

Represents what the system is allowed or preferred to do.

Suggested fields:

- `scale`
- `time_budget`
- `cost_tier`
- `max_temperature`
- `air_sensitive_ok`
- `moisture_sensitive_ok`
- `allowed_solvents`
- `banned_solvents`
- `allowed_reagent_classes`
- `banned_reagents`
- `inventory_only`
- `optimizer_mode`
- `hardware_limits`
- `policy_flags`

### `EvidenceHit`

Represents one ranked piece of source-specific evidence.

Suggested fields:

- `source`
- `score`
- `confidence`
- `match_summary`
- `matched_features`
- `proposed_fragments`
- `warnings`
- `provenance`

### `ConditionCandidate`

Represents one composed recommendation.

Suggested fields:

- `rank`
- `chemicals`
- `conditions`
- `confidence`
- `source_contributions`
- `reason`
- `warnings`
- `uncertainties`

## Functional Boundaries

### 1. `build_reaction_context(...)`

Responsibility:

- canonicalize the reaction input
- extract core deterministic features
- produce one reusable reaction context artifact

This function should internally handle:

- normalization
- motif extraction
- reagent role inference
- reaction key generation
- family guess
- quality flag generation

It should not handle:

- user preference constraints
- inventory restrictions
- cost policy
- literature search
- final recommendation ranking

### 2. `resolve_constraints(...)`

Responsibility:

- combine user, platform, policy, inventory, and hardware limits into one constraint object

It should not perform chemistry featurization.

### 3. Evidence retrieval functions

Each retrieval function should be source-specific and relatively fast.

#### `get_rule_evidence(...)`

Purpose:

- retrieve rule-based condition fragments
- apply deterministic pattern logic

#### `get_motif_evidence(...)`

Purpose:

- retrieve condition fragments based on motif-level matches

#### `get_similarity_evidence(...)`

Purpose:

- retrieve precedent conditions based on reaction similarity or featurized distance

#### `get_literature_evidence(...)`

Purpose:

- retrieve or summarize literature-backed condition fragments

All four should return the same `EvidenceHit` structure so the composer does not need source-specific logic.

### 4. `compose_condition_candidates(...)`

Responsibility:

- merge evidence from all active sources
- assemble a small ranked set of full condition candidates
- score candidates deterministically
- attach confidence, warnings, and source contributions

This function should return a small number of outputs, such as:

- one primary recommendation
- two backups

### 5. Optional `review_or_refine_candidates(...)`

Responsibility:

- perform heuristic chemistry sense-checking
- help when evidence is sparse or conflicting
- improve explanation quality
- suggest backup ideas or next experiments

The LLM should sit on top of the deterministic pipeline, not replace it.

## What Should Be Simplified

The main simplification decision is:

- do not expose every chemistry subroutine as a public tool

Specifically:

- `normalize_reaction` should not remain a primary public orchestration API if it is only a substep
- `extract_motifs` should usually be internal
- `extract_reactant_roles` should usually be internal
- family classification can remain internal if only `ReactionContext` needs it

If some internal helper later proves valuable for independent reuse, it can be promoted to public API later.

## Migration Strategy

Use an adapter-based migration, not a big-bang rewrite.

### Phase 1: Define new contracts

Create a small schema module for:

- `ReactionContext`
- `RecommendationConstraints`
- `EvidenceHit`
- `ConditionCandidate`

Deliverables:

- stable type definitions
- basic serialization support

### Phase 2: Introduce `build_reaction_context(...)`

Implement one new deterministic context builder.

Initial version can wrap existing internals or existing `chemtools` functions.

Deliverables:

- new context-building API
- unit tests for representative reaction inputs
- quality flags for uncertain inputs

### Phase 3: Introduce `resolve_constraints(...)`

Add one function that merges:

- user request preferences
- inventory limitations
- hardware/platform limitations
- policy restrictions

Deliverables:

- one merged constraints object
- clear precedence rules

### Phase 4: Wrap existing recommendation engines as evidence providers

Do not rewrite everything yet.

Instead, wrap current rule, motif, similarity, and literature engines behind:

- `get_rule_evidence(...)`
- `get_motif_evidence(...)`
- `get_similarity_evidence(...)`
- `get_literature_evidence(...)`

Deliverables:

- common `EvidenceHit` output format
- adapters for existing implementations

### Phase 5: Add deterministic composition

Implement `compose_condition_candidates(...)`.

This is the main orchestration target for the agent layer.

Deliverables:

- condition assembly
- deterministic ranking
- confidence and warning generation
- source contribution summaries

### Phase 6: Update LLM orchestration

Refactor the LLM layer so it calls:

- `build_reaction_context(...)`
- `resolve_constraints(...)`
- selected evidence providers
- `compose_condition_candidates(...)`

Then optionally calls:

- `review_or_refine_candidates(...)`

Deliverables:

- simpler orchestration flow
- better compatibility with tool-calling agents
- clearer source attribution

### Phase 7: Deprecate coarse legacy entrypoints

After the new path is stable, deprecate the most monolithic recommendation entrypoints.

Deliverables:

- compatibility wrappers
- deprecation notes
- migration guide for callers

## Suggested Module Layout

Keep it compact.

Suggested package layout:

```text
chemtools/
  recommendation/
    context.py
    constraints.py
    evidence.py
    compose.py
    adapters.py
    llm_review.py
```

Suggested responsibilities:

- `context.py`: `ReactionContext`, `build_reaction_context`
- `constraints.py`: `RecommendationConstraints`, `resolve_constraints`
- `evidence.py`: `EvidenceHit` and source-specific evidence functions
- `compose.py`: `ConditionCandidate`, composition and ranking
- `adapters.py`: wrappers over existing legacy recommenders
- `llm_review.py`: optional LLM refinement and explanation

## Backward Compatibility

Backward compatibility is important.

Recommended approach:

- keep current coarse recommenders working during migration
- implement them internally via the new APIs where possible
- avoid breaking old caller payload shapes immediately
- add thin compatibility wrappers instead of duplicate logic

## Testing Plan

Focus tests on the new boundaries.

### Unit tests

- `build_reaction_context(...)` returns consistent fields
- `resolve_constraints(...)` merges sources correctly
- each evidence function returns valid `EvidenceHit` objects
- `compose_condition_candidates(...)` produces ranked candidates

### Contract tests

- all evidence providers return the same schema
- composer accepts mixed evidence bundles without source-specific branching errors

### Regression tests

- legacy recommendation entrypoints still produce acceptable outputs
- representative reactions still return recommendations after migration

## Risks

### Risk: Hidden coupling inside legacy recommenders

Mitigation:

- use adapters first
- avoid premature internal rewrites

### Risk: Too many public abstractions

Mitigation:

- keep the stable public API small
- hide most chemistry subroutines inside `build_reaction_context(...)`

### Risk: LLM layer still receives coarse inputs

Mitigation:

- pass context, constraints, and evidence artifacts separately
- let the LLM operate after deterministic retrieval, not instead of it

## Out of Scope

This plan does not include:

- full redesign of the external `chemtools` data model
- retraining similarity models
- building a complete new literature retrieval backend
- major UI changes in ChemRobox

## Recommended First Slice

The recommended first implementation slice is:

1. add the new schema models
2. implement `build_reaction_context(...)`
3. implement `resolve_constraints(...)`
4. wrap one existing source as `get_similarity_evidence(...)`
5. add a minimal `compose_condition_candidates(...)`

This gives a usable end-to-end path with minimal disruption.

## Success Criteria

The refactor is successful when:

- the public recommendation API is smaller and easier to understand
- chemistry understanding and constraints are clearly separated
- evidence sources are independently callable and reusable
- the composer can combine multiple evidence types consistently
- the LLM layer can orchestrate recommendation steps more intelligently
- legacy callers still function during migration
