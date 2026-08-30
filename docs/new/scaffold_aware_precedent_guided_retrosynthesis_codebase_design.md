# Scaffold-Aware Precedent Guidance for the Existing Retrosynthesis System

**Status:** Proposed post-release-gate design

**Version:** 0.1

**Date:** 2026-08-28

**Primary implementation owner:** `core_retrosynthesis`

**Supporting owners:** `reactive_taxonomy`, `condition_recommender`,
`condition_registry`, and `chem_coworker`

## 1. Purpose and authority

This document defines a scaffold-aware, precedent-guided extension that fits
the architecture and evidence rules already implemented in this repository. It
narrows the broader research proposal in
[`scaffold_aware_precedent_guided_retrosynthesis_design.md`](scaffold_aware_precedent_guided_retrosynthesis_design.md)
to a staged implementation over existing molecular, route, operator, search,
condition, and assistance contracts.

This document does not replace the current canonical references:

- [`type_agnostic_retrosynthesis_design_and_status.md`](type_agnostic_retrosynthesis_design_and_status.md)
  remains the current retrosynthesis design-and-status reference;
- [`type_agnostic_reaction_recommendation_implementation.md`](type_agnostic_reaction_recommendation_implementation.md)
  remains the primary condition-recommendation roadmap; and
- [`route_action_replay_benchmark_design_and_status.md`](route_action_replay_benchmark_design_and_status.md)
  remains the authority for observed route-action supervision and replay.

The current release gates take precedence. Scaffold-aware work may proceed as
an isolated experiment, but it must not change the canonical ranking path,
untouched evaluation set, converted corpus, or public application path until
the applicable release gates pass.

## 2. Decision summary

The system will add scaffold and route-precedent evidence as a **soft,
deterministic search prior** over the existing validated generic operator
system.

The first implementation will:

1. derive several structural-core observations from the target graph;
2. retrieve source-disjoint route-core precedents from existing local route
   artifacts;
3. group retrieved routes by graph-backed construction strategy;
4. compile the evidence into bounded validation-scheduling and ranking priors;
5. preserve a baseline branch that receives no scaffold prior;
6. generate concrete disconnections only through the existing operator ladder;
7. preserve every existing hard application and forward-signature gate; and
8. expose the evidence to the LLM only after deterministic proposal IDs exist.

The first implementation will not add a mandatory A-B reaction taxonomy, an
LLM-authored chemistry policy, free-form precursor generation, document vector
search, a graph database, or a second public retrosynthesis path.

## 3. Why this extension is useful

The present single-step system is strongest at generating, validating, ranking,
and grouping local disconnections. The route corpus adds a different kind of
evidence: which verified product sites, retained edits, synthon partitions, and
reaction sequences occurred within observed multistep routes.

Scaffold-aware guidance can use that evidence to answer bounded questions:

- Which target regions resemble the constructed region of observed routes?
- Which validated operator/site/synthon strategies recur around those regions?
- Which route motifs tend to construct rather than decorate the matched core?
- Can this evidence schedule validation across strategically distinct
  candidates more efficiently?
- Does the guidance improve `SITE1` or `STRAT1` recovery under a fixed
  validation or route-search budget?

It cannot establish that an observed strategy is optimal, transferable, or
experimentally successful on the new target. Patent occurrence is positive
route evidence, not a controlled comparison against unchosen alternatives.

## 4. Governing principles

### 4.1 Molecular and route graphs remain authoritative

Core observations, matches, route relationships, and candidate identities must
come from parsed molecular graphs, route occurrence graphs, atom lineage, graph
edits, topology, and registered operators.

Text labels, named reactions, source classifications, embeddings, and LLM
interpretations may annotate or rank compatible graph evidence. They must not
override contradictory graph evidence.

### 4.2 Type-agnostic identity remains mandatory

Primary retrieval and strategy identities use the existing generic vocabulary:

- `ReactionSignature` and `ProductTransformation`;
- reaction topology, ring changes, and local environments;
- `OP1` retained-operator identity;
- `SITE1` product-site identity;
- `SYN1` synthon-partition identity;
- `STRAT1` verified strategy identity; and
- exact realization identity for concrete precursor sets.

Electrophile/nucleophile or A-B role labels are optional interpretations. They
must not become required retrieval keys, operator identities, or hard routing
conditions.

### 4.3 Observation, interpretation, and preference stay separate

The extension has three explicit layers:

1. **Observation:** target subgraphs, attachment edges, route cores, lineage,
   ring changes, reaction signatures, and observed route actions.
2. **Interpretation:** structural-core roles, possible synthetic-core roles,
   construction-strategy families, and adaptation cautions.
3. **Preference:** validation scheduling, bounded ranking contributions,
   branch quotas, and explanations.

An interpretation failure must not discard valid graph observations. A weak or
missing precedent must return the baseline search behavior rather than reject a
candidate.

### 4.4 Hard chemistry precedes scaffold preference

Scaffold evidence cannot admit an invalid operator, rescue a failed source
round trip, bypass graph applicability, override reaction compatibility, or
turn an unverified precursor set into a result.

The order remains:

```text
operator retrieval
    -> structural applicability
    -> reverse application
    -> forward signature validation
    -> compatibility and realism gates
    -> scaffold/precedent preference
    -> diversity and presentation
```

Scaffold evidence may eventually help schedule which candidates reach an
expensive validation attempt first. Final `STRAT1` identity still requires the
existing verified graph evidence.

### 4.5 Precedent is a soft prior, not an answer key

Every scaffold-aware search must preserve an explicitly measured baseline
branch. Precedent may reorder or allocate bounded work, but it cannot consume
the entire search budget.

The baseline branch is the existing operator search under the same target,
library, chemistry gates, and total budget. It is not a weaker or different
chemistry engine.

### 4.6 Ambiguity remains typed

The system must preserve:

- multiple plausible target cores;
- symmetry-equivalent or bounded atom-lineage alternatives;
- unavailable route evidence;
- conflicting construction interpretations;
- source overlap and leakage exclusions;
- insufficient independent support; and
- unsupported adaptation dimensions.

No component may silently force one core, route family, or atom correspondence.

## 5. Existing assets to reuse

The extension is intentionally based on current repository contracts.

| Existing asset | Reuse in this design |
| --- | --- |
| `reactive_taxonomy.target_audit` | Target identity, stereochemistry, motifs, and reactive-site observations |
| `reactive_taxonomy.strategic_complexity` | Graph-derived ring construction, scaffold split, stereochemistry installation, peripheral attachment, and tactical classes |
| `ReactionSignature` and `ReactionTopology` | Generic reaction, edit, environment, and ring-change evidence |
| `ReactionRouteTree` | Occurrence-preserving observed and planned route representation |
| `RouteCoreProjection` | Route steps enriched with generic signatures, cores, minimized structures, lineage, and motifs |
| `ObservedRouteActionLabel` | Independent `SITE1`, `OP1`, `SYN1`, `STRAT1`, realization, and round-trip supervision facets |
| `CoupledRouteStrategyOccurrence` | Existing graph-backed two-step route relationship evidence |
| `GenericDisconnectionCandidate` | Concrete forward-validated one-step realization |
| `StrategyProposal` | User-facing verified strategy plus alternate realizations |
| hierarchical ranking and route-action policy | Existing deterministic scheduling and preference mechanisms |
| `RetrosynthesisConditionEvidence` | Structure-backed condition support after a concrete reaction exists |
| `chem_coworker.assistance` | Bounded evidence ledger, stable aliases, proposal-ID actions, verification, and rollback |

The first prototype should use the existing JSONL, JSON, Parquet, and local
SQLite artifacts. New external stores are not required.

## 6. Package ownership and dependency boundaries

```text
reactive_taxonomy
    structural-core observations
    target atom references
    topology and subgraph descriptors
    core-to-molecule graph matching
              |
              v
core_retrosynthesis
    route-core precedent index
    precedent matches
    strategy-family extraction
    adaptation assessment
    scaffold search priors
    one-step scheduling and multistep integration
              |
              +------------------+
              |                  |
              v                  v
condition_recommender       chem_coworker
    condition evidence          bounded LLM inspection,
    for concrete steps          selection, and explanation
```

### 6.1 `reactive_taxonomy`

This package owns molecular facts and reusable graph operations:

- structural-core candidate generation;
- target-scoped atom references;
- induced-subgraph and attachment descriptors;
- ring-system, topology, stereo-density, and complexity observations;
- exact, typed, and shape-level core keys;
- deterministic core matching; and
- validation of versioned core definitions.

It must not import `core_retrosynthesis`, `condition_recommender`, or
`chem_coworker`. It must not decide that a structural core is the preferred
synthetic strategy.

### 6.2 `core_retrosynthesis`

This package owns retrosynthetic interpretation and search behavior:

- candidate synthetic-core hypotheses over structural observations;
- route-core precedent conversion and indexing;
- route and strategy-family retrieval;
- precedent adaptation assessments;
- scaffold-aware validation scheduling;
- bounded ranking contributions;
- branch accounting and baseline reservation;
- multistep search integration; and
- scaffold-aware evaluation artifacts.

### 6.3 `condition_registry` and `condition_recommender`

These packages retain their current ownership. Scaffold or strategy labels may
be included as evidence, but condition identity and retrieval continue to use
canonical recipes, reaction signatures, edits, handles, environments,
compatibility, and provenance.

The scaffold layer must not create a condition, infer a condition role from a
route-family label, or declare a condition compatible without the existing
condition compatibility path.

### 6.4 `chem_coworker`

This package remains an application and advisory shell. It may:

- expose compact core and precedent evidence;
- let the model inspect or compare deterministic hypotheses;
- let the model select a code-generated scaffold-prior proposal ID;
- run the selected deterministic search;
- verify the retained result; and
- explain evidence and uncertainty.

It must not contain subgraph algorithms, route clustering rules, retrieval
weights, chemistry compatibility rules, or direct graph edits.

### 6.5 `app`

The application layer renders typed results and submits requests. It does not
own scaffold chemistry, retrieval, ranking, or LLM repair rules.

## 7. Core observation contracts

The exact names below are proposed. They should be implemented as immutable,
typed dataclasses with explicit `schema_version` and definition provenance.

### 7.1 Target-scoped atom reference

Do not describe atom identifiers assigned to an unmapped query molecule as
observed atom mapping. Use a target-scoped reference:

```python
@dataclass(frozen=True)
class MoleculeAtomReference:
    molecule_id: str
    component_index: int
    atom_index: int
    element: str
    formal_charge: int
    aromatic: bool
    hybridization: str
    environment_id: str
    schema_version: str
```

`atom_index` is scoped to the normalized molecule identified by `molecule_id`.
Cross-molecule retrieval uses graph-derived core keys, not equality of atom
indices.

### 7.2 Structural core observation

```python
StructuralCoreKind = Literal[
    "ring_linker_framework",
    "ring_system",
    "topological_complexity_region",
    "stereo_dense_region",
    "rare_local_subgraph",
]

@dataclass(frozen=True)
class StructuralCoreObservation:
    core_observation_id: str
    molecule_id: str
    kind: StructuralCoreKind
    atom_references: tuple[MoleculeAtomReference, ...]
    attachment_bond_keys: tuple[str, ...]
    structural_exact_key: str
    structural_typed_key: str
    structural_shape_key: str
    topology_tokens: tuple[str, ...]
    descriptor_values: tuple[tuple[str, float], ...]
    evidence_codes: tuple[str, ...]
    definition_id: str
    definition_version: str
    schema_version: str
```

Identity should bind normalized chemistry, atom membership, attachment
environment, schema version, and definition version. Display labels and model
rationales must not participate in identity.

These structural keys are distinct from the existing reaction-core
`exact_core_key`, `typed_core_key`, and `shape_core_key`. Reaction-core keys
describe a transformation; structural-core keys describe a molecule subgraph.

### 7.3 Structural candidate generation

The first definition should produce a bounded union of deterministic candidates:

1. the Bemis-Murcko scaffold backbone, which is the union of ring systems and
   their linkers, when it does not approximate the entire target;
2. explicit minimal linker paths between directly connected ring systems;
3. local interfaces around balanced non-ring graph bridges, retaining the
   central focus bond without asserting that it should be disconnected;
4. connected carbon frameworks, including branch and stereocentre evidence;
5. local stereochemical-backbone regions;
6. distinctive heterocyclic, nonaromatic, fused, spiro, or bridged ring systems,
   rather than every isolated simple ring; and
7. the largest target subgraph matching an admitted route-core key once the
   source-disjoint route index exists.

Candidates are canonicalized and deduplicated by atom membership plus attachment
environment. The first policy should cap output at five observations using
versioned deterministic tie-breaks.

Commercial availability, route counts, and LLM judgments do not belong in this
observation object.

## 8. Synthetic-core hypothesis contracts

`core_retrosynthesis` converts structural observations and route evidence into
synthetic interpretations.

```python
SyntheticCoreRole = Literal[
    "possible_constructed_core",
    "possible_inherited_scaffold",
    "possible_stereochemical_bottleneck",
    "possible_late_stage_platform",
    "unresolved",
]

@dataclass(frozen=True)
class SyntheticCoreHypothesis:
    hypothesis_id: str
    core_observation_id: str
    role: SyntheticCoreRole
    evidence_ids: tuple[str, ...]
    support_level: str
    support_score: float | None
    conflict_codes: tuple[str, ...]
    uncertainty_codes: tuple[str, ...]
    definition_id: str
    definition_version: str
    schema_version: str
```

`support_score` is a deterministic ranking score unless a separate calibration
study proves probabilistic meaning. It must not be labeled a probability by
default.

Multiple hypotheses may refer to the same structural core. A missing or
conflicting interpretation must preserve the underlying observation.

## 9. Precedent corpus and index

### 9.1 Admissible source records

The initial index is built only from records that already satisfy the relevant
route contracts:

- a valid `ReactionRouteTree`;
- a valid `RouteCoreProjection`;
- explicit source route, dataset, patent, and split provenance where available;
- per-step chemistry and evidence-quality status;
- route atom-lineage status and ambiguity; and
- independent route-action evidence facets.

An observed route may contribute at different evidence levels. Failure of
operator round-trip admission does not erase a verified product site or synthon
label, but those facets cannot be presented as an executable operator.

### 9.2 Index records

```python
@dataclass(frozen=True)
class ScaffoldPrecedentRecord:
    record_id: str
    route_tree_id: str
    route_core_id: str
    source_dataset_id: str
    source_route_id: str
    patent_id: str | None
    split: str | None
    target_molecule_id: str
    structural_exact_keys: tuple[str, ...]
    structural_typed_keys: tuple[str, ...]
    structural_shape_keys: tuple[str, ...]
    topology_tokens: tuple[str, ...]
    route_motif_ids: tuple[str, ...]
    reaction_node_ids: tuple[str, ...]
    strategy_ids: tuple[str, ...]
    evidence_quality: str
    warnings: tuple[str, ...]
    source_schema_versions: tuple[tuple[str, str], ...]
    schema_version: str
```

The primary artifact should be nested JSONL or Parquet with a manifest binding
source checksums, conversion version, definition versions, split, and counts.
A local SQLite index may accelerate exact and token retrieval. It is a derived
artifact, not the evidence source of truth.

### 9.3 Retrieval ladder

The initial retrieval ladder is:

| Level | Required evidence | Purpose |
| --- | --- | --- |
| `exact_target_route` | Exact normalized target identity, source-disjoint | Detect an existing route without conflating analogues |
| `exact_core` | Exact core graph and compatible attachment topology | Same structural core |
| `typed_core` | Element/bond-aware relaxed core plus critical-match checks | Close structural analogue |
| `shape_core` | Topology-preserving core shape plus explicit mismatches | Topological analogue |
| `route_motif` | Compatible ring-change or route-core motif | Construction-pattern analogue |
| `local_reaction` | Existing reaction-signature ladder | Validate concrete implementation |

Salts, protected forms, deprotected forms, and oxidation-state variants are not
`exact_target_route` matches. They belong in separately identified relaxed
levels with explicit transformation and mismatch evidence.

### 9.4 Query contract

```python
@dataclass(frozen=True)
class ScaffoldPrecedentQuery:
    query_id: str
    target_molecule_id: str
    core_observation_ids: tuple[str, ...]
    allowed_levels: tuple[str, ...]
    source_exclusion: "PrecedentSourceExclusion"
    maximum_records_per_level: int
    definition_id: str
    definition_version: str
    schema_version: str
```

The caller may choose bounds and source exclusions. It may not supply arbitrary
SMARTS, executable code, or unvalidated core identities.

### 9.5 Precedent match contract

```python
@dataclass(frozen=True)
class ScaffoldPrecedentMatch:
    match_id: str
    query_id: str
    record_id: str
    core_observation_id: str
    retrieval_level: str
    matching_structural_key: str
    matched_route_motif_ids: tuple[str, ...]
    matched_strategy_ids: tuple[str, ...]
    critical_matches: tuple[str, ...]
    critical_mismatches: tuple[str, ...]
    ambiguity_codes: tuple[str, ...]
    evidence_quality: str
    score_components: tuple[tuple[str, float], ...]
    score: float
    schema_version: str
```

Compatibility and critical mismatch checks run before similarity ranking.
Whole-molecule fingerprint similarity may be one score component but cannot
override a core topology, reactive-atom, stereo, or attachment mismatch.

## 10. Construction-strategy families

### 10.1 Deterministic family identity

The first strategy-family implementation clusters route segments using
graph-backed fields only:

- core atom retention and coverage;
- ring changes and reaction topology;
- ordered or partially ordered strategic-complexity classes;
- verified `SITE1`, `OP1`, `SYN1`, and `STRAT1` identities where available;
- fragment-union and ring-closure relationships;
- route-core motif identities;
- lineage-supported ordering; and
- explicit evidence-quality and ambiguity status.

Named reactions and LLM-generated descriptions may annotate a cluster after
identity construction. They do not define cluster membership.

### 10.2 Family contract

```python
@dataclass(frozen=True)
class PrecedentStrategyFamily:
    family_id: str
    identity_key: str
    structural_core_key: str
    architectural_class: str
    route_motif_keys: tuple[str, ...]
    verified_site_keys: tuple[str, ...]
    verified_operator_keys: tuple[str, ...]
    verified_synthon_keys: tuple[str, ...]
    representative_record_ids: tuple[str, ...]
    independent_source_ids: tuple[str, ...]
    total_record_count: int
    evidence_level: str
    strengths: tuple[str, ...]
    limitations: tuple[str, ...]
    ambiguity_codes: tuple[str, ...]
    definition_id: str
    definition_version: str
    schema_version: str
```

`architectural_class` should initially reuse or carefully extend the existing
graph-derived strategic classes. Any vocabulary extension belongs in a
validated, versioned definition. A free-form family name is display metadata.

### 10.3 Evidence limits

A family occurrence supports the statement that a strategy was observed. It
does not prove:

- that the route was optimized;
- that unchosen strategies failed;
- that reported yield is reproducible;
- that the strategy transfers to the query target;
- that a condition is compatible with the query; or
- that a route-action label is executable as an operator.

Independent support counts must deduplicate by an appropriate source identity,
such as patent family or document, rather than counting repeated examples as
independent evidence.

## 11. Precedent adaptation assessment

Adaptation assessment is deterministic and dimension-specific. It should avoid
one opaque transfer score.

```python
AdaptationDisposition = Literal[
    "supported",
    "caution",
    "conflict",
    "unavailable",
]

@dataclass(frozen=True)
class PrecedentAdaptationAssessment:
    assessment_id: str
    match_id: str
    disposition: AdaptationDisposition
    topology_status: str
    attachment_status: str
    reactive_environment_status: str
    stereochemistry_status: str
    functional_group_status: str
    oxidation_state_status: str
    condition_status: str
    source_independence_status: str
    evidence_ids: tuple[str, ...]
    conflict_codes: tuple[str, ...]
    caution_codes: tuple[str, ...]
    definition_id: str
    definition_version: str
    schema_version: str
```

`condition_status` consumes condition-recommender evidence only after a
concrete reaction or admitted reaction hypothesis exists. It must not be
inferred from a strategy-family label alone.

The first policy should abstain on mechanistic dependencies, scale transfer,
hazard, and unreported negative results unless typed evidence exists.

## 12. Scaffold search prior

### 12.1 Prior, not chemistry policy

The compiled object is deliberately named a prior. It affects bounded
preference and scheduling; it does not define valid chemistry.

```python
@dataclass(frozen=True)
class ScaffoldSearchPrior:
    prior_id: str
    target_molecule_id: str
    source_hypothesis_ids: tuple[str, ...]
    source_family_ids: tuple[str, ...]
    favored_site_keys: tuple[str, ...]
    favored_operator_keys: tuple[str, ...]
    favored_synthon_keys: tuple[str, ...]
    favored_topology_tokens: tuple[str, ...]
    candidate_score_adjustments: tuple[tuple[str, float], ...]
    validation_queue_quotas: tuple[tuple[str, int], ...]
    minimum_baseline_validation_slots: int
    maximum_precedent_fraction: float
    fallback_triggers: tuple[str, ...]
    evidence_ids: tuple[str, ...]
    definition_id: str
    definition_version: str
    schema_version: str
```

All fields are compiled by deterministic code from known hypotheses and
families. An LLM may select among complete prior proposal IDs; it cannot provide
site keys, operator keys, weights, or quotas.

### 12.2 Allowed first effects

The first implementation may:

- interleave provisional candidates across precedent-supported and baseline
  strategy partitions before expensive validation;
- reserve validation slots for distinct product sites and synthons;
- add a bounded ranking component after hard verification;
- reserve a fixed minimum number of baseline validation attempts; and
- expose why a candidate received precedent support or a mismatch caution.

The first implementation may not:

- remove structurally applicable baseline operators;
- bypass a score-band guard;
- assign `STRAT1` before forward verification;
- change operator compilation or source admission;
- convert route occurrence into experimental success probability; or
- use a named family as a hard filter.

### 12.3 Fallback behavior

The prior is ignored or relaxed when:

- no admissible source-disjoint precedents exist;
- every match has a critical conflict;
- core hypotheses remain unresolved;
- supported partitions yield no forward-validated candidates;
- the baseline queue is starved; or
- a configured stagnation threshold is reached.

Fallback must be reported in diagnostics. It must not silently change retrieval
or search semantics.

## 13. Single-step integration

The first product-facing experiment should be single-step because its identities,
hard validation gates, and evaluation metrics are strongest.

```text
target
  -> structural-core observations
  -> local route-core precedent retrieval
  -> deterministic strategy families
  -> zero or more scaffold-prior proposals
  -> existing L2/L1/L0 operator retrieval
  -> scaffold-aware validation scheduling plus baseline reservation
  -> existing reverse application and forward validation
  -> existing compatibility, realism, ranking, and STRAT1 grouping
  -> scaffold evidence attached to retained StrategyProposal objects
```

### 13.1 Candidate contract impact

Do not add dozens of flattened scaffold fields to
`GenericDisconnectionCandidate`. Prefer one optional typed assessment or a
parallel result projection:

```python
@dataclass(frozen=True)
class ScaffoldCandidateAssessment:
    candidate_realization_id: str
    matched_hypothesis_ids: tuple[str, ...]
    matched_family_ids: tuple[str, ...]
    adaptation_assessment_ids: tuple[str, ...]
    queue_partition: str
    scheduling_rank: int | None
    ranking_adjustment: float
    evidence_ids: tuple[str, ...]
    warnings: tuple[str, ...]
    schema_version: str
```

The structural candidate remains valid or invalid independently of this
assessment.

### 13.2 Comparison baseline

Every experiment records both paths:

- `baseline`: current search and validation scheduling;
- `scaffold_prior`: identical chemistry and total budget with the prior enabled.

The comparison must distinguish generation coverage from validation scheduling.
If the prior only finds a result because it consumed more validation attempts,
that is not an efficiency gain.

## 14. Multistep integration

Multistep integration follows only after the single-step prior passes its
ablation gate.

The observed `ReactionRouteTree` continues to represent one concrete route.
Alternative search choices remain in the planner's search-state contract and
must not be serialized as an observed route tree.

The planner may carry a compact annotation:

```python
@dataclass(frozen=True)
class RouteScaffoldState:
    core_hypothesis_id: str
    strategy_family_id: str | None
    construction_status: str
    matched_route_motif_ids: tuple[str, ...]
    adaptation_assessment_ids: tuple[str, ...]
    evidence_ids: tuple[str, ...]
    uncertainty_codes: tuple[str, ...]
    schema_version: str
```

Initially, `construction_status` is an interpretation over observed graph
changes. Protection, oxidation, and stereochemical plans remain annotations
unless registered operators can construct and forward-validate each transition.

Separate bounded queues may be tested for:

- baseline search;
- exact- or typed-core precedent adaptation;
- topology-level route motifs; and
- distinct construction-strategy families.

All queues share one explicit total expansion budget. A partial route remains
partial, regardless of precedent strength.

## 15. Condition integration

Scaffold evidence can narrow the explanation of why a precedent is relevant,
but it does not replace reaction-condition retrieval.

For each concrete validated step:

1. build or reuse its `ReactionSignature` and structural context;
2. run compatibility before similarity;
3. retrieve canonical recipes through the existing fallback ladder;
4. preserve raw and resolved condition provenance;
5. attach scaffold/route precedent only as supplemental evidence; and
6. disclose when route-level and local-reaction precedents disagree.

Condition support from a different route step or a topology-only analogue must
not be represented as direct support for the query reaction.

## 16. Bounded LLM role

### 16.1 Allowed actions

Future assistance capabilities may include:

- `inspect_core_hypothesis`;
- `compare_core_hypotheses`;
- `inspect_precedent_family`;
- `compare_precedent_families`;
- `apply_scaffold_prior`; and
- existing strategy, route, repair, and verification actions.

Each mutating action accepts only a stable proposal ID generated by code.

### 16.2 Prohibited actions

The LLM cannot provide:

- atom indices or atom maps;
- core atom membership;
- SMARTS or substructure queries;
- bond-preservation or bond-cut lists;
- reaction or precursor SMILES;
- operator, site, or synthon identities not present in evidence;
- retrieval weights or search quotas;
- conditions or protecting groups;
- hard constraints without user confirmation; or
- claims that a precedent transferred successfully.

### 16.3 Agent loop

```text
inspect target/core evidence
    -> inspect or compare precedent families
    -> select one actionable prior proposal ID
    -> deterministic search with baseline reservation
    -> compare issue, validation, and coverage diagnostics
    -> verify retained strategies or routes
    -> repeat within repair/search budgets or abstain
```

The LLM may explain why one deterministic option appears preferable. Acceptance,
rollback, cycle prevention, and budget accounting remain code-owned.

LLM assistance should be introduced only after the deterministic prior has a
measured benefit. Otherwise an LLM ablation cannot distinguish reasoning value
from a weak retrieval subsystem.

## 17. Definitions and configuration

New chemistry and ranking behavior must use validated, versioned definitions.
The likely initial files are:

```text
reactive_taxonomy/definitions/
  structural_core_observations.v1.json
  structural_core_matching.v1.json

core_retrosynthesis/definitions/
  scaffold_precedent_retrieval.v1.json
  precedent_strategy_family.v1.json
  precedent_adaptation.v1.json
  scaffold_search_prior.v1.json
```

Definitions may contain:

- descriptor names and bounded thresholds;
- allowed retrieval levels;
- critical mismatch codes;
- deterministic tie-break rules;
- minimum independent support;
- queue partitions and maximum quotas;
- score-band guards;
- fallback triggers; and
- schema and definition versions.

Executable graph behavior remains in explicit Python registries. Definitions
must not name arbitrary modules for dynamic import.

## 18. Provenance and stable identity

Every serialized result records:

- normalized target identity;
- source route, dataset, patent, document, and split identifiers when available;
- source record and route-core schema versions;
- core observation, matching, retrieval, clustering, adaptation, and prior
  definition versions;
- operator-library identity;
- search and validation budgets;
- all source exclusions;
- evidence IDs and independent-support identities;
- ambiguity, conflict, and fallback codes; and
- application or assistance lineage when a prior was selected.

Timestamps may be operational metadata. They must not participate in chemical
identity or deterministic result IDs.

Stable IDs use normalized chemistry and definition versions rather than display
labels, LLM prose, or source reaction names.

## 19. Leakage control

Precedent-guided evaluation has a higher leakage risk than ordinary ranking.
Every benchmark query must carry a source exclusion contract.

```python
@dataclass(frozen=True)
class PrecedentSourceExclusion:
    route_ids: tuple[str, ...]
    patent_ids: tuple[str, ...]
    document_ids: tuple[str, ...]
    exact_target_ids: tuple[str, ...]
    split_ids: tuple[str, ...]
    cutoff_date: str | None
    schema_version: str
```

Minimum safeguards are:

1. exclude the query route and all its reactions;
2. exclude the same patent family or document lineage;
3. bind all retrieval artifacts to the training split;
4. use scaffold-disjoint and patent-disjoint panels where applicable;
5. add forward-time evaluation when dates are trustworthy;
6. prevent the untouched answer route from appearing in model-visible evidence;
7. freeze definitions before consulting the untouched set; and
8. report exclusions and post-exclusion support counts.

An exact known route is useful for retrospective analysis but cannot count as a
held-out recovery if the same route or source lineage was retrievable.

## 20. Evaluation design

### 20.1 Single-step primary metrics

Under identical operator libraries and validation budgets, compare:

- `SITE1@k`;
- `OP1@k`;
- `SYN1@k`;
- `STRAT1@k`;
- exact realization recall;
- distinct verified strategies returned;
- candidates proposed per verified strategy;
- validation attempts per returned strategy;
- time to first verified useful strategy;
- invalid and unresolved rates;
- fallback-level distribution; and
- scaffold-prior coverage and abstention.

Exact precursor recovery remains a regression metric, not the only usefulness
metric.

### 20.2 Multistep metrics

After single-step success, measure:

- solved route found under a fixed expansion budget;
- time and expansions to first solved route;
- solved starting-material fraction for partial routes;
- strategy-family diversity;
- route-tree distance and duplicate rate;
- strongest and weakest step evidence;
- compatibility and condition coverage;
- protection and tactical-step burden;
- source-disjoint precedent support; and
- blinded chemist preference.

### 20.3 Required ablations

At minimum compare:

1. current baseline;
2. structural-core diversity only;
3. route-core retrieval only;
4. scaffold retrieval plus deterministic prior;
5. prior without baseline reservation;
6. reaction-only versus route-motif retrieval;
7. exact/typed core versus topology-level retrieval; and
8. deterministic prior versus deterministic prior plus LLM selection.

All variants use the same chemistry library, hard gates, split, and total
validation or expansion budget.

### 20.4 Review packet

The blind chemist packet should hide the source answer and include:

- target and structural-core candidates;
- retrieved precedent structures with source IDs masked consistently;
- matched and mismatched core features;
- baseline and scaffold-aware strategies in randomized order;
- validation and compatibility evidence;
- uncertainty and adaptation cautions; and
- review questions on usefulness, novelty, evidence quality, and misleading
  precedent.

## 21. Implementation sequence

### Gate 0 — Protect the current release process

Before implementation:

- record the current canonical baseline and pending release gates;
- register a separate experiment definition and result directory;
- prohibit tuning against the untouched set; and
- keep the scaffold prior disabled in public paths.

### Phase 1 — Structural-core observation prototype

Implement in `reactive_taxonomy`:

- target-scoped atom references;
- deterministic ring-system and complexity-region candidates;
- exact, typed, and shape core keys;
- attachment descriptors;
- ordering and serialization invariance tests; and
- a small chemist review artifact.

Exit criterion: candidate generation is deterministic and invariant, and a
blinded review judges at least one of the bounded observations useful on the
registered development panel.

Implementation status (2026-08-28): the deterministic prototype is implemented
in `reactive_taxonomy.structural_cores`, with versioned observation and matching
definitions, public serialization contracts, and a self-contained blind review
renderer. The v1.1 observation policy incorporates the initial drug-panel
feedback: whole-molecule complexity regions and repetitive simple-ring slots
were replaced by scaffold backbones, focused balanced bridges, explicit
ring-linker paths, continuous carbon frameworks, local stereochemical regions,
and distinctive rings. Unit coverage includes serialization and component-order
invariance, coverage and diversity limits, target-scoped atom and focus-bond
resolution, invalid-input abstention, and exact/typed/shape key behavior. It is
intentionally disconnected from route generation and ranking. The
route-core-derived candidate listed in Section 7.3 remains deferred until the
source-disjoint Phase 2 index exists. The human-review portion of the exit
criterion remains pending; passing machine tests alone does not activate the
scaffold prior.

### Phase 2 — Existing-corpus precedent index

Implement in `core_retrosynthesis`:

- conversion from route-core artifacts;
- source-exclusion contracts;
- exact, typed, shape, motif, and local-reaction lookup;
- critical mismatch checks before similarity;
- manifest and snapshot validation; and
- compact precedent review output.

Exit criterion: retrieval is reproducible, leakage checks pass, and false
neighbors are characterized before ranking integration.

### Phase 3 — Deterministic strategy families and adaptation

Implement:

- graph-backed route segment identity;
- deterministic family clustering;
- independent source support;
- typed adaptation assessments; and
- ambiguity and conflict preservation.

Exit criterion: known construction families are recovered without using source
reaction names or LLM descriptions as identity.

### Phase 4 — Shadow single-step prior

Implement:

- prior compilation;
- validation queue partitions;
- baseline reservation;
- bounded post-verification ranking evidence;
- diagnostics and explanations; and
- fixed-budget baseline comparison.

Run in shadow mode only.

Exit criterion: the prior improves a preregistered `SITE1`/`STRAT1` efficiency or
diversity metric without material regression in coverage, invalid rate, or
baseline chemistry strata.

### Phase 5 — Bounded multistep experiment

Only after Phase 4 passes:

- add `RouteScaffoldState` annotations;
- test separate bounded queues;
- preserve the existing route-tree and partial-route semantics;
- recompute condition and compatibility evidence for every retained step; and
- compare fixed-budget route outcomes.

Exit criterion: source-disjoint route solvability or expert preference improves
without increasing false solved routes or hiding unsupported steps.

### Phase 6 — Assistance integration

Only after deterministic value is demonstrated:

- add evidence projections and stable aliases;
- add inspect and compare actions;
- expose only actionable prior proposal IDs;
- reuse verification, rollback, cycle prevention, and repair budgets; and
- evaluate deterministic prior with and without LLM selection.

Exit criterion: LLM assistance improves expert-rated selection or evidence use
without inventing chemistry, reducing diversity, or increasing unsupported
claims.

### Phase 7 — Canonical integration or removal

Promote only after leakage-safe evaluation and blind review. If the prior does
not provide reproducible value, remove the experiment rather than maintaining a
parallel public path.

## 22. Testing requirements

### 22.1 Unit tests

Add positive, negative, ambiguous, and conflicting cases for:

- core candidate generation;
- atom-order and serialization invariance;
- exact, typed, and shape core identity;
- attachment mismatch;
- ring topology and stereochemistry mismatch;
- source exclusions;
- route-core conversion;
- independent support counting;
- family clustering;
- adaptation disposition;
- prior compilation and quota bounds;
- baseline reservation;
- fallback and abstention; and
- LLM proposal-ID validation.

### 22.2 Regression tests

Required chemistry strata include:

- ring construction;
- scaffold splitting;
- peripheral attachment;
- stereochemistry installation;
- functional-group interconversion;
- protection-state change;
- multi-event reactions;
- ambiguous lineage;
- incomplete reaction evidence; and
- exact route, same-core, close-core, and topology-only precedents.

### 22.3 Architecture tests

Tests must prevent:

- `reactive_taxonomy` importing retrosynthesis or application packages;
- app-layer core chemistry;
- direct `Chem.MolFromSmarts()` outside the centralized cache;
- named-family or display-label fields entering stable identity;
- LLM-supplied molecular fields in mutating actions; and
- scaffold retrieval bypassing hard validation.

### 22.4 Snapshot changes

Index and evaluation snapshots are chemistry artifacts. A change must explain:

- source and definition changes;
- coverage changes;
- retrieval-level distribution;
- false-neighbor changes;
- support-count changes;
- validation scheduling changes;
- metric changes by chemistry stratum; and
- any new rejection, ambiguity, or fallback reason.

## 23. Initial file layout

The implementation should remain small until Phase 2 demonstrates useful
retrieval.

```text
reactive_taxonomy/
  structural_cores.py
  definitions/
    structural_core_observations.v1.json
    structural_core_matching.v1.json

core_retrosynthesis/
  scaffold_precedent_models.py
  scaffold_precedent_index.py
  scaffold_precedent_retrieval.py
  scaffold_strategy_families.py
  scaffold_search_prior.py
  definitions/
    scaffold_precedent_retrieval.v1.json
    precedent_strategy_family.v1.json
    precedent_adaptation.v1.json
    scaffold_search_prior.v1.json

tests/reactive_taxonomy/
  test_structural_cores.py

tests/core_retrosynthesis_tests/
  test_scaffold_precedent_index.py
  test_scaffold_precedent_retrieval.py
  test_scaffold_strategy_families.py
  test_scaffold_search_prior.py
```

Do not create API routes, external stores, or assistance actions during the
first retrieval prototype.

## 24. Initial decision set

The initial implementation adopts these decisions:

1. Maximum five deterministic structural-core observations.
2. Retain up to three non-duplicate observations for retrieval.
3. Use query-scoped atom references, not invented atom maps.
4. Use generic graph identities as primary strategy and retrieval keys.
5. Treat named families and A-B roles as optional annotations only.
6. Build the first index from existing route-core and route-action artifacts.
7. Use local versioned artifacts and optional SQLite, not new services.
8. Run compatibility and critical mismatch checks before similarity.
9. Compile a soft prior; never let precedent define valid chemistry.
10. Reserve baseline validation capacity in every guided experiment.
11. Keep all scoring terms non-probabilistic unless calibrated.
12. Introduce no template-free generator in the MVP.
13. Add no LLM action until deterministic retrieval passes its ablation gate.
14. Require patent/source exclusion and fixed-budget evaluation.
15. Promote through one canonical path or remove the experiment.

## 25. Open questions that require evidence

The following should remain experiment questions rather than early design
assumptions:

1. Which deterministic core candidates best predict useful route precedents?
2. Does an exact/typed core key add value beyond existing reaction and route-core
   motifs?
3. How much baseline validation capacity is necessary to prevent precedent
   lock-in?
4. Is the best first integration validation scheduling, post-validation ranking,
   or multistep branch allocation?
5. Which route segment length captures strategy without fragmenting support?
6. How should independent source identity be defined across patent families?
7. Which topology relaxations produce useful analogues rather than false
   neighbors?
8. Can adaptation dimensions be evaluated reliably without text or condition
   evidence?
9. Does LLM selection improve on deterministic family scores once evidence is
   compact and typed?
10. Which metrics correlate with blinded chemist preference without leaking the
    reported route?

## 26. Final position

The codebase should become scaffold-aware by enriching the existing generic
operator search with source-disjoint route-core evidence, not by replacing it
with a new reaction taxonomy or a free-form agent.

The intended architecture is:

```text
deterministic target-core observations
    -> leakage-controlled route-core retrieval
    -> graph-backed construction-strategy families
    -> typed adaptation evidence
    -> bounded soft search prior with baseline reservation
    -> existing operator generation and hard forward validation
    -> existing compatibility, conditions, ranking, and STRAT1 grouping
    -> optional proposal-ID LLM assistance
```

This sequence preserves the repository's chemistry-first source of truth,
reuses implemented route and strategy contracts, provides a measurable
single-step MVP, and creates a controlled path toward more strategic multistep
planning without maintaining a parallel retrosynthesis system.
