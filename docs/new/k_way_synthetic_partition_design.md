---
title: "General Synthetic Partition Landscape and Strategic Route Realization"
version: "0.2"
status: "Design and implementation plan"
encoding: "UTF-8"
---

# General Synthetic Partition Landscape and Strategic Route Realization

## 1. Purpose

The first purpose of this system is to provide a broad, inspectable picture of
how a target molecule may be divided into synthetic modules. The second purpose
is to test whether selected divisions can be realized as chemically credible
routes.

These are separate tasks:

```text
Target molecular graph
        ↓
Synthetic partition landscape
        ↓
Strategic partition selection
        ↓
Route realization
        ↓
Precedent, condition, and feasibility assessment
```

The landscape is an analysis of the target. A route realization is an
evidence-backed sequence of transformations. A useful partition may have many
realizations, and failure to realize a partition must not erase the underlying
molecular analysis.

The molecular graph remains the source of truth. Reaction names, module names,
and LLM interpretations are optional annotations and must not override graph,
atom-correspondence, or transformation evidence.

## 2. Main decisions

### 2.1 Any positive module count is valid

For target atom set `A(T)`, a synthetic partition is:

\[
P(T)=\{M_1, M_2, \ldots, M_k\}, \qquad k \ge 1
\]

with:

\[
\bigcup_i M_i=A(T), \qquad M_i\cap M_j=\varnothing \text{ for } i\ne j
\]

Interpretation of `k`:

| Value | Meaning |
| --- | --- |
| `k = 1` | Preserve the target-derived framework, use unary transformations, or report that no useful fragmentation was found. |
| `k = 2` | Conventional two-part disconnection. |
| `k >= 3` | Multi-module or convergent construction view. |

The system must not manufacture extra modules to reach a preferred value of
`k`. Returning `k = 1`, `k = 2`, or an explicit abstention is a valid result.

### 2.2 Modules are symmetric

The base partition has no required core, appendage, linker, or donor. Module
order has no chemical meaning, and partition identity must be invariant to
module permutation.

Modules may still differ in measured properties such as:

- atom count and connectivity;
- graph complexity;
- stereochemical burden;
- number and type of attachment interfaces;
- precedent and availability evidence;
- expected synthetic difficulty.

Terms such as `core`, `appendage`, and `donor` may be added as optional,
evidence-backed interpretations. They are not routing keys and do not
participate in partition identity.

### 2.3 A target module is not a concrete precursor

For:

```text
Target = A + B + C
```

the expression means that the target-derived atoms are partitioned into
`A`, `B`, and `C`. It is not a proposed three-component reaction.

Each target module can have many latent precursor states:

```text
A ← A*halide, A*boronate, A*carbonyl, A*protected
B ← B*1, B*2, ...
C ← C*1, C*2, ...
```

A route may then be:

```text
B* + C* → BC*
BC* + A* → protected target
protected target → Target
```

or:

```text
A* + B* → AB*
AB* + C* → Target
```

One partition therefore has zero, one, or many route realizations.

### 2.4 Analysis, interpretation, and realization remain separate

```text
Observation
  target graph, atom properties, functional groups, stereochemistry,
  candidate edit sites, retrons, and local environments

Interpretation
  candidate partitions, optional module roles, strategic interfaces,
  latent-state hypotheses, evidence, and uncertainty

Realization
  explicit precursor structures, validated one-step operators,
  dependencies, route order, precedents, and conditions
```

An interpretation failure must not delete valid observations. An unrealized
partition must be labeled as such rather than presented as a route.

## 3. Scope

The first version should:

1. accept a normalized, stereochemical target structure;
2. generate several role-neutral partitions across useful values of `k`;
3. retain exact target-atom membership for every module;
4. show the strategic interfaces associated with each partition;
5. project observed or planned route frontiers into the same partition model;
6. attempt bounded realization with already validated one-step operators;
7. preserve latent handles and non-target atoms separately;
8. report realization status, evidence, warnings, and unresolved gaps;
9. return strategically diverse results, including valid `k = 1` and `k = 2`
   outcomes;
10. support human review through JSON and HTML/SVG output.

The first version will not claim to:

- enumerate every chemically possible synthesis;
- prove experimental feasibility from the molecular graph alone;
- guarantee commercial availability, yield, scale, or safety;
- infer mechanism from overall reaction records;
- turn several physical reactions into one composite reaction;
- allow an LLM to create unvalidated atoms, maps, or graph edits.

The landscape is bounded by the available chemistry definitions, operator
library, route corpus, and search limits. Those coverage limits must be visible
in the result.

## 4. Conceptual model

### 4.1 Target module

A `TargetModule` is an immutable set of target atom references. It contains only
target-derived atoms.

```python
@dataclass(frozen=True)
class TargetModule:
    module_id: str
    target_atom_maps: frozenset[int]
    attachment_atom_maps: frozenset[int]
    graph_complexity: float | None
```

`module_id` is content-derived. A display index such as `M1` is assigned only
after canonical sorting.

### 4.2 Optional module annotation

```python
@dataclass(frozen=True)
class ModuleAnnotation:
    module_id: str
    label: str
    proposed_role: str
    confidence: float | None
    evidence: tuple[str, ...]
    warnings: tuple[str, ...]
```

Examples of optional roles include `complex_fragment`, `linker`,
`functional_carbon_source`, and `heteroatom_source`. `unresolved` is always
valid. These annotations must not affect module identity or hard validation.

### 4.3 Latent module state

A latent state describes how a module could appear in a precursor.

```python
@dataclass(frozen=True)
class LatentModuleState:
    latent_state_id: str
    module_ids: tuple[str, ...]
    mapped_smiles: str
    target_atom_maps: frozenset[int]
    non_target_atoms: tuple["NonTargetAtom", ...]
    state_annotations: tuple[str, ...]
    evidence_status: str
```

Non-target atoms are classified explicitly:

```text
TACTICAL_HANDLE
PROTECTING_GROUP
AUXILIARY
COUNTERION_OR_SALT
UNKNOWN
```

A latent state may represent one module or an already joined set such as
`AB*` or `BC*`.

### 4.4 Strategic interface

A strategic interface connects modules whose atoms are related by one or more
target-side edits or by an explicitly defined transformation retron.

```python
@dataclass(frozen=True)
class StrategicInterface:
    interface_id: str
    module_ids: tuple[str, ...]
    target_atom_maps: frozenset[int]
    target_bond_edits: tuple["BondEdit", ...]
    candidate_operator_ids: tuple[str, ...]
    interpretation_annotations: tuple[str, ...]
    evidence_status: str
```

Most interfaces connect two modules, but the contract permits a hyperedge for a
validated multicenter or multicomponent transformation. A hyperedge must not be
used to collapse a multistep sequence.

### 4.5 Synthetic partition

```python
@dataclass(frozen=True)
class SyntheticPartition:
    partition_id: str
    target_id: str
    modules: tuple[TargetModule, ...]
    interfaces: tuple[StrategicInterface, ...]
    annotations: tuple[ModuleAnnotation, ...]
    generation_evidence: tuple[str, ...]
    warnings: tuple[str, ...]
    schema_version: str
    definition_version: str
```

Partition identity uses:

```text
normalized target identity
+ canonical unordered collection of target-atom module sets
+ normalized target-side interface edits
+ schema and definition versions
```

It does not use module order, display labels, optional roles, source reaction
names, or LLM prose.

### 4.6 Partition landscape

```python
@dataclass(frozen=True)
class SyntheticPartitionLandscape:
    target_id: str
    target_smiles: str
    partitions: tuple[SyntheticPartition, ...]
    searched_k_values: tuple[int, ...]
    generation_coverage: tuple[str, ...]
    unresolved_motifs: tuple[str, ...]
    diagnostics: "PartitionSearchDiagnostics"
    schema_version: str
```

The landscape is the main target-analysis result. It may contain competing and
partially supported views.

### 4.7 Strategic route realization

```python
@dataclass(frozen=True)
class StrategicRouteRealization:
    realization_id: str
    partition_id: str
    route_tree_id: str | None
    frontier_states: tuple[LatentModuleState, ...]
    interface_realizations: tuple["InterfaceRealization", ...]
    dependency_edges: tuple[tuple[str, str], ...]
    status: str
    evidence_summary: "RealizationEvidenceSummary"
    warnings: tuple[str, ...]
    schema_version: str
```

Allowed statuses are:

```text
fully_realized
partially_realized
unrealized_but_plausible
contradicted
not_attempted
```

Only `fully_realized` means that every required interface has a coherent,
validated transformation path. It still does not guarantee laboratory success.

## 5. Architecture and ownership

The feature should extend existing package boundaries rather than create a
parallel retrosynthesis system.

```text
reactive_taxonomy
  owns molecular observations, atom environments, graph edits,
  functional groups, reaction signatures, and graph-derived complexity

core_retrosynthesis
  owns partition contracts, landscape generation, operator application,
  bounded realization search, route trees, ranking, and route evidence

condition_registry
  owns condition-substance and canonical-recipe identity

condition_recommender
  provides structure-backed precedent and condition evidence where needed

chem_coworker / app
  renders and orchestrates; contains no partition or route chemistry rules
```

Existing contracts to reuse:

| Existing capability | Use in this design |
| --- | --- |
| Generic operator ladder and `GenericDisconnectionCandidate` | Generate and validate explicit reverse transformations. |
| `ReactionRouteTree` | Canonical observed/planned route realization graph. |
| `MultistepRetrosynthesisRoute` | Existing bounded planned-route result and evidence. |
| Strategic-complexity assessment | Inspectable partition and expansion evidence, not a feasibility claim. |
| Strategy and synthon identities | Cluster tactical variants after validation. |
| Reaction compatibility and condition evidence | Route-level cautions and later ranking evidence. |

Suggested new modules:

```text
core_retrosynthesis/
  synthetic_partition.py       immutable contracts and canonical identity
  partition_projection.py      project route trees into target partitions
  partition_landscape.py       enumerate and rank target-level views
  partition_realization.py     realize a selected partition
  partition_evaluation.py      benchmark and review records

core_retrosynthesis/definitions/
  synthetic_partition_policy.v1.json
```

Do not duplicate molecule parsing, atom mapping, operator application, route
trees, condition resolution, or precedent retrieval in these modules.

## 6. Landscape generation

### 6.1 Candidate sources

Generate candidate interfaces from several evidence sources:

1. product-side sites from applicable, validated generic operators;
2. interfaces projected from structurally valid observed routes;
3. versioned retron or strategic-seam definitions;
4. low-priority graph seams used only as proposals;
5. optional model or LLM proposals that reference existing atom and candidate
   identifiers.

Every source retains its provenance and evidence class. A structural seam is
not automatically a credible reaction interface.

### 6.2 Partition enumeration

For a configurable range such as `k = 1..6`:

1. normalize and atom-map the target;
2. generate candidate interfaces;
3. build compatible interface combinations;
4. calculate the induced target-atom components;
5. canonicalize the unordered module sets;
6. reject duplicate atom ownership or incomplete coverage;
7. retain supported, ambiguous, and unresolved partitions separately;
8. cluster tactical variants sharing the same target partition and strategic
   edits;
9. select a diverse landscape across `k`, interface sets, ring strategies, and
   evidence classes.

`k = 1` is emitted directly and may accumulate unary-transformation or
no-disconnection evidence.

### 6.3 Symmetric evaluation

No module receives a privileged position in the base algorithm. Partition
features should include:

- exact atom coverage and exclusivity;
- connectedness of each target module;
- interface support and weakest-interface evidence;
- module complexity distribution;
- stereochemical and ring-boundary cautions;
- trivial-fragment burden;
- number of unresolved motifs;
- observed-route support;
- diversity relative to already retained partitions.

A core-retention score may be used by an optional named search policy, but it is
not a universal field or validation rule.

### 6.4 Ranking

Initial ranking should be evidence-tiered rather than probability-like:

```text
1. structural validity
2. interface evidence class
3. weakest-interface evidence
4. realization coverage, when attempted
5. heuristic strategic utility
6. diversity
```

Heuristic component scores must be labeled as heuristics. They must not be
called probabilities until calibrated on held-out review data.

## 7. Route realization

### 7.1 Realization objective

Given a selected partition, find an explicit route tree whose frontier and atom
provenance expose the requested modules, possibly through latent precursor
states.

For desired modules `M_i` and frontier components `F_j`, projection compares
target-derived atoms only. Tactical atoms do not change module ownership.

A realization match must consider:

- atom-weighted module overlap;
- strategic-boundary agreement;
- interface coverage;
- latent-state provenance;
- unresolved or implicit donor atoms.

Mean unweighted Jaccard similarity alone is insufficient because it gives a
one-atom module the same importance as a large structural module.

### 7.2 Symmetric bounded search

The search starts from the target and may expand any unresolved frontier
component. Expansion priority is dynamic and may use:

- unresolved target complexity;
- applicable validated transformations;
- terminal or literature evidence;
- expected progress toward the selected partition;
- route depth and branching cost;
- compatibility warnings.

Detached components are not automatically frozen. A component is terminal only
when supported by explicit terminal evidence or by a clearly labeled bounded
search policy. This permits convergent routes and independent preparation of
multiple complex modules.

An optional `anchor_guided` policy may prefer one annotated region for
medicinal-chemistry workflows. It must produce the same contracts and expose
the policy in diagnostics.

### 7.3 Accepted route edges

Every executable route edge must be one of:

- a validated generic one-step operator application;
- an admitted, independently validated coupled action that retains its physical
  reaction steps; or
- an explicitly non-executable strategic hypothesis.

Executable edges require graph validity, mapping consistency, source-defined
operator identity, and forward-signature agreement. Strategic hypotheses never
count as completed realization.

### 7.4 Dependencies and latent states

The realization tracks:

- which step creates or consumes each handle;
- protection and deprotection state;
- oxidation state;
- preservation or transformation of defined stereochemistry;
- target-atom provenance through every precursor;
- hard and soft step-order constraints;
- at least one valid forward topological order.

A cycle in hard dependencies rejects the realization. Unresolved compatibility
produces review evidence rather than silent acceptance.

### 7.5 Stopping and abstention

Stop when one of the following is true:

- all requested interfaces are realized and frontier modules are terminal;
- the requested partition is only partially exposed within the search bound;
- no validated expansion remains;
- a hard contradiction is found;
- the configured resource limit is reached.

Always distinguish chemical rejection from bounded-search failure and missing
operator coverage.

## 8. Evidence and validation

### 8.1 Evidence levels

Use explicit categorical evidence:

| Level | Meaning |
| --- | --- |
| `E4` | Exact observed route or intermediate realization. |
| `E3` | Same validated interface with close local environments. |
| `E2` | Same reaction-center pattern with relevant structural similarity. |
| `E1` | Broad strategic analogy or defined retron only. |
| `E0` | Unsupported proposal. |

Evidence level and realization status are independent. An exact operator
application may still have weak precedent transfer, and a strategically
attractive partition may remain unrealized.

### 8.2 Structural validation

Every accepted partition must have:

- a valid normalized target graph;
- exact coverage of target atoms;
- no duplicated atom ownership;
- deterministic identity under atom serialization and module permutation;
- valid interface atom references;
- retained stereochemical and aromaticity information.

### 8.3 Realization validation

Every fully realized plan must have:

- valid precursor molecular graphs;
- consistent atom correspondence;
- no unexplained creation or deletion of target atoms;
- validated one-step transformations;
- an acyclic hard dependency graph;
- at least one forward ordering;
- no unresolved hard handle conflict;
- explicit evidence for every interface;
- an identified weakest interface and warnings.

## 9. Output

The main API should return the landscape even when realization is not requested:

```http
POST /v1/synthetic-partitions
```

```json
{
  "target_smiles": "...",
  "k_range": {"min": 1, "max": 6},
  "return_count": 12,
  "include_optional_annotations": true
}
```

A selected partition can be realized separately:

```http
POST /v1/synthetic-partitions/realize
```

```json
{
  "target_smiles": "...",
  "partition_id": "PART1-...",
  "max_depth": 6,
  "maximum_expansions": 1000
}
```

The user-facing view should show:

1. a target depiction colored by role-neutral module membership;
2. optional role labels visually separated from structural facts;
3. strategic interfaces and their evidence;
4. alternative latent states for each module;
5. zero or more route realization trees;
6. status, weakest link, warnings, and search coverage.

## 10. Implementation plan

Implementation must follow the repository's current validation gates. New
partition work begins as a review-only projection and must not alter production
ranking until the existing baseline, blind review, disagreement resolution, and
untouched evaluation gates are complete.

### Current status — 2026-09-01

| Phase | Status | Implementation |
| --- | --- | --- |
| Phase 0 | Complete | The pre-implementation commit, runtime, contract versions, panel hashes, and 1,193-test baseline are frozen in `synthetic_partition_phase0_baseline.v1.json`. |
| Phase 1 | Implemented, review-only | Immutable role-neutral contracts, deterministic identities, mapped reaction projection, route-frontier projection, latent states, JSON round trips, and static HTML review. |
| Phase 2 | Implemented, review-only | Validated-operator partition seeds, tactical-variant clustering, bounded interface combinations, evidence tiers, diversity limits, and explicit abstention. |
| Phase 3 | Not started | Partition-constrained multistep realization remains separate from the Phase 2 landscape. |

Phase 2 combinations show that independently supported interfaces induce a
target partition. They are labeled `operator_combination_unrealized` and do not
claim that the interfaces coexist in one route.

### Phase 0 — Freeze scope and baseline

1. Record current route-tree, operator, reaction-signature, and ranking versions.
2. Run the complete deterministic test suite and existing definition validators.
3. Select a fixed development panel and an untouched test panel.
4. Define review questions before examining system output.

Exit criterion: reproducible baseline artifacts and hash-bound panels.

### Phase 1 — Contracts and observed-route projection

1. Implement immutable partition, module, interface, latent-state, and
   realization contracts.
2. Implement permutation-invariant IDs and serialization.
3. Project target-atom membership from existing `ReactionRouteTree` frontiers.
4. Produce partitions at every observed route depth, including `k = 1` and
   `k = 2`.
5. Export JSON and a simple module-colored HTML review.

This phase uses observed routes and does not generate new chemistry.

Exit criterion: exact atom coverage, deterministic round trips, and chemist-
readable projections on the development panel.

### Phase 2 — Operator-derived partition landscape

1. Generate candidate interfaces from validated one-step operator results.
2. Add compatible interface-set enumeration across a bounded `k` range.
3. Add evidence tiers, duplicate clustering, and strategic diversity.
4. Retain structural seams and retron-only candidates in separate, non-
   executable evidence classes.
5. Add explicit abstention and unresolved-coverage output.

Exit criterion: the system produces valid, symmetric landscapes without forcing
fragmentation and without treating heuristic seams as reactions.

### Phase 3 — Partition-constrained realization

1. Reuse bounded multistep search and the canonical planned route tree.
2. Use symmetric frontier expansion and do not introduce a core-only path.
3. Project target atoms through every generated precursor.
4. Add latent-state and module-membership tracking.
5. Add partition-match, interface-coverage, dependency, and realization status.
6. Retain partial routes and diagnose missing operators or resource limits.

Exit criterion: fully realized results round-trip through validated physical
steps; partial and contradicted results cannot be confused with them.

### Phase 4 — Precedent, conditions, and review

1. Retrieve structure-backed precedents for every realized interface.
2. Add recipe and condition evidence after structural validation.
3. Assess handle coexistence, selectivity cautions, and protection burden.
4. Generate a blind chemist review packet with negative and abstention cases.
5. Calibrate ranking only on the development and validation partitions.

Exit criterion: reviewed evidence supports the displayed ordering and no hard-
incompatible realization is promoted.

### Phase 5 — Optional strategic proposal layer

Only after deterministic validation:

1. allow an LLM or learned model to propose candidate partitions or latent
   states using existing atom and candidate IDs;
2. validate every proposal with deterministic code;
3. measure incremental coverage and false-proposal rate;
4. keep model-only proposals out of executable realization status.

Exit criterion: the proposal layer improves held-out strategic coverage without
reducing validity or abstention quality.

## 11. Evaluation

### 11.1 Required baselines

Compare:

```text
single-step operator portfolio
ordinary bounded multistep search
observed-route frontier projection
symmetric partition landscape
partition-constrained realization
optional anchor-guided realization
```

Use ablations for interface evidence, precedent evidence, compatibility checks,
and optional annotations.

### 11.2 Metrics

Partition metrics:

- exact target-atom coverage and exclusivity;
- permutation and serialization invariance;
- expert usefulness at top 1 and top k;
- strategic-boundary precision and recall;
- useful coverage across `k`;
- trivial-module rate;
- abstention correctness;
- pairwise strategic diversity;
- inter-reviewer agreement.

Realization metrics:

- full and partial realization rates with denominators;
- validated interface coverage;
- precursor-to-product round-trip validity;
- dependency-graph validity;
- search nodes and time to first valid realization;
- weakest-interface evidence distribution;
- unsupported-step count;
- hard-incompatible result count;
- expert preference against ordinary multistep search.

Report results separately for medicinal-chemistry appendage installation,
convergent fragment coupling, late core formation, ring construction,
rearrangement or cascade chemistry, protection-heavy routes, and targets with no
useful multi-module partition.

### 11.3 Leakage control

Partition and route evaluation must be patent-, source-, and scaffold-disjoint
where possible. A precedent used to extract an operator cannot also be counted
as independent evidence for held-out transfer without explicit disclosure.

## 12. Tests

Add deterministic tests for:

- `k = 1`, `k = 2`, and `k >= 3` partitions;
- exact atom coverage and duplicate-ownership rejection;
- module-order and reactant-order invariance;
- symmetric modules with no required role annotations;
- optional annotations not changing partition identity;
- multiple latent states for one target module;
- one frontier component containing several joined modules;
- unary transformations that preserve `k`;
- convergent expansion of more than one frontier component;
- atom donors and precursor-only handles;
- ring-forming and stereochemical interfaces;
- ambiguous, unsupported, conflicting, and contradicted evidence;
- complete, partial, failed, and unattempted realization;
- deterministic JSON round trips and definition validation.

## 13. First-release acceptance criteria

A review-ready first release should:

1. accept a stereochemical target structure;
2. return a valid landscape across the configured `k` range;
3. permit `k = 1`, `k = 2`, and abstention without penalty;
4. preserve exact target-atom provenance;
5. remain invariant to module order and molecular serialization;
6. keep optional roles outside structural identity;
7. reuse validated operators and the canonical route-tree contract;
8. distinguish every realization status clearly;
9. expose evidence, limitations, and search diagnostics;
10. demonstrate improved expert strategic coverage or route selection on a
    leakage-controlled held-out panel.

Production promotion additionally requires completion of the repository's
blind-review and untouched-evaluation gates.

## 14. Summary

The system maintains three linked but distinct objects:

```text
SyntheticPartitionLandscape
  all retained role-neutral views of the target across k >= 1

SyntheticPartition
  one unordered partition of target-derived atoms and its interfaces

StrategicRouteRealization
  one dependency-aware sequence of physical transformations using latent
  precursor states to realize that partition
```

The design deliberately avoids selecting a core too early. It first asks:

> What credible construction views are visible in this molecular graph?

It then asks:

> Which of those views can be realized by validated chemistry, in which latent
> forms and assembly orders, and with what evidence and uncertainty?

This separation provides a broad molecular planning picture while keeping route
claims explicit, testable, and chemically grounded.
