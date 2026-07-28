# Connectivity-First Reaction Grammar: Phase 1 Implementation Plan

## Status

Implemented on 2026-07-28.

Implementation result:

- added the typed connectivity contracts in `reaction_models.py`;
- added canonical graph construction and the `CEG1` shadow identity in
  `reaction_connectivity.py`;
- dual-wrote mapped, reconstructed, and inferred observations inside
  `EditNormalizationResult`;
- preserved current `ReactionEdit`, `ReactionAnalysis`, signature identity,
  admission, retrieval, and recommendation behavior;
- added `connectivity_shadow_report.json` to reaction-edit evaluation;
- added exact, projected, unresolved, hydrogen, charge, inference, symmetry,
  aromatic-domain, topology, and canonicalization-overflow regressions; and
- passed the complete repository suite with 493 tests.

The 12-case reaction-edit benchmark passed its machine gate with 100% shadow
coverage and 100% compatibility parity. It observed 21
`observed_product`, two `main_product_projection`, and three
`exact_reconstruction` scoped changes, with no unsupported bond-domain or
canonicalization-overflow cases.

This plan implements the observation foundation from
[`connectivity_first_reaction_grammar_design.md`](connectivity_first_reaction_grammar_design.md).
It intentionally does not replace the current reaction grammar operators,
signature identity, admission logic, retrieval behavior, or recommendation
indices.

## 1. Phase objective

Introduce an evidence-scoped connectivity observation that can distinguish:

- exact bond transitions between corresponded atoms;
- attachment loss caused by an endpoint being absent from the reported product;
- exact main-product reconstruction from a registered operator;
- conservative correspondence inference;
- schema-level hydrogen changes;
- observed formal-charge changes; and
- unresolved bond states.

The new observation is dual-written beside the current `ReactionEdit` result.
Existing public chemistry behavior remains authoritative during this phase.

Phase 1 answers:

```text
What connectivity and atom-state evidence do we actually have?
```

It does not yet answer:

```text
Can every grammar be executed by one generic rewrite engine?
```

That is the next phase after the observation contract passes its gates.

## 2. Required outcomes

Phase 1 is complete when:

1. Mapped bond changes with both endpoints present receive exact before/after
   bond states.
2. A reactant attachment whose endpoint is absent from the reported product is
   represented as `endpoint_absent` with `main_product_projection` scope, not as
   an observed product-side `no_bond`.
3. Operator-predicted changes retain `exact_reconstruction` scope and are not
   upgraded to mapped observations.
4. Hydrogen-count changes are aggregated into typed `HydrogenDelta` objects.
5. Formal-charge changes between mapped atoms are retained as
   `AtomStateTransition` objects.
6. A canonical connectivity edit graph and shadow key are deterministic,
   atom-order invariant, and reactant-component-order invariant.
7. Existing `ReactionEdit`, `ReactionSignature.signature_id`, serialized
   `ReactionAnalysis`, admission, retrieval, and recommendation behavior remain
   unchanged.
8. The full repository test suite passes.

## 3. Explicit non-goals

Phase 1 does not:

- add `apply_connectivity_rewrite`;
- change reaction grammar JSON;
- migrate or remove `center_replacement`, `join_two_anchors`,
  `pair_addition`, `pair_elimination`, or `change_bond_order`;
- merge existing site types;
- add new reaction-family support;
- change `EditArchetype`;
- use the shadow key in production signature identity;
- bump `REACTION_SIGNATURE_SCHEMA_VERSION`;
- change recommendation admission, retrieval, scoring, or index formats;
- infer omitted byproducts;
- normalize aromatic bond arithmetic;
- model coordination or dative bonds;
- infer radical or isotope changes; or
- promote agents to atom-contributing reactants.

These boundaries prevent an observation-contract change from becoming a
simultaneous chemistry migration.

## 4. Phase 1 contracts

### 4.1 `BondState`

Add to `reactive_taxonomy/reaction_models.py`:

```python
BondStateKind = Literal[
    "bond",
    "no_bond",
    "endpoint_absent",
    "unknown",
]


@dataclass(frozen=True)
class BondState:
    state_kind: BondStateKind
    order: Optional[str]
```

Validation:

- `bond` requires an allowed typed bond order.
- `no_bond`, `endpoint_absent`, and `unknown` require `order=None`.
- Phase 1 arithmetic supports `SINGLE`, `DOUBLE`, and `TRIPLE`.
- `AROMATIC`, `DATIVE`, and other domains may be stored as before/after typed
  states only when exact, but do not receive `delta_units`.
- Unsupported transitions receive a warning rather than approximate arithmetic.

### 4.2 Observation scope

```python
ConnectivityObservationScope = Literal[
    "observed_product",
    "main_product_projection",
    "exact_reconstruction",
    "correspondence_inference",
    "unresolved",
]
```

Scope meanings:

| Scope | Meaning |
|---|---|
| `observed_product` | Before and after states are supported by validated atom correspondence |
| `main_product_projection` | An attachment is absent from the reported main product, but full-product fate is not observed |
| `exact_reconstruction` | A current registered operator exactly reconstructed the reported main product |
| `correspondence_inference` | A conservative unmapped correspondence supports the transition |
| `unresolved` | The available evidence does not establish a definite transition |

Scopes are ordered by evidence policy, but code must not silently convert one
scope to another.

### 4.3 `BondTransition`

```python
@dataclass(frozen=True)
class BondTransition:
    atom_1: ReactionAtomReference
    atom_2: ReactionAtomReference
    before_state: BondState
    after_state: BondState
    delta_units: Optional[int]
    observation_scope: ConnectivityObservationScope
    evidence: str
    confidence: float
```

Rules:

- Endpoint ordering is canonical.
- `delta_units` is derived, not independently authored.
- `delta_units` is present only for localized covalent states with definite
  before and after orders.
- `endpoint_absent` never means `no_bond`.
- A transition retains the reactant-origin atom references used by the current
  edit schema. Product correspondence remains available through atom-map
  evidence during Phase 1.
- Product-only atoms remain a completeness error; Phase 1 does not invent
  reactant provenance.

### 4.4 `HydrogenDelta`

```python
@dataclass(frozen=True)
class HydrogenDelta:
    atom: ReactionAtomReference
    before_count: int
    after_count: int
    delta_count: int
    observation_scope: ConnectivityObservationScope
    evidence: str
    confidence: float
```

One object records the complete count change on one atom. The compatibility
converter may expand `abs(delta_count)` into the existing per-hydrogen
`ReactionEdit` representation.

A hydrogen delta does not identify a transferred physical hydrogen atom.

### 4.5 `AtomStateTransition`

```python
@dataclass(frozen=True)
class AtomStateTransition:
    reactant_atom: ReactionAtomReference
    product_atom: ReactionAtomReference
    before_formal_charge: int
    after_formal_charge: int
    before_radical_electrons: Optional[int]
    after_radical_electrons: Optional[int]
    before_isotope: Optional[int]
    after_isotope: Optional[int]
    observation_scope: ConnectivityObservationScope
    evidence: str
    confidence: float
```

Phase 1 populates formal charge only. Radical and isotope fields are retained
as `None` unless exact values are deliberately added with tests. This avoids a
future incompatible contract change while keeping the first implementation
bounded.

### 4.6 `ConnectivityEditGraph`

```python
@dataclass(frozen=True)
class ConnectivityEditGraph:
    bond_transitions: Tuple[BondTransition, ...]
    hydrogen_deltas: Tuple[HydrogenDelta, ...]
    atom_state_transitions: Tuple[AtomStateTransition, ...]
    edit_component_keys: Tuple[str, ...]
    shadow_key: str
    evidence: str
    confidence: float
    warnings: Tuple[str, ...]
    schema_version: str = "1.0"
```

`edit_component_keys` are graph-connected edit observations, not claims about
mechanistic event count or temporal order.

The graph is attached to the internal `EditNormalizationResult` with a default
of `None`. It is not attached to `ReactionSignature` or `ReactionAnalysis`
during Phase 1.

## 5. Implementation work packages

### WP1: typed models and validation

Files:

- `reactive_taxonomy/reaction_models.py`
- new `reactive_taxonomy/reaction_connectivity.py`
- `reactive_taxonomy/__init__.py`
- `tests/reactive_taxonomy/test_reaction_connectivity_models.py`

Tasks:

1. Add the typed contracts from Section 4.
2. Add constructors or validators for legal state combinations.
3. Add deterministic endpoint ordering helpers.
4. Add localized bond-order rank conversion in one place.
5. Reject contradictory states and invalid confidence values.
6. Export only stable public types; keep construction helpers module-local
   unless another package needs them.

Acceptance:

- model construction tests cover all state kinds;
- invalid state/order combinations fail deterministically;
- no new imports from legacy `chemtools`;
- no direct `Chem.MolFromSmarts()` calls are introduced.

### WP2: mapped-transition extraction

Files:

- `reactive_taxonomy/reaction_edits.py`
- `reactive_taxonomy/reaction_connectivity.py`
- `tests/reactive_taxonomy/test_mapped_connectivity_transitions.py`

Tasks:

1. Reuse the current validated `_MappedSide` data.
2. Partition mapped atoms into:
   - present on both sides;
   - reactant-only;
   - product-only.
3. For a bond whose endpoints are both present on both sides:
   - create exact `bond` or `no_bond` before/after states;
   - set `observation_scope="observed_product"`;
   - derive localized `delta_units`.
4. For a reactant bond with at least one reactant-only endpoint:
   - retain the reactant `bond` state;
   - set the after-state to `endpoint_absent`;
   - set `observation_scope="main_product_projection"`;
   - leave `delta_units=None`.
5. Do not construct a transition for product-only atoms as though their
   reactant state were `no_bond`; retain the existing completeness and mapping
   warnings.
6. Aggregate mapped hydrogen-count changes.
7. Compare formal charge for atoms present on both sides.
8. Preserve invalid-map, duplicate-map, element-mismatch, and partial-map
   behavior.

Acceptance:

- mapped C=C to C–C produces one exact transition with `delta_units=-1`;
- mapped C–Br with Br absent from products produces a projected
  `endpoint_absent` transition;
- mapped C–N formation produces an exact `delta_units=+1` transition;
- product-only atom maps never become invented formed-bond provenance;
- current `ReactionEdit` output remains byte-for-byte equivalent for the
  existing regression fixtures.

### WP3: reconstruction and inferred-transition adapters

Files:

- `reactive_taxonomy/reaction_edits.py`
- `reactive_taxonomy/reaction_connectivity.py`
- `tests/reactive_taxonomy/test_reconstructed_connectivity_transitions.py`

Tasks:

1. Convert current operator `BondChange` results into transitions with
   `observation_scope="exact_reconstruction"`.
2. Preserve current operator evidence strings separately from scope.
3. Convert conservative scaffold/global correspondence results with
   `observation_scope="correspondence_inference"`.
4. Never upgrade reconstructed or inferred transitions to
   `observed_product`.
5. Reconcile mapped and reconstructed transition keys:
   - exact agreement;
   - agreement modulo projection scope;
   - definite conflict;
   - not comparable because one state is unknown.
6. Add typed warnings rather than treating every non-identical set as a hard
   chemical contradiction.

Proposed reconciliation warnings:

```text
MAPPING_RECONSTRUCTION_TRANSITION_CONFLICT
MAPPING_RECONSTRUCTION_SCOPE_DIFFERENCE
PROJECTED_ATTACHMENT_NOT_FULLY_OBSERVED
UNSUPPORTED_BOND_DOMAIN
ATOM_STATE_RECONSTRUCTION_UNSUPPORTED
```

Acceptance:

- Suzuki reconstruction remains reconstruction-scoped even when the main
  product matches exactly;
- Br2 addition with all heavy endpoints retained can agree exactly with mapped
  transitions;
- Si–H addition produces heavy-bond transitions plus hydrogen deltas without a
  fictitious mapped H atom;
- beta elimination preserves projected leaving-group loss;
- existing candidate selection does not change.

### WP4: compatibility projection

Files:

- `reactive_taxonomy/reaction_connectivity.py`
- `reactive_taxonomy/reaction_edits.py`
- `tests/reactive_taxonomy/test_connectivity_edit_compatibility.py`

Tasks:

1. Add a one-way converter from the connectivity observation to current
   `ReactionEdit` views.
2. Exact localized states map to `formed`, `broken`, or `order_changed`.
3. Hydrogen deltas expand to the current per-unit `hydrogen_change` edits.
4. Projected attachment losses reproduce current `broken` compatibility edits
   only for parity, while the connectivity graph retains their projection
   scope.
5. Unknown states produce no definite compatibility edit.
6. Compare converted edits with the existing normalizers during tests and
   optional audit runs.

Do not make `ReactionEdit` the input source for mapped transition construction:
doing so would lose the distinction between `no_bond` and
`endpoint_absent`.

Acceptance:

- conversion is deterministic;
- existing edit ordering remains stable;
- existing display labels and archetype inference remain unchanged;
- any parity exception is documented as a chemistry correction rather than
  accepted by snapshot replacement.

### WP5: canonical edit graph and shadow key

Files:

- `reactive_taxonomy/reaction_connectivity.py`
- `tests/reactive_taxonomy/test_connectivity_edit_graph.py`

Tasks:

1. Build graph nodes from atoms incident to bond, hydrogen, or atom-state
   transitions.
2. Build colored edges from before-state, after-state, transition scope, and
   localized delta.
3. Partition disconnected graph components as `edit_component_keys`.
4. Canonicalize without using:
   - reactant component order;
   - atom serialization index;
   - atom-map number as chemical identity;
   - grammar role name;
   - reaction label;
   - named family.
5. Use atom-map numbers only to establish correspondence within one record.
6. Canonicalize symmetric graphs with deterministic color refinement and a
   bounded tie-resolution strategy.
7. Emit an explicit overflow warning rather than using input order as a
   fallback.
8. Hash the canonical token as a non-production `CEG1:` shadow key.

Suggested node labels:

```text
element
formal charge
aromatic flag
hybridization
hydrogen delta
atom-state transition summary
```

Suggested edge labels:

```text
before state
after state
delta units or NONE
observation scope
```

Local-environment identity should be available for a more specific shadow view
but excluded from the base graph topology key.

Acceptance:

- partner order and component order do not change the key;
- swapping equivalent alkene endpoints does not change the key;
- graphs with the same edit counts but different atom-sharing topology have
  different keys;
- projected and fully observed attachment losses have different keys;
- multi-event edit-component ordering is deterministic;
- canonicalization overflow is explicit and deterministic.

### WP6: internal dual-write integration and audit

Files:

- `reactive_taxonomy/reaction_edits.py`
- existing reaction-edit evaluation code under `reactive_taxonomy/`
- `tests/reactive_taxonomy/test_reaction_connectivity_integration.py`

Tasks:

1. Add `connectivity_edit_graph: Optional[ConnectivityEditGraph] = None` to
   `EditNormalizationResult`.
2. Populate it for mapped, exactly reconstructed, and conservatively inferred
   results.
3. Keep `ReactionSignature` construction based on existing `ReactionEdit`
   fields.
4. Do not include the shadow key in `signature_id`, serialized
   `ReactionAnalysis`, recommendation records, or indices.
5. Extend the reaction-edit evaluation artifact with an optional internal
   comparison report:
   - existing edit token;
   - shadow connectivity key;
   - evidence scopes;
   - parity status;
   - warnings.
6. Record counts of:
   - exact observed transitions;
   - projected attachment losses;
   - reconstructed transitions;
   - correspondence-inferred transitions;
   - unsupported bond domains;
   - canonicalization overflow;
   - current-edit compatibility mismatches.

Acceptance:

- ordinary API and CLI serialized output is unchanged;
- signature IDs are unchanged for the existing deterministic-ID fixtures;
- the audit can be regenerated deterministically;
- no recommendation behavior changes.

## 6. Required chemistry fixtures

The Phase 1 fixture set should include at least:

| Case | Required observation |
|---|---|
| Mapped C=C hydrogenation | Exact double-to-single transition and H gains |
| Mapped C=O reduction | Exact multiplicity transition; charge/H state preserved |
| Suzuki coupling | Exact formed C–C; Br/B attachment losses remain projected when endpoints disappear |
| C–N coupling | Exact formed C–N; C–X projection and N–H delta distinguished |
| Br2 alkene addition | Exact C=C decrement, Br–Br loss, and two C–Br formations when all atoms are present |
| R3SiH alkene addition | C–Si formation plus H deltas without fictitious H mapping |
| Beta-halo elimination | C–C increment, projected C–X loss, H loss |
| Intramolecular closure | Formed bond scope and ring topology remain unchanged |
| Formal-charge change | `AtomStateTransition` produced without inventing a bond edit |
| Product-only atom | Completeness/mapping warning; no invented reactant provenance |
| Invalid duplicate map | Existing invalid-map behavior retained |
| Aromatic ring edit | Exact typed state or explicit unsupported-domain warning; no 1.5 arithmetic |
| Symmetric addition | One canonical shadow key |
| Multi-event substitution | Deterministic edit-component multiset |

Every category needs positive, negative, ambiguous, and conflicting-evidence
coverage where applicable.

## 7. Test and validation commands

During development:

```powershell
pytest -q tests/reactive_taxonomy/test_reaction_connectivity_models.py
pytest -q tests/reactive_taxonomy/test_mapped_connectivity_transitions.py
pytest -q tests/reactive_taxonomy/test_reconstructed_connectivity_transitions.py
pytest -q tests/reactive_taxonomy/test_connectivity_edit_compatibility.py
pytest -q tests/reactive_taxonomy/test_connectivity_edit_graph.py
pytest -q tests/reactive_taxonomy/test_reaction_connectivity_integration.py
pytest -q tests/reactive_taxonomy
```

Before handoff:

```powershell
pytest -q
git diff --check
```

Capture the complete-suite baseline at implementation start. Do not hard-code a
stale expected test count into the plan.

## 8. Promotion gates

Phase 1 must pass all gates before Phase 2 begins.

### Gate A: observation correctness

- No absent product endpoint is encoded as an exactly observed `no_bond`.
- Product-only atoms do not receive invented reactant provenance.
- Hydrogen and formal-charge changes are retained independently.
- Unsupported bond domains remain explicit.

### Gate B: compatibility parity

- Existing `ReactionEdit` fixtures are unchanged.
- Existing display labels and `EditArchetype` results are unchanged.
- Current mapping/reconstruction conflicts remain visible.
- Any intentional chemistry correction has a dedicated regression test and
  documented impact.

### Gate C: identity stability

- Existing `signature_id` values remain unchanged.
- The shadow key is deterministic and order invariant.
- Shadow-key collisions and splits are reported and chemically reviewed.
- The shadow key is not used in retrieval.

### Gate D: system isolation

- No legacy `chemtools` imports.
- No `condition_registry` or `condition_recommender` dependency is added to
  `reactive_taxonomy`.
- No network or mutable global state is introduced.
- No arbitrary executable code is loaded from definitions.

### Gate E: regression and performance

- The complete suite passes.
- Transition extraction remains bounded by mapped bonds and current candidate
  counts.
- Canonicalization has an explicit size bound and overflow behavior.
- Batch featurization performance is measured against the starting baseline.

## 9. Suggested pull-request sequence

### PR 1: contracts and canonical state validation

- typed models;
- state validators;
- endpoint ordering;
- unit tests;
- no integration.

### PR 2: mapped observation extraction

- exact versus projected transition extraction;
- hydrogen aggregation;
- formal-charge transitions;
- mapped fixtures.

### PR 3: reconstruction and compatibility adapters

- current operator and correspondence adapters;
- reconciliation rules;
- compatibility parity tests.

### PR 4: edit graph and shadow canonicalization

- graph construction;
- edit components;
- canonical token and key;
- symmetry and ordering tests.

### PR 5: internal dual-write audit

- `EditNormalizationResult` integration;
- evaluation report;
- full regression and performance report;
- documentation of Phase 1 findings.

Each PR should be independently testable and should not mix public retrieval
changes with chemistry observation changes.

## 10. Phase 1 exit report

The final Phase 1 report should state:

- contract and internal schema versions;
- number of evaluated reactions;
- counts by observation scope;
- compatibility parity rate;
- signature-ID parity result;
- shadow-key collision and split counts;
- unsupported bond-domain count;
- canonicalization overflow count;
- performance change;
- known unresolved chemistry;
- full test results; and
- recommendation on proceeding to Phase 2.

## 11. Phase 2 handoff

After all Phase 1 gates pass, Phase 2 may introduce the bounded declarative
rewrite compiler and `apply_connectivity_rewrite`.

The first executor pilots should be:

1. Suzuki coupling;
2. C–N release-and-connect;
3. Br2 and R3SiH addition to an alkene through one pair-addition rewrite; and
4. beta elimination.

They collectively test:

- exact and projected bond transitions;
- explicit and virtual-hydrogen links;
- addition and elimination bond multiplicity;
- endpoint permutation;
- product projection;
- operator parity; and
- mapped/reconstruction reconciliation.

Current operators should be removed only after the generic executor reproduces
their products, edits, ambiguities, warnings, and ordering across the complete
regression corpus.
