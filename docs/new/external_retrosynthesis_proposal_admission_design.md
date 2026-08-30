# External retrosynthesis proposal admission design

## 1. Purpose

This document defines a clean public tool for checking a retrosynthetic proposal
that originated outside the canonical search system, including proposals from:

- an LLM;
- a web or literature search;
- a chemist;
- an external single-step retrosynthesis model; or
- an imported route-planning system.

The tool accepts proposed precursor and target structures, reconstructs all
internal chemistry evidence deterministically, and reports whether the proposal
is invalid, unresolved, structurally verified, operator-supported,
precedent-supported, or independently forward-reproduced.

External claims are never treated as validation evidence. In particular, an
external source cannot supply a trusted internal operator ID, strategy ID,
signature, forward-validation status, compatibility disposition, or confidence
class.

## 2. Scope

The first public contract evaluates one concrete single-step proposal:

```text
precursor_1.precursor_2... >> target
```

It answers separate questions:

1. Are all structures parseable and valence-valid?
2. Do the displayed reaction sides agree with any supplied mapped reaction?
3. Can the molecular graph transformation be resolved without inventing atom
   correspondence?
4. Does the transformation yield a complete generic reaction signature?
5. Does that signature match an admitted internal operator?
6. Are admitted precedents available for the operator and local environment?
7. Are there stereochemical, precursor, reaction-regime, or selectivity
   conflicts?
8. Can the proposed precursors independently reproduce the intended product in
   the forward operator system?
9. Are compatible conditions supported by admitted precedents?

These questions must not be collapsed into a single opaque probability.

## 3. Non-goals

The tool does not:

- prove laboratory feasibility, yield, scalability, or safety;
- trust an LLM-generated atom mapping;
- create precursor SMILES, reaction SMARTS, or protecting groups;
- silently repair an invalid proposal;
- infer commercial availability from molecular weight;
- use web popularity as chemical evidence;
- treat a named reaction label as graph evidence;
- make external models part of `reactive_taxonomy`; or
- bypass the canonical transition-provider admission boundary.

Repair, alternate-strategy search, protection-state enumeration, condition
recommendation, and route assembly remain separate actions invoked after the
assessment.

## 4. Package ownership

The public composition tool belongs in:

```text
core_retrosynthesis/external_proposal_assessment.py
```

It may compose existing standalone packages but must not move their chemistry:

| Responsibility | Owning package |
|---|---|
| Parsing, valence, molecular observation | `reactive_taxonomy` |
| Atom-mapping validation and graph edits | `reactive_taxonomy` |
| Generic reaction signature | `reactive_taxonomy` |
| Operator identity and library matching | `core_retrosynthesis` |
| Precursor and reaction compatibility | existing chemistry policy modules |
| Forward reconstruction and competition | `forward_synthesis` |
| Condition identity and recipes | `condition_registry` |
| Condition retrieval and ranking | `condition_recommender` |
| Agent-facing selection and explanation | `chem_coworker` |

The core assessment must be deterministic and network-free. Web retrieval,
reference downloading, and LLM calls happen outside this module and contribute
only typed source metadata or proposed structures.

## 5. Public input contract

```python
@dataclass(frozen=True)
class ExternalProposalSource:
    source_id: str
    source_kind: Literal[
        "chemist",
        "literature",
        "web",
        "llm",
        "external_ssr",
        "imported_route",
    ]
    reference_id: str | None = None
    source_uri: str | None = None
    model_id: str | None = None
    retrieved_at: str | None = None
    quoted_claim: str | None = None


@dataclass(frozen=True)
class ExternalRetrosynthesisProposal:
    target_smiles: str
    precursor_smiles: str
    mapped_reaction_smiles: str | None = None
    proposed_conditions: Mapping[str, Any] | None = None
    sources: tuple[ExternalProposalSource, ...] = ()
    external_proposal_id: str | None = None
    claimed_transformation: str | None = None
    schema_version: str = "1.0"
```

`claimed_transformation` is advisory text. It cannot route admission or override
contradictory structural evidence.

The input deliberately excludes:

```text
operator_id
operator_signature
strategy_id
template_id
forward_validation_status
compatibility_disposition
admission_status
internal confidence
```

Those fields are outputs derived by deterministic code.

## 6. Assessment dependencies

The main entry point receives explicit dependencies and limits:

```python
def assess_external_retrosynthesis_proposal(
    proposal: ExternalRetrosynthesisProposal,
    operator_library: GenericTemplateLibrary,
    *,
    forward_library: ForwardOperatorLibrary | None = None,
    mapping_provider: AtomMappingProvider | None = None,
    condition_evaluator: ConditionEvidenceEvaluator | None = None,
    limits: ExternalProposalAssessmentLimits | None = None,
) -> ExternalRetrosynthesisAssessment:
    ...
```

Dependency injection keeps the chemistry packages network-free and makes the
assessment reproducible. An absent optional dependency produces an explicit
`not_run` gate; it never produces a silent pass.

Suggested bounded limits include:

```python
@dataclass(frozen=True)
class ExternalProposalAssessmentLimits:
    maximum_operator_matches: int = 20
    maximum_precedent_matches: int = 10
    maximum_forward_operators: int = 300
    maximum_forward_products: int = 20
    mapping_timeout_seconds: float | None = None
```

## 7. Evidence-gate contract

Every stage returns a typed gate rather than mutating a candidate status:

```python
GateStatus = Literal[
    "pass",
    "warn",
    "fail",
    "unresolved",
    "not_run",
    "out_of_scope",
]


@dataclass(frozen=True)
class ExternalProposalGate:
    gate_id: str
    status: GateStatus
    summary: str
    evidence_ids: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    definition_id: str | None = None
```

Required gate IDs are:

```text
input_structure
reaction_side_consistency
atom_correspondence
reaction_completeness
reaction_signature
operator_support
precedent_support
stereochemistry
precursor_compatibility
reaction_compatibility
functional_group_selectivity
forward_reproduction
forward_competition
condition_support
```

Gate order is fixed and serialized. Interpretation failure must not erase valid
lower-level observations.

## 8. Overall assessment contract

```python
ExternalProposalStatus = Literal[
    "invalid",
    "ambiguous",
    "structurally_verified",
    "operator_supported",
    "precedent_supported",
    "forward_reproduced",
    "verified_with_cautions",
    "rejected_by_compatibility",
]


@dataclass(frozen=True)
class ExternalRetrosynthesisAssessment:
    assessment_id: str
    proposal_id: str
    canonical_target_smiles: str | None
    canonical_precursor_smiles: str | None
    normalized_mapped_reaction_smiles: str | None
    status: ExternalProposalStatus
    actionable: bool
    gates: tuple[ExternalProposalGate, ...]
    reaction_signature: ReactionSignature | None
    operator_matches: tuple[ExternalOperatorMatch, ...]
    precedent_matches: tuple[StepPrecedentMatch, ...]
    forward_assessment: RouteStepForwardAssessment | None
    compatibility_evidence: ExternalCompatibilityEvidence | None
    condition_evidence: RetrosynthesisConditionEvidence | None
    sources: tuple[ExternalProposalSource, ...]
    warnings: tuple[str, ...]
    schema_versions: Mapping[str, str]
    schema_version: str = "1.0"
```

`actionable=True` means that the proposal may be converted into an admitted
transition candidate for downstream route search. It does not mean that the
reaction is experimentally proven.

## 9. Deterministic pipeline

### 9.1 Canonicalize input structures

The tool first checks that:

- the target is one molecular graph;
- every precursor component parses;
- all structures pass RDKit sanitization and valence checks;
- target and precursor collections canonicalize deterministically; and
- the proposal contains no agents or products hidden in malformed notation.

Failure produces `status="invalid"` and stops structure-dependent gates.

### 9.2 Validate or obtain atom correspondence

Evidence order is:

1. supplied mapped reaction whose structure is independently validated;
2. injected mapper output validated by
   `reactive_taxonomy.validate_external_atom_mapping`;
3. unique internal correspondence inference;
4. ambiguous or unresolved correspondence retained for review.

A mapped reaction must preserve the map-stripped canonical precursor and target
structures exactly. Duplicate maps, incomplete active-center mapping, changed
bond orders caused by the mapper, or changed stereochemistry are failures or
explicit unresolved evidence.

Mapper confidence is provenance, not chemical confidence. Valid external maps
retain the warning that independent validation was required.

### 9.3 Build observation, edits, and signature

The validated mapped reaction, or the unambiguously inferred reaction when no
map is available, is processed through `reactive_taxonomy.featurize_reaction`.

The tool retains:

- parsed components;
- atom correspondence evidence;
- formed, broken, order-changed, and hydrogen edits;
- reaction completeness;
- reaction core;
- stereochemical changes;
- generic reaction signature;
- ambiguity candidates;
- warnings and definition versions.

No signature means the proposal cannot reach `structurally_verified`.

### 9.4 Derive internal operator identity

`core_retrosynthesis.analyze_generic_reaction` derives operator identity from
the validated structural evidence. The tool searches the admitted operator
library by normalized operator signature.

Operator matching levels are:

```text
exact_operator_signature
same_edit_archetype_and_site_environment
same_edit_archetype_only
no_operator_support
```

Only `exact_operator_signature` establishes `operator_supported`. Broader
matches are analogical evidence and must not become trusted internal IDs.

### 9.5 Retrieve admitted precedents

Precedents are searched only among source-round-tripped, admitted library
records. Ranking uses separate evidence fields:

- exact operator agreement;
- product local-environment similarity;
- precursor local-environment similarity;
- stereochemical agreement;
- independent reference support;
- source quality and provenance;
- condition availability.

Exact product or precursor identity is reported separately from similarity.
Similarity alone never verifies a proposal.

### 9.6 Run compatibility and selectivity checks

The proposal is evaluated with existing deterministic policies:

- precursor intramolecular reactive-pair compatibility;
- reaction-regime versus functional-group compatibility;
- functional-group competition at the proposed reaction center;
- stereochemical preservation, creation, loss, or ambiguity;
- strategic-complexity classification.

Any policy disposition of `reject` results in
`status="rejected_by_compatibility"`, even when graph reconstruction succeeds.
Warnings and demotions remain visible without being converted into hard failure.

### 9.7 Independently challenge the forward direction

When a forward library is supplied, the intended product is hidden from the
blind search. The tool then reports:

- whether the intended product is structurally reproduced by the matched
  operator;
- whether blind forward search recovers it;
- intended-product rank;
- the best competing product;
- score margin;
- condition incompatibility; and
- whether enumeration limits were reached.

Forward dispositions retain their existing meanings:

```text
clear
competitive
unsupported
structurally_inconsistent
condition_incompatible
out_of_scope
```

Forward reproduction is strong independent evidence, not a yield prediction.

### 9.8 Assess conditions when available

Conditions are normalized through `condition_registry` and assessed through the
canonical recommender. Raw web or LLM condition strings do not become resolved
substance identities without registry evidence.

Condition results distinguish:

```text
direct precedent support
generic fallback support
incompatible recipe
insufficient evidence
not assessed
```

The absence of condition evidence normally produces a caution rather than
invalidating a structurally verified proposal.

## 10. Status resolution

Overall status is derived from gates in a fixed order:

```text
invalid structure or reaction sides
    -> invalid

unresolved/ambiguous atom correspondence or signature
    -> ambiguous

complete signature, no exact admitted operator
    -> structurally_verified

exact admitted operator, no admitted precedent match
    -> operator_supported

exact admitted operator plus admitted precedent support
    -> precedent_supported

forward target structurally reproduced
    -> forward_reproduced

any hard compatibility rejection
    -> rejected_by_compatibility

otherwise supported but warnings remain
    -> verified_with_cautions
```

The serialized report also includes the strongest achieved evidence tier so a
warning does not erase the fact that lower gates passed.

Suggested actionability policy for v1:

| Status | Actionable by default |
|---|---:|
| `invalid` | no |
| `ambiguous` | no |
| `structurally_verified` | review only |
| `operator_supported` | review only |
| `precedent_supported` | yes, unless a hard gate rejects |
| `forward_reproduced` | yes, unless a hard gate rejects |
| `verified_with_cautions` | configurable, with cautions retained |
| `rejected_by_compatibility` | no |

## 11. Internal-confidence comparison

An external proposal is not less trusted merely because it came from an LLM or
web page. It is less trusted when it lacks the evidence automatically retained
by an internal provider.

The tool may declare `internal_confidence_parity=True` only when:

1. mapped or uniquely inferred atom correspondence is structurally verified;
2. a complete generic signature is reproduced;
3. the signature exactly matches an admitted internal operator;
4. at least one admitted precedent supports the operator;
5. no hard compatibility gate rejects the proposal; and
6. the same admission checks used by `TransitionProviderOrchestrator` pass.

Forward reproduction and condition support are reported as additional evidence.
They are not required for parity with the current operator-ladder provider,
because the provider itself does not prove experimental execution.

## 12. Conversion into a transition candidate

External assessments must not construct `GenericDisconnectionCandidate` before
validation. After admission, a separate function may convert a supported
assessment:

```python
def admitted_transition_from_external_assessment(
    assessment: ExternalRetrosynthesisAssessment,
    selected_operator_match_id: str,
) -> AdmittedExternalTransition:
    ...
```

The conversion copies only deterministic internal fields. External source claims
remain provenance and never populate trusted identity fields.

The resulting transition can enter the provider orchestrator under a registered
provider such as:

```text
external_proposal_assessment
```

Its provider-local score is not compared directly with scores from SSR or
operator-ladder providers.

## 13. Agent-facing behavior

`chem_coworker` may expose two bounded tools:

```text
assess_external_proposal
admit_external_proposal
```

The LLM may supply the external precursor and target strings exactly as obtained
from the source, select a returned operator-match ID, or request further
inspection. It may not edit the structures during admission or invent missing
internal IDs.

If the assessment is ambiguous, the next actions may include:

- request deterministic atom mapping;
- inspect correspondence alternatives;
- search for exact precedents;
- request a registered repair action; or
- retain the proposal as non-actionable advisory evidence.

## 14. Serialization and identity

Proposal identity uses canonical precursor and target structures, not source
wording:

```text
EXTPROP1(canonical_precursors, canonical_target)
```

Assessment identity additionally includes:

- normalized mapped-reaction identity when available;
- operator-library definition/version;
- compatibility-policy versions;
- forward-library version;
- condition-recommender definition/version;
- assessment schema version; and
- bounded-limit definition.

Source URIs, LLM model names, and quoted claims are provenance and do not alter
chemical proposal identity.

## 15. Security and robustness

- Never execute code, SMARTS, or tool names supplied by an external proposal.
- Never dynamically import an external operator.
- Enforce molecule, atom, component, mapping, and enumeration bounds.
- Reject malformed reaction notation before expensive work.
- Preserve raw input separately from canonical structures.
- Avoid embedding full proprietary source documents in assessment artifacts.
- Treat URLs and citations as unverified provenance until resolved by an
  authorized external layer.
- Cache only by normalized chemistry and versioned configuration.

## 16. Tests

### 16.1 Unit cases

Required positive cases:

- valid mapped C-N, C-O, C-S, C-C, reduction, oxidation, and protection-state
  proposals;
- unmapped proposal with uniquely inferred correspondence;
- exact admitted operator match;
- analogical precedent match;
- forward target recovery;
- condition-supported proposal; and
- reactant-order invariance.

Required negative and ambiguous cases:

- invalid valence;
- multiple target products;
- mapper changes structure;
- duplicate or incomplete map numbers;
- mapped/display structure contradiction;
- ambiguous correspondence;
- incomplete reaction signature;
- operator-signature mismatch;
- stereochemical contradiction;
- hard precursor compatibility rejection;
- hard reaction-regime rejection;
- intended product not forward-reproduced;
- competing forward product; and
- unsupported conditions.

### 16.2 Parity cases

Take admitted internal provider candidates, remove their trusted identity fields,
submit only precursor/target structures plus validated mapped reactions, and
require the external tool to recover the same:

- canonical reaction identity;
- generic signature;
- operator signature;
- compatibility disposition;
- admitted operator match; and
- precedent IDs.

The test must also demonstrate that fabricated external operator IDs have no
effect.

### 16.3 Held-out evaluation

Run the tool on:

- held-out observed route steps;
- internal SSR proposals;
- external SSR proposals;
- LLM-generated proposals;
- literature-extracted proposals; and
- a chemist-reviewed challenge panel covering protection, stereochemistry,
  chemoselectivity, rearrangements, and multicomponent chemistry.

Report gate-level coverage, false admission, false rejection, ambiguity,
operator recovery, forward recovery, compatibility adjudication, latency, and
reviewer agreement. Do not tune on the untouched panel.

## 17. Minimal implementation plan

### Milestone 1: structural assessment

- Add immutable input, gate, result, and source contracts.
- Canonicalize structures and validate reaction sides.
- Validate supplied mapped reactions.
- Run observation, completeness, signature, and compatibility gates.
- Return deterministic statuses without constructing a search candidate.

### Milestone 2: internal evidence recovery

- Build an operator-signature index for the loaded generic library.
- Return exact and analogical operator matches separately.
- Retrieve admitted precedents with local-environment and stereochemical
  evidence.
- Add internal-provider parity tests.

### Milestone 3: independent challenge

- Add optional target-hidden forward reconstruction and competition assessment.
- Add optional canonical condition evidence.
- Freeze status resolution and actionability policy from held-out review.

### Milestone 4: orchestration integration

- Add the registered external-proposal provider.
- Expose identifier-based assessment/admission tools to `chem_coworker`.
- Record proposal source, gate outcomes, selected operator match, and route
  contribution in replay telemetry.

## 18. Release gate

The public tool remains advisory until it:

1. reproduces internal provider admission on a frozen parity panel;
2. rejects every mapped/display contradiction case;
3. preserves ambiguous correspondence rather than forcing a signature;
4. demonstrates acceptable false-admission and false-rejection rates under
   blinded chemist review;
5. reports every policy and schema version;
6. passes the complete repository test suite; and
7. does not introduce a second reaction-validation or recommendation path.

Only after these gates pass may an external proposal be inserted into canonical
route search as an actionable transition.

