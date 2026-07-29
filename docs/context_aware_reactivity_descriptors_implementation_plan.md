# Context-Aware Reactivity Descriptors: Implementation Plan

## Implementation status

Completed on 2026-07-28. The typed profile is now the sole active environment
descriptor path. The intentional identity migration is active:

- profile schema `1.0`;
- reaction-signature schema `3.0` and `RS3` identifiers;
- `signature_features.v3.json` and `taxonomy_manifest.v3.json`;
- recommendation-record schema `3.3`;
- generic converter `generic_conversion.v2.3`;
- generic-index, concise-review, and condition-rule schemas `2.0`.

GUI, CLI, review exports, rule facts, signatures, indexing, retrieval, and
similarity consume the same typed profile. The former open `steric` and
`electronic` environment dictionaries and their rule aliases have been removed.

## 1. Objective

Implement a deterministic, graph-derived reactivity profile for every selected
reaction site. The profile should use a consistent outer structure while
allowing aromatic, alkyl, unsaturated, activated-center, and heteroatom contexts
to expose different chemistry-specific factors.

The target organization is:

```text
context identity | steric accessibility | electronic activation |
reactive-center state | modifiers and liabilities
```

The implementation must:

- simplify the chemist-facing reaction summary;
- improve the structured descriptors used by rules and retrieval;
- preserve atom, site, and evidence provenance;
- remain independent of named reaction families;
- distinguish missing, unresolved, and not-applicable descriptors;
- avoid silently changing existing signature IDs, admission, retrieval, or
  rankings during the initial implementation;
- migrate all consumers before removing legacy environment fields.

This plan implements the design in
[`aromatic_reactivity_descriptor_proposal.md`](aromatic_reactivity_descriptor_proposal.md)
and extends it to non-aromatic contexts.

## 2. Scope and non-goals

### In scope

- Typed, immutable reactivity-profile contracts in `reactive_taxonomy`.
- Context dispatch derived from the molecular graph and existing context facets.
- Aromatic, alkyl, alkenyl, alkynyl, acyl, sulfonyl/phosphoryl, and heteroatom
  descriptor variants.
- Orthogonal leaving-group, transfer-carrier, acid/base, and reactivity-liability
  modifiers.
- A shared chemist-facing renderer.
- Canonical serialization into reaction signatures and recommendation records.
- Versioned signature tokens and similarity features after shadow evaluation.
- Condition-rule fact migration.
- Unit, chemistry-regression, invariance, dataset, and performance testing.

### Not in scope for the first activation

- Quantum-chemical electron density, partial charges, or transition-state
  calculations.
- Learned embeddings as replacements for deterministic descriptors.
- Predicting a mechanism or named family from descriptors alone.
- Universal numerical comparison between unrelated contexts.
- Automatic recalibration of condition recommendations without outcome evidence.
- Maintaining both legacy and new descriptor paths permanently.

Any graph-derived electronic score must be described as an interpretable
chemistry prior, not a physical electron-density measurement.

## 3. Architectural decisions

### 3.1 Uniform envelope, typed context variants

All sites receive one `SiteReactivityProfile`, but `context` is a tagged union:

```python
ContextDescriptor = (
    AromaticContextDescriptor
    | AlkylContextDescriptor
    | AlkenylContextDescriptor
    | AlkynylContextDescriptor
    | ActivatedCenterContextDescriptor
    | HeteroatomContextDescriptor
    | OtherContextDescriptor
)
```

The same outer field names make serialization and rendering predictable. The
variant contracts prevent unrelated chemistry from being forced into identical
fields.

### 3.2 Scaffold context and handle modifiers remain orthogonal

Do not model `transfer_partner` or `leaving_group` as scaffold contexts. An
aryl-boron partner is still aromatic at the transferred carbon; boron is an
orthogonal transfer-carrier modifier. Likewise, an alkyl bromide has an alkyl
context plus a leaving-group modifier.

This separation permits:

```text
aromatic context + B(OH)2 transfer carrier
alkyl context + Br leaving group
acyl context + Cl leaving group
heteroatom context + removable H
```

It also prevents named-family assumptions from entering the observation layer.

### 3.3 Observation and interpretation are separate

The base profile stores:

- graph topology;
- atom state;
- bounded local steric measurements;
- recognized electronic contributions;
- handle or carrier modifiers;
- uncertainty and method provenance.

Optional interpretation records may state:

- electrophile activation;
- lone-pair availability;
- acidity class;
- elimination liability;
- transfer liability.

An interpretation must identify its rule and confidence. It must not overwrite
the underlying observation.

### 3.4 Descriptor state is explicit

Every optional descriptor block uses:

```python
DescriptorStatus = Literal[
    "observed",
    "derived",
    "unresolved",
    "not_applicable",
    "not_computed",
]
```

`None` alone must not mean all five states. Retrieval and rules must treat
`unresolved` as uncertainty and `not_applicable` as a legitimate absence.

### 3.5 Signature features use bins, not raw scores

Canonical records retain normalized continuous scores and contribution details.
Reaction-signature identity uses only stable categorical tokens defined by a
versioned vocabulary.

Examples:

```text
context:aromatic
aromatic_family:pyridine
heteroatom:N@ortho:pyridine_like
steric_access:hindered
ortho_occupancy:2_of_2
ortho_burden:high
electronic_demand:slightly_poor
alkyl_center:secondary
alpha_branching:branched
beta_hydrogen:present
lone_pair:pyrrole_like
```

Display labels, source reaction names, raw atom indices, and floating-point
scores must not enter signature identity.

## 4. Target contracts

Add the public contracts under `reactive_taxonomy/descriptors/models.py`.
All public fields require type annotations and docstrings.

### 4.1 Common profile

```python
@dataclass(frozen=True)
class DescriptorEvidence:
    source: str
    method: str
    confidence: float
    contributing_atom_indices: tuple[int, ...] = ()
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class ReactiveCenterProfile:
    atom_index: int
    element: str
    hybridization: str
    formal_charge: int
    radical_electrons: int
    hydrogen_count: int
    heavy_atom_attachment_count: int
    aromatic: bool
    in_ring: bool
    ring_sizes: tuple[int, ...]
    substitution_class: str | None
    conjugation_class: str | None
    lone_pair_class: str | None
    lone_pair_availability: str | None
    acidity_class: str | None
    evidence: DescriptorEvidence


@dataclass(frozen=True)
class StericProfile:
    accessibility_class: Literal["open", "moderate", "hindered", "severe"]
    accessibility_score: float
    approach_burden_class: Literal["none", "low", "medium", "high"]
    branch_contributions: tuple["StericContribution", ...]
    context_metrics: dict[str, int | float | str | bool]
    evidence: DescriptorEvidence


@dataclass(frozen=True)
class ElectronicProfile:
    activation_axis: str
    activation_class: str
    activation_score: float
    contributions: tuple["ElectronicContribution", ...]
    evidence: DescriptorEvidence


@dataclass(frozen=True)
class ReactivityModifier:
    modifier_type: Literal[
        "leaving_group",
        "transfer_carrier",
        "removable_hydrogen",
        "coordination",
        "redox",
        "elimination",
        "strain",
        "other",
    ]
    modifier_id: str
    class_name: str
    attributes: tuple[tuple[str, str], ...]
    evidence: DescriptorEvidence


@dataclass(frozen=True)
class SiteReactivityProfile:
    site_id: str
    center_atom_index: int
    context_kind: str
    context: ContextDescriptor
    reactive_center: ReactiveCenterProfile
    steric: StericProfile
    electronic: ElectronicProfile
    modifiers: tuple[ReactivityModifier, ...]
    flags: tuple[str, ...]
    status: DescriptorStatus
    definition_versions: tuple[tuple[str, str], ...]
    schema_version: str = "1.0"
```

The exact `context_metrics` fields must be validated by context kind. It is a
temporary serialization convenience, not permission for arbitrary keys. If the
fields stabilize during the vertical slices, replace the dictionary with typed
context-specific steric detail records before activation.

### 4.2 Context variants

#### Aromatic

```python
@dataclass(frozen=True)
class AromaticContextDescriptor:
    context_kind: Literal["aromatic"]
    system_class: str
    ring_family: str
    ring_sizes: tuple[int, ...]
    aromatic_ring_count: int
    fused: bool
    anchor_in_ring: bool
    heteroatoms: tuple[AromaticHeteroatom, ...]
    ortho_occupancy_count: int
    ortho_capacity: int
    ortho_burden_class: str
    ortho_burden_score: float
```

#### Alkyl

```python
@dataclass(frozen=True)
class AlkylContextDescriptor:
    context_kind: Literal["alkyl"]
    carbon_substitution: Literal["methyl", "primary", "secondary", "tertiary"]
    alpha_carbon_neighbor_count: int
    alpha_branched: bool
    beta_branch_count: int
    beta_hydrogen_count: int
    cyclic: bool
    ring_sizes: tuple[int, ...]
    benzylic: bool
    allylic: bool
    propargylic: bool
    adjacent_heteroatoms: tuple[str, ...]
```

#### Alkenyl and alkynyl

```python
@dataclass(frozen=True)
class AlkenylContextDescriptor:
    context_kind: Literal["alkenyl"]
    endpoint_substitution: tuple[int, int]
    alkene_class: str
    stereochemistry: str | None
    cyclic: bool
    ring_size: int | None
    conjugation_class: str
    allylic_branch_count: int


@dataclass(frozen=True)
class AlkynylContextDescriptor:
    context_kind: Literal["alkynyl"]
    terminal: bool
    endpoint_substitution: tuple[int, int]
    conjugation_class: str
    propargylic_branch_count: int
```

#### Activated centers

```python
@dataclass(frozen=True)
class ActivatedCenterContextDescriptor:
    context_kind: Literal["acyl", "sulfonyl", "phosphoryl"]
    center_class: str
    attached_group_classes: tuple[str, ...]
    conjugation_class: str
    heteroatom_substitution: tuple[str, ...]
    enolizable: bool | None
```

#### Heteroatoms

```python
@dataclass(frozen=True)
class HeteroatomContextDescriptor:
    context_kind: Literal["heteroatom"]
    element: Literal["N", "O", "S", "P", "other"]
    substitution_class: str
    attached_contexts: tuple[str, ...]
    resonance_class: str
    lone_pair_class: str
    proton_count: int
```

`OtherContextDescriptor` preserves center properties and evidence without
inventing a classification.

## 5. Context-specific factors

| Context | Identity | Steric factors | Electronic factors | Important modifiers |
|---|---|---|---|---|
| Aromatic | Ring family, size, fusion, heteroatoms and relative positions | Ortho occupancy, substituent burden, adjacent fusion | Anchor electron demand, charge, resonance/inductive contributors | Leaving group, transfer carrier, coordination |
| Alkyl | Methyl/1°/2°/3°, cyclic, benzylic/allylic/propargylic | α substitution, α/β branching, ring approach | C-center polarization, adjacent heteroatom or stabilizing group | Leaving group, β-H/elimination, radical stability |
| Alkenyl | Endpoint substitution, E/Z, cyclicity, conjugation | Vinylic and allylic crowding | Alkene polarization and conjugation | Leaving group, addition orientation |
| Alkynyl | Terminal/internal, substitution, conjugation | Propargylic crowding | Terminal acidity and polarization | Transfer carrier, metal coordination |
| Acyl | Acyl derivative, attached group, conjugation | Carbonyl approach shielding | Carbonyl electrophilic activation and resonance donation | Leaving group, enolizability, chelation |
| Sulfonyl/phosphoryl | Center class and attached heteroatoms | Approach shielding | Center activation and charge distribution class | Leaving group, redox and coordination |
| N/O/S center | Element, H count, substitution, attached contexts | Attached-group and α/ortho burden | Lone-pair availability, resonance, acidity/basicity | Deprotonation, coordination, oxidation |

The axes in the `ElectronicProfile` are context-specific:

```text
aromatic: electronic_demand
alkyl: center_polarization
alkenyl/alkynyl: pi_polarization
acyl: carbonyl_activation
heteroatom: lone_pair_availability
```

Scores on different axes must not be directly compared.

## 6. Definitions and validation

Add these versioned definition files:

```text
reactive_taxonomy/definitions/reactivity_descriptor_rules.v1.json
reactive_taxonomy/definitions/aromatic_systems.v1.json
reactive_taxonomy/definitions/reactivity_rendering.v1.json
```

### `reactivity_descriptor_rules.v1.json`

Owns:

- allowed context kinds and descriptor statuses;
- steric distance limits and attenuation rules;
- size and branching weights;
- score normalization and bin boundaries;
- electronic contribution vocabulary and weights;
- context-specific activation axes;
- allowed modifier types and liability flags;
- method identifiers and schema version.

### `aromatic_systems.v1.json`

Owns:

- graph or SMARTS definitions for validated named aromatic systems;
- generic fallback classes;
- heteroatom electronic-role annotations;
- display names.

All SMARTS compilation must use
`reactive_taxonomy.chemistry.smarts_cache.compile_smarts()` or
`compile_smarts_batch()`. No descriptor module may call
`Chem.MolFromSmarts()` directly.

### `reactivity_rendering.v1.json`

Owns:

- fixed display-field order;
- readable category labels;
- omission rules for zero, not-applicable, and low-value diagnostic fields;
- compact and expanded rendering templates.

Python owns graph traversal, ring perception, contribution calculations,
normalization, validation, and deterministic rendering behavior.

### Definition activation

The activation is complete. Descriptor definition hashes are embedded in each
`SiteReactivityProfile` and included in the identity-bearing
`taxonomy_manifest.v3.json`. This intentionally changed signature IDs from the
`RS2` namespace to `RS3`.

## 7. Proposed package changes

### `reactive_taxonomy`

Add:

```text
reactive_taxonomy/descriptors/__init__.py
reactive_taxonomy/descriptors/models.py
reactive_taxonomy/descriptors/registry.py
reactive_taxonomy/descriptors/common.py
reactive_taxonomy/descriptors/aromatic.py
reactive_taxonomy/descriptors/alkyl.py
reactive_taxonomy/descriptors/unsaturated.py
reactive_taxonomy/descriptors/activated_centers.py
reactive_taxonomy/descriptors/heteroatom.py
reactive_taxonomy/descriptors/modifiers.py
reactive_taxonomy/descriptors/rendering.py
reactive_taxonomy/descriptors/tokens.py
```

Modify:

- `models.py`: add `reactivity_profile` to `SiteEnvironment` as an optional
  shadow field, then make it required after migration.
- `reaction_models.py`: add `reactivity_profile` to `ReactionPartner` and
  `ReactionPartnerEnvironment`.
- `environments.py`: delegate new calculations to the descriptor dispatcher;
  retain legacy dictionary construction temporarily.
- `reaction_environments.py`: copy typed profiles without family-specific
  reinterpretation.
- `reaction_signatures.py`: serialize profiles but exclude them from signature
  tokens until the signature activation phase.
- `api.py`, `cli.py`, `validation.py`, and package exports: expose and validate
  the new contracts.

The descriptor package must not import `condition_registry`,
`condition_recommender`, or legacy `chemtools`.

### `condition_recommender`

Modify:

- `conversion/generic.py`: serialize the nested profile.
- `conversion/concise_review.py`: use the shared renderer.
- `signature_features.py`: introduce versioned context-aware tokens during
  signature activation.
- `similarity.py`: compare context-compatible features and continuous values.
- `generic_indexing.py`: index the new categorical environment tokens.
- `rules/models.py`: add typed facts for activated profile fields.
- `rules/facts.py`: project verified profile observations into rule facts.
- `rules/matcher.py`: support only allowlisted declarative predicates.
- rule-definition validation and review output: validate and explain new facts.

### `app`

Modify `app/generic_recommender_gui.py` to consume compact and expanded strings
from the shared renderer. Remove chemistry-specific interpretation from the
application layer.

## 8. Implementation phases

### Phase 0: Baseline and migration inventory

Tasks:

1. Capture current environment JSON for a representative chemistry corpus.
2. Capture L0-L4 keys, `signature_id`, admission status, retrieval level,
   candidate counts, and recommendation ranks.
3. Add characterization tests for every current legacy consumer.
4. Record all direct accesses to `steric`, `electronic`,
   `local_heavy_atoms_r2`, `qualitative_sum`, `attached_groups`, and
   `center_substitution_class`.
5. Select a frozen evaluation set spanning all target contexts.

Acceptance:

- The complete test suite passes before implementation.
- Baseline artifacts are reproducible.
- Every legacy consumer has an owner and migration phase.

### Phase 1: Typed profile foundation in shadow mode

Tasks:

1. Add common models, statuses, evidence, context union, and serialization.
2. Add definition loaders and full schema/vocabulary validation.
3. Implement the context dispatcher using existing graph-derived context facets.
4. Attach an optional shadow `reactivity_profile` to `SiteEnvironment`.
5. Preserve existing `steric` and `electronic` dictionaries unchanged.
6. Ensure `_partner_token()` and signature definition versions remain unchanged.

Acceptance:

- Every usable reactive site receives either a typed profile or an explicit
  unresolved profile.
- Existing L0-L4 keys and `signature_id` values remain byte-for-byte identical.
- Partner-order and SMILES-serialization invariance pass.
- No recommendation, admission, or rule behavior changes.

### Phase 2: Aromatic vertical slice

Tasks:

1. Perceive connected aromatic ring systems and normalized anchor-relative
   positions.
2. Classify named ring families with generic fallbacks.
3. Record ring sizes, fusion, heteroatom elements, charges, and electronic roles.
4. Implement ortho occupancy and size-weighted ortho burden.
5. Implement auditable aromatic electronic-demand contributions.
6. Add reactive heteroatom lone-pair classification for aromatic N/O/S sites.
7. Add compact and expanded aromatic rendering.

Acceptance:

- Benzene, pyridine isomers, diazines, five-membered heteroaromatics, and fused
  systems are correctly distinguished.
- Ring numbering and symmetric serialization do not change the normalized
  descriptor.
- Two equally occupied ortho environments with different substituent sizes
  receive different burden classes.
- The concise output contains ring identity and high-value reactivity factors,
  but omits raw shell atom counts.

### Phase 3: Alkyl vertical slice

Tasks:

1. Classify reactive carbon substitution independently of partner role.
2. Calculate α and β branching, cyclicity, ring size, and β-H availability.
3. Detect benzylic, allylic, propargylic, and adjacent-heteroatom activation.
4. Implement bounded graph accessibility and branch contributions.
5. Emit elimination, strain, and radical-stabilization modifiers as observations.
6. Add compact and expanded alkyl rendering.

Acceptance:

- Methyl, primary, secondary, tertiary, neopentyl-like, cyclic, benzylic,
  allylic, and propargylic examples are separated.
- Primary amine identity remains independent of steric class on its attached
  alkyl group.
- β-H absence and presence are deterministic and atom-provenanced.
- Existing alkyl substitution grammars and product reconstructions retain parity.

### Phase 4: Remaining context variants

Implement in this order:

1. alkenyl;
2. alkynyl;
3. acyl;
4. sulfonyl and phosphoryl;
5. heteroatom-centered N/O/S/P;
6. generic fallback.

For each variant:

- add one typed context record;
- add one deterministic calculator;
- add one compact renderer;
- add one positive, negative, ambiguous, charged, and invariance test set;
- add explicit `not_applicable` and `unresolved` behavior;
- preserve legacy outputs until the migration gate.

Acceptance:

- All context IDs currently emitted by `context_facets.v2.json` have a typed
  profile or intentional fallback.
- No calculator uses reaction names or source labels as structural truth.
- Unrecognized chemistry remains serialized with evidence rather than guessed.

### Phase 5: Unified rendering and record serialization

Tasks:

1. Finalize one domain-owned compact renderer and one expanded evidence renderer.
2. Add `reactivity_profile` to `ReactionPartner` and generic canonical records.
3. Replace GUI, CLI, and concise-review chemistry formatting with the shared
   renderer.
4. Rename the user-facing section to `Reactivity profile`.
5. Keep legacy dictionaries in canonical JSON for compatibility, but stop
   displaying their low-value fields.
6. Bump review/export schema versions and document the added nested field.

Target compact examples:

```text
Electrophile Ar-Br — benzene; ortho burden high (2/2);
electron demand balanced

N partner AromN-H — carbazole; access open;
N lone pair pyrrole-like/low availability

Electrophile R-Br — secondary alkyl, beta-branched;
access hindered; beta-H present
```

Acceptance:

- GUI, CLI, and concise CSV use the same category names and field ordering.
- Application code contains no descriptor interpretation rules.
- Default output omits method names, raw scores, atom indices, and shell counts.
- Expanded output retains scores, contributors, confidence, and provenance.

### Phase 6: Condition-rule fact migration

Tasks:

1. Add stable profile facts to `PartnerRuleFacts`.
2. Add allowlisted predicates for context kind, ring family, ortho burden,
   accessibility, alkyl substitution, branching, β-H status, activation class,
   lone-pair availability, and selected modifiers.
3. Project facts only from structurally verified profiles.
4. Run old and new rule matchers in shadow comparison.
5. Convert rule definitions only after matches are explained and reviewed.
6. Bump the condition-rule schema when new predicates become active.

Acceptance:

- Every existing rule match has an explained parity result or approved chemistry
  correction.
- Missing profile facts cannot accidentally satisfy a rule.
- Low-confidence or unresolved facts do not trigger hard structural overrides.
- Rules remain declarative and contain no executable code references.

### Phase 7: Signature and retrieval activation

This is the intentional identity-breaking phase.

Tasks:

1. Define `environment_tokens_v2()` from stable binned profile fields.
2. Decide which tokens participate in L0 and L4; keep L1 handle identity
   environment-neutral.
3. Replace full open-dictionary hashing with canonical typed environment tokens.
4. Add the new descriptor definitions to the taxonomy identity manifest.
5. Bump the reaction-signature schema and signature-feature definition versions.
6. Generate new signature IDs and rebuild all generic indexes.
7. Add context-aware similarity:
   - categorical agreement for identities and bins;
   - normalized numerical distance within the same activation axis;
   - explicit missingness penalties;
   - no cross-axis numerical comparison.
8. Keep hard chemistry compatibility ahead of the new similarity score.
9. Record retrieval and score traces using readable feature names.

Acceptance:

- Signature IDs are deterministic and partner-order invariant.
- Index loaders reject incompatible legacy environment schemas with a clear
  migration error.
- Unknown-family mapped reactions remain retrievable through the generic ladder.
- Established-family baseline performance does not materially regress.
- Retrieval changes are attributable to specific context-aware features.

### Phase 8: Calibration, cleanup, and removal

Tasks:

1. Evaluate coverage, accuracy, retrieval-level distribution, and ranking changes.
2. Calibrate bins and weights only where outcome data provide adequate support.
3. Label uncalibrated values as chemistry priors.
4. Remove legacy environment reads from the GUI, CLI, conversion, rules,
   signatures, indexing, and similarity.
5. Remove `local_heavy_atoms_r2`, `qualitative_sum`, and duplicate legacy class
   paths from public contracts.
6. Remove compatibility adapters after all supported artifacts have migrated.
7. Update READMEs and the primary type-agnostic implementation document.

Acceptance:

- Repository search finds no active production consumer of removed legacy keys.
- Full tests pass with only the typed path enabled.
- Dataset snapshot changes include chemistry explanations.
- One canonical profile, rendering, rule-fact, and retrieval-feature path remains.

## 9. Testing plan

### 9.1 Unit tests

Add:

```text
tests/reactive_taxonomy/test_descriptor_models.py
tests/reactive_taxonomy/test_aromatic_descriptors.py
tests/reactive_taxonomy/test_alkyl_descriptors.py
tests/reactive_taxonomy/test_unsaturated_descriptors.py
tests/reactive_taxonomy/test_activated_center_descriptors.py
tests/reactive_taxonomy/test_heteroatom_descriptors.py
tests/reactive_taxonomy/test_reactivity_modifiers.py
tests/reactive_taxonomy/test_reactivity_rendering.py
tests/reactive_taxonomy/test_descriptor_definitions.py
```

Test:

- type and range validation;
- deterministic serialization;
- status semantics;
- contribution sorting;
- score bin boundaries;
- invalid definition rejection;
- SMARTS cache usage;
- omission and expanded-rendering rules.

### 9.2 Chemistry regression matrix

Include at minimum:

- benzene, pyridines, diazines, pyrrole, imidazole, indole, quinoline,
  naphthalene, and carbazole;
- zero, one, and two ortho substituents with small and bulky groups;
- methyl, primary, secondary, tertiary, neopentyl, cyclic, benzylic, allylic,
  and propargylic centers;
- terminal/internal alkenes and alkynes, E/Z, ring-constrained, and conjugated
  examples;
- acid chlorides, anhydrides, esters, amides, aldehydes, ketones, sulfonyl
  halides, and phosphoryl centers;
- primary/secondary amines, aniline, diarylamine, amide, sulfonamide, alcohol,
  phenol, thiol, and aromatic N-H;
- neutral, charged, protonated, deprotonated, and radical centers;
- leaving-group and transfer-carrier variants.

For every context include:

- positive detection;
- negative near-neighbor;
- ambiguous or unsupported structure;
- conflicting source-label evidence;
- explicit charge case;
- alternative valid SMILES;
- atom-map and reactant-order invariance.

### 9.3 Reaction and recommendation regression

Retain the required reaction-signature regressions:

- Suzuki;
- C-N, C-O, and C-S;
- mapped unknown-family reactions;
- invalid maps;
- deterministic IDs;
- partner-order invariance.

Add:

- aromatic versus alkyl versions of the same leaving group;
- benzylic versus unactivated alkyl;
- hindered versus open aryl partners;
- pyridine position isomers;
- N-H partners with different lone-pair availability;
- activated versus resonance-deactivated acyl centers.

Before Phase 7, assert no changes to:

- L0-L4 keys;
- `signature_id`;
- admission status;
- retrieval level and candidate counts;
- rule matches;
- recommendation order.

During Phase 7, approve changes only with per-record feature traces.

### 9.4 Performance tests

Measure:

- molecule featurization time by heavy-atom and ring count;
- aromatic-system classification time;
- descriptor-cache hit rate;
- dataset conversion throughput;
- index size and build time;
- query retrieval latency.

Targets:

- no network calls;
- no unbounded ring or subgraph search;
- SMARTS compiled once through the centralized cache;
- descriptor calculation linear or bounded-polynomial for supported graph sizes;
- no material regression in interactive query latency.

## 10. Migration gates and versioning

| Gate | Legacy dictionaries | Typed profile | Signature identity | Retrieval/rules |
|---|---|---|---|---|
| Baseline | Active | Absent | v2 unchanged | Unchanged |
| Shadow | Active | Serialized at site level | v2 unchanged | Legacy active |
| Display | Active but hidden | Partner and records active | v2 unchanged | Legacy active |
| Rule shadow | Active | Active | v2 unchanged | Old/new compared |
| Signature activation | Compatibility only | Active | Intentionally bumped | New features active |
| Cleanup | Removed | Sole path | New version required | Sole typed path |

No phase may delete a legacy field until all consumers in the repository have
migrated and the next artifact version is active.

## 11. Pull-request sequence

Keep changes reviewable with this order:

1. `test(reactive_taxonomy): freeze environment migration baseline`
2. `feat(reactive_taxonomy): add typed reactivity profile contracts`
3. `feat(reactive_taxonomy): add aromatic context descriptors`
4. `feat(reactive_taxonomy): add alkyl context descriptors`
5. `feat(reactive_taxonomy): add remaining context descriptors`
6. `feat(reactive_taxonomy): add shared reactivity renderer`
7. `refactor(condition-recommender): serialize typed reactivity profiles`
8. `refactor(app): use shared reactivity profile rendering`
9. `feat(condition-recommender): add typed rule facts in shadow mode`
10. `feat(condition-recommender): activate context-aware signature features`
11. `refactor: remove legacy environment descriptor path`

Each PR must state:

- chemistry motivation;
- changed contracts and schemas;
- definition-version changes;
- signature or index impact;
- migration/removal impact;
- tests and dataset evaluation results.

## 12. Definition of done

The work is complete when:

- every usable selected site has one typed context-aware reactivity profile;
- aromatic profiles capture ring family, sizes, fusion, and relative heteroatoms;
- alkyl profiles capture substitution, branching, activation, cyclicity, and
  β-hydrogen status;
- remaining supported contexts expose their highest-value reactivity factors;
- electronic axes are context-specific and never mislabeled as physical electron
  density;
- GUI, CLI, CSV, JSON, rules, signatures, and retrieval use the same canonical
  typed observations;
- hard chemistry filters still precede similarity;
- signatures and indexes are deterministic and version-compatible;
- ambiguity, confidence, contributors, and provenance are retained;
- all legacy descriptor consumers and duplicate rendering paths are removed;
- `pytest -q` passes.

## 13. Recommended first implementation slice

Start with Phases 0-2 only:

1. freeze baseline signatures and consumer behavior;
2. add typed common contracts and definition validation;
3. attach profiles in shadow mode;
4. implement the aromatic vertical slice;
5. render the proposed clean aromatic summary in a test-only or opt-in view;
6. prove that current signature IDs and recommendation behavior are unchanged.

This slice directly addresses the confusing `Ar`/`HeteroAr` output while
establishing the reusable contract needed for alkyl and other contexts. It does
not prematurely activate new retrieval or condition-rule behavior.
