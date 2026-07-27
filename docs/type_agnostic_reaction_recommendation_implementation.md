# Type-Agnostic Reaction Recommendation System

## System Design and Implementation Plan

**Status:** Proposed implementation plan  
**Scope:** `reactive_taxonomy`, `condition_registry`, and `condition_recommender`  
**Primary goal:** Recommend reaction conditions from reaction datasets even when a reaction has no reliable named type.

---

## 1. Purpose

Named reaction families such as Suzuki, Buchwald–Hartwig, and Ullmann are useful chemistry annotations. They provide strong retrieval filters and make explanations easier for chemists. They are not, however, a complete representation of a reaction:

- many useful dataset records have no reaction name;
- source reaction names may be missing, inconsistent, or incorrect;
- related transformations may use different mechanisms or catalyst systems;
- an unfamiliar transformation may still share reactive handles, bond changes, local environments, and condition requirements with known precedents.

The new system must therefore be **reaction-type assisted, not reaction-type dependent**.

Every valid reaction should receive a chemistry-based representation. A named family is added only when supported by evidence. Recommendation must remain possible when the family is unknown or ambiguous.

```text
Named family (optional):       Suzuki
Generic transformation:       sp2 C–X + C–B -> sp2 C–C
Reaction signature (required): bond edits + reactive sites + environments + spectators
```

“Type-agnostic” does not mean chemistry-agnostic. Retrieval must still enforce bond-change, reactive-handle, valence, functional-group, and condition-compatibility constraints.

---

## 2. Design Principles

### 2.1 The molecular graph is the source of truth

Reaction labels and family names are derived annotations. Machine identity must come from parsed structures, participating atoms, bond edits, and their environments.

### 2.2 Separate observation, interpretation, and recommendation

The system should have three explicit layers:

1. **Observation:** components, functional groups, reactive sites, atom correspondence, bond edits, local descriptors, and retained groups.
2. **Interpretation:** transformation class, role assignment, product connection, compatible named families, and confidence.
3. **Recommendation:** precedent admission, retrieval, compatibility filtering, scoring, recipe aggregation, and explanation.

A failure to assign a named family must not discard valid observations.

### 2.3 Confidence controls behavior

Family classification and role assignment must include evidence and confidence. Family-specific filters may be used at high confidence; lower-confidence records fall back to generic transformation and reaction-signature retrieval.

### 2.4 Exact evidence is preferred, uncertainty is preserved

The preferred evidence order is:

1. supplied valid atom mapping and observed bond edits;
2. exact product reconstruction from a taxonomy operator;
3. uniquely supported reactive-site grammar without exact reconstruction;
4. generic graph comparison or unresolved candidate set.

The first implementation must not silently invent atom correspondence or force a named family. Ambiguous candidates should be retained as alternatives or routed to review.

### 2.5 Definitions and code have different responsibilities

Versioned JSON definitions should hold declarative chemistry data: handles, grammar requirements, compatible family names, rendering rules, feature vocabularies, and retrieval weights. Python operators should implement graph edits, validation, descriptor calculation, and other logic that cannot be safely expressed as data.

---

## 3. Target Architecture

```text
reaction SMILES / dataset row
          |
          v
reactive_taxonomy
  parse components
  featurize molecules and sites
  determine observed/predicted bond edits
  build generic ReactionSignature
  optionally classify family
          |
          v
condition_registry
  resolve condition identities
  assign contextual roles
  construct canonical condition recipe
          |
          v
condition_recommender
  generic dataset conversion and admission
  hierarchical chemistry retrieval
  compatibility filtering
  similarity/ranking and recipe aggregation
          |
          v
recommendations + precedents + explanations + cautions
```

Package ownership must remain clean:

- `reactive_taxonomy` owns molecular and reaction chemistry features.
- `condition_registry` owns substance identity, aliases, roles, families, and condition normalization.
- `condition_recommender` owns dataset admission, indices, retrieval, scoring, and recommendations.
- None of these standalone packages should depend on legacy `chemtools` modules.

---

## 4. Canonical Reaction Representation

### 4.1 `ReactionEdit`

Replace family-specific assumptions about a single new connection with a generic edit model.

```python
@dataclass(frozen=True)
class ReactionAtomReference:
    side: str
    component_index: int
    atom_index: int
    atom_map_number: int | None
    element: str
    formal_charge: int
    aromatic: bool
    hybridization: str
    local_environment_id: str


@dataclass(frozen=True)
class ReactionEdit:
    edit_type: Literal["formed", "broken", "order_changed", "hydrogen_change"]
    atom_1: ReactionAtomReference
    atom_2: ReactionAtomReference | None
    old_order: str | None
    new_order: str | None
    evidence: str
    confidence: float
```

Hydrogen changes may initially be implicit when hydrogens are absent from reaction SMILES, but the schema must allow them for reductions, oxidations, eliminations, and proton-transfer-sensitive transformations.

### 4.2 `ReactionPartner`

Use neutral partner terminology instead of hard-coding `electrophile`, `transfer_partner`, or `nucleophile` into the base schema.

```python
@dataclass(frozen=True)
class ReactionPartner:
    partner_id: str
    component_index: int
    role: str | None
    role_confidence: float
    reactive_site_ids: tuple[str, ...]
    handle_tokens: tuple[str, ...]
    anchor_contexts: tuple[str, ...]
    chemist_label: str
    steric: dict[str, Any]
    electronic: dict[str, Any]
    nearby_groups: tuple[dict[str, Any], ...]
    spectator_group_ids: tuple[str, ...]
    flags: tuple[str, ...]
```

Family overlays may still expose convenient roles such as `electrophile`, `nucleophile`, `transfer_partner`, `acyl_donor`, or `hydrogen_acceptor`.

### 4.3 `ProductTransformation`

The current `ProductConnection` is valuable for single-bond-forming reactions and should be retained as a compatibility view. The general model must support multiple edits.

```python
@dataclass(frozen=True)
class ProductTransformation:
    edits: tuple[ReactionEdit, ...]
    formed_connection_labels: tuple[str, ...]
    concise_label: str | None
    exact_product_verified: bool
    evidence: str
```

This supports substitutions and couplings now, and later supports hydrogenation, oxidation, elimination, cyclization, rearrangement, and multi-bond transformations.

### 4.3.1 `ReactionEvent` and multi-event reactions

`ProductTransformation.edits` remains the aggregate observed transformation.
The edits are additionally partitioned into connected `ReactionEvent` objects:

```python
@dataclass(frozen=True)
class ReactionEvent:
    event_id: str
    event_signature_key: str
    edits: tuple[ReactionEdit, ...]
    partner_ids: tuple[str, ...]
    reactive_site_ids: tuple[str, ...]
    formed_bond_types: tuple[str, ...]
    broken_bond_types: tuple[str, ...]
    order_changes: tuple[str, ...]
    topology: ReactionTopology
    transformation_class: str | None
    named_family: str | None
    evidence: str
    confidence: float
```

Edits sharing an atom reference belong to the same event. Disconnected groups
are separate events, even when they occur on one substrate. Event order is
canonical and does not depend on reactant serialization. A signature with two
or more events has `event_scope="multi_event"`; equivalent event labels may
be grouped by multiplicity for display. Cross-event relations record shared
components without inferring an unsupported temporal sequence.

Balanced unmapped handle substitutions may use exact composite reconstruction:
operators must consume distinct reactive sites and distinct partner instances,
and the combined product must exactly match an observed product. Missing
stoichiometric partner copies are never synthesized in memory to force a match.

### 4.4 `ReactionTopology`

Intramolecularity is an orthogonal topology dimension, not a named reaction
family. It is derived from atom-provenanced formed bonds and observed reactant
and product graphs.

```python
@dataclass(frozen=True)
class ReactionTopology:
    reaction_scope: Literal[
        "intramolecular", "intermolecular", "mixed", "unimolecular", "unresolved"
    ]
    participating_component_indices: tuple[int, ...]
    role_component_indices: dict[str, int]
    same_component_role_groups: tuple[tuple[str, ...], ...]
    formed_bond_scopes: tuple[str, ...]
    reactant_tether_distances: tuple[int, ...]
    formed_ring_sizes: tuple[int, ...]
    ring_count_delta: int | None
    evidence: str
    confidence: float
```

Grammar definitions use typed role relationships with `same`, `different`, or
`same_or_different`. This lets one chemistry grammar and one operator represent
both intermolecular reactions and intramolecular closures. Unknown mapped
reactions receive topology directly from their observed edits.

### 4.5 `ReactionSignature`

This is the stable interface between featurization, dataset conversion, and recommendation.

```python
@dataclass(frozen=True)
class ReactionSignature:
    signature_id: str
    edit_signature: str
    formed_bond_types: tuple[str, ...]
    broken_bond_types: tuple[str, ...]
    order_changes: tuple[str, ...]
    events: tuple[ReactionEvent, ...]
    event_count: int
    event_scope: str
    event_relations: tuple[ReactionEventRelation, ...]
    partners: tuple[ReactionPartner, ...]
    product_transformation: ProductTransformation | None
    topology: ReactionTopology
    transformation_class: str | None
    transformation_confidence: float
    named_family: str | None
    family_confidence: float
    compatible_named_families: tuple[str, ...]
    spectator_groups: tuple[ReactionSpectatorGroup, ...]
    global_descriptors: dict[str, Any]
    warnings: tuple[str, ...]
    evidence_quality: str
    definition_versions: dict[str, str]
    schema_version: str
```

The `signature_id` must be deterministic and generated from normalized, versioned chemistry fields—not display labels or source reaction names.

### 4.6 Signature levels

Precompute several keys for hierarchical retrieval:

Reaction topology participates in L0–L2 so ring closures do not collide with
intermolecular assembly at chemistry-specific retrieval tiers. L3 deliberately
retains only normalized bond edits, providing a topology-agnostic fallback.
Fallback ranking must therefore score reaction scope, ring formation, and ring
size, and explanations must disclose any intra/intermolecular mismatch rather
than presenting the precedent as topology-equivalent.

- **L0 exact signature:** edit topology + handle tokens + anchor contexts + key local features.
- **L1 handle signature:** edit topology + handle families.
- **L2 transformation signature:** formed/broken bond types + transformation class.
- **L3 bond-edit signature:** formed/broken/order-changed atom-pair types.
- **L4 environment signature:** local steric/electronic/group fingerprints without a family requirement.

These keys make fallback behavior explicit and auditable.

---

## 5. Optional Reaction-Family Layer

### 5.1 `ReactionFamilySpec`

Introduce a registry contract so adding a family does not require scattered conditionals.

```python
@dataclass(frozen=True)
class ReactionFamilySpec:
    family_id: str
    display_name: str
    transformation_classes: tuple[str, ...]
    grammar_ids: tuple[str, ...]
    required_edit_patterns: tuple[str, ...]
    forbidden_edit_patterns: tuple[str, ...]
    operator_id: str | None
    environment_builder_id: str | None
    admission_policy_id: str | None
    condition_regime_ids: tuple[str, ...]
    definition_version: str
```

The registry should combine versioned definitions with explicit Python plugin tables for graph operators and environment builders. Do not dynamically import code specified by arbitrary JSON.

### 5.2 Family classification output

Family classification should return ranked evidence rather than one unconditional name:

```text
family_candidates:
  - family_id: Suzuki
    confidence: 0.98
    evidence: exact grammar and product reconstruction
  - family_id: generic_sp2_c_c_coupling
    confidence: 1.00
    evidence: observed C-C formation and compatible handles
```

Recommended policy:

- `>= 0.90`: family may be used as a strict first retrieval tier;
- `0.60-0.90`: family may contribute to scoring but not exclude generic precedents;
- `< 0.60`: store as a candidate only and retrieve type-agnostically.

Thresholds must be configuration values and validated on held-out data.

---

## 6. Featurization Workflow

### 6.1 Input parsing

- Parse `reactants>agents>products` and `reactants>>products`.
- Preserve source SMILES and atom-map numbers.
- Canonicalize a separate unmapped representation.
- Identify multiple products and invalid components explicitly.

### 6.2 Atom correspondence and edit extraction

Use this order:

1. Validate and use supplied atom maps when present.
2. Run taxonomy grammar candidates and exact product reconstruction.
3. Reconcile mapped edits with operator-predicted edits when both exist.
4. When neither source is usable, allow a narrow conserved-scaffold
   correspondence fallback: exactly one substantial reactant scaffold, exactly
   one product, every product heavy atom accounted for, and all best mappings
   producing the same normalized edit set.
5. Preserve ambiguous, unresolved, or conflicting evidence for review.

The initial release should borrow deterministic mapping validation and graph-comparison ideas from the old system where useful, but it must not depend on RXNMapper or legacy modules. A future optional mapper may implement a narrow interface and record its provenance.

The conserved-scaffold fallback is not a general reaction mapper. It rejects
multi-substrate assembly, insufficient scaffold conservation, candidate-limit
overflow, and chemically distinct minimal correspondences. Symmetry-equivalent
atom assignments may be accepted only when their normalized edit sets agree.

Every valid parsed reaction also receives a typed product-atom completeness
assessment. It records heavy-atom and element counts, product and reactant
mapping coverage, product-heavy-atom coverage, side-specific map inconsistencies,
and suspected missing reactants or insufficient partner multiplicity. Reactant
atoms absent from the reported main product are retained without rejection
because byproducts are commonly omitted. Definite product-heavy-atom excess
blocks signature generation and indexing; unresolved provenance and
inconsistent product mapping are retained for review. Missing reactants and
stoichiometric partner copies must never be synthesized.

### 6.3 Reactive-site and environment features

For every atom participating in an edit, capture:

- element, charge, aromaticity, hybridization, ring membership, and hydrogen count;
- reactive handle and leaving/transfer group;
- attachment-side context such as `Ar`, `HeteroAr`, `Alkenyl`, `Alkyl`, or `Acyl`;
- local steric class and numerical descriptors;
- local electronic class and interpretable contributing groups;
- nearby functional groups by graph distance;
- coordinating, acidic, basic, oxidizable, reducible, or catalyst-poisoning flags;
- alternative/competing reactive sites.

Descriptor calculations must be mechanism-neutral at the base layer. Family-specific interpretations belong in overlays.

### 6.4 Spectator groups

Spectators must be derived relative to the selected or observed edits, not treated as a fixed whole-molecule list. Store atom provenance, distance to each reaction center, and whether the group is unchanged, potentially coordinating, condition-sensitive, or competing.

### 6.5 Chemist-facing labels

Labels are rendered after the structured transformation is built. Selection
uses an evidence-ordered ladder:

1. expose conflicts when supplied mapping and grammar reconstruction disagree;
2. use an exact grammar label only when reconstruction and normalized edits
   agree;
3. otherwise match a versioned generic edit pattern;
4. when one validated bond-order event has sufficient local context, render a
   contextual before/after overlay while retaining the generic pattern label;
5. otherwise render the normalized edits literally;
6. use a reactant-only grammar label only when no verified edit evidence exists;
7. report the transformation as unresolved when neither edits nor exact
   reconstruction are available.

Generic label patterns are declarative and identify reusable edit combinations,
not named reaction families. A pattern such as C–X cleavage plus C–N formation
may support `C–N substitution`, but it does not distinguish Buchwald–Hartwig,
Ullmann, or SNAr chemistry. Named-family labels remain optional overlays.

Examples:

```text
Ar-Cl + Ar-B(OH)2 -> Ar-Ar
Alkenyl-Cl + Ar-B(OH)2 -> Ar-Alkenyl
HeteroAr-Br + R-NH2 -> HeteroAr-NH-R
unknown family: C-N substitution
unknown family: Ar-CH=CH2 -> Ar-CH2-CH3
reductive amination: HeteroAr-CH=O + Ar-NH2 -> HeteroAr-CH2-NH-Ar
unmatched edit pattern: C-O bond cleavage; C-N bond formation
```

The initial reductive-amination implementation follows the same evidence
ladder. A declarative grammar identifies an available aldehyde, ketone, or
formaldehyde carbonyl and a free primary or secondary amine. A registered
operator predicts carbonyl-oxygen removal, C–N formation, carbon hydrogen gain,
and nitrogen hydrogen loss. `reductive_amination` is assigned only after exact
product reconstruction; a reactant-only candidate is retained without forcing
the family when the product disagrees or does not resolve the assignment.

Contextual text is a display overlay, not reaction identity. Store the generic
interpretation (`C=C hydrogenation`), literal structural edits, contextual
reactant motif, and contextual product motif separately. Rendering definitions
and display style must not affect `signature_id`.

When the product connection is not verified, leave the product label empty:

```text
Ar-X + N-H ->
```

Do not present a predicted product label as observed fact.

---

## 7. Unified Dataset Conversion

### 7.1 Replace family-specific conversion engines

The current Suzuki, C-N, C-O, and C-S converters duplicate parsing, admission, flattening, and report generation. Replace them incrementally with:

```text
condition_recommender/conversion/
  engine.py              generic streaming conversion
  input_schema.py        column aliases and source adapters
  admission.py           common and plugin policies
  flatten.py             stable CSV/Parquet views
  reports.py             coverage and rejection reports
  policies/
    suzuki.py
    c_n_coupling.py
    acyl_transfer.py
```

Source adapters should only map dataset columns into a common raw record. They must not perform chemistry classification.

### 7.2 Converted record schema

Evolve `RecommendationRecord` to contain:

- raw source identity and row number;
- reaction SMILES;
- complete `ReactionSignature` or a versioned serialized form;
- optional named family and confidence;
- transformation class and confidence;
- admission tier and reasons;
- normalized condition recipe;
- yield, temperature, time, concentration, atmosphere, and other available outcomes;
- raw source fields for traceability;
- taxonomy, registry, and converter definition versions.

Keep nested JSON or Parquet as the canonical data artifact. CSV should be a review/export view because it cannot naturally represent multiple edits, partners, roles, or family candidates.

Persisted recommendation indices must record and validate the reaction-signature
schema, taxonomy definition versions, converted-record schema, and converter
version. A stale or mixed artifact must fail with a regeneration instruction;
it must not silently miss L0-L2 keys and degrade to a generic fallback.

### 7.3 Admission policy

Common admission gates:

- parseable reaction and products;
- valid outcome range;
- usable condition identities;
- sufficient edit or grammar evidence;
- no unresolved contradiction between mapping and reconstruction;
- provenance and schema versions present.

Suggested tiers:

- **verified:** exact or validated edit evidence, usable signature, and usable conditions;
- **review:** chemically plausible but ambiguous family/site/product or incomplete condition roles;
- **rejected:** invalid structure/outcome, no usable transformation evidence, or irreconcilable record.

Unknown `named_family` is not itself a reason for review or rejection.

Family plugins may add stricter gates for family-specific indices, while the same record can remain valid in the generic index.

---

## 8. Condition Registry Integration

Recommendation requires a structured recipe rather than source column labels alone.

### 8.1 Contextual role resolution

A substance can have multiple possible roles. Resolve its role using:

- canonical identity;
- source field;
- reaction signature and high-confidence family, if available;
- co-occurring condition components;
- curated role/family definitions.

For example, triethylamine may be a base, reagent, or additive; DMF may be a solvent or a chemically participating reagent. Preserve both the resolved role and the evidence.

### 8.2 `ResolvedConditionRecipe`

Add a versioned recipe with components such as:

```text
catalysts[]
ligands[]
bases[]
oxidants[]
reductants[]
activators[]
additives[]
solvents[]
other_components[]
temperature / time / concentration / atmosphere
```

Each component should retain raw identifier, canonical substance ID, resolved roles, role confidence, amount when available, and provenance.

Do not block generic recommendation while registry curation is incomplete. Unresolved substances remain searchable by normalized identity and receive an uncertainty flag.

---

## 9. Type-Agnostic Retrieval

### 9.1 Hierarchical retrieval ladder

Use the narrowest pool with adequate support:

1. high-confidence same family + exact partner signature;
2. same generic transformation + compatible handles and contexts;
3. same edit signature + compatible reactive-site environments;
4. same formed bond/order change + relaxed handles;
5. environment-neighbor retrieval across transformation labels;
6. broad condition-regime prior, returned with a strong caution.

Family is a retrieval feature and optional filter, not a mandatory partition key.

### 9.2 Hard compatibility before similarity

Before ranking, exclude or strongly penalize precedents with incompatible chemistry, for example:

- wrong net bond edit;
- incompatible substrate valence or reactive-handle class;
- condition-sensitive group conflicts;
- catalyst poison or strong coordination conflicts;
- oxidant/reductant incompatibility;
- missing mandatory regime components;
- incompatible solvent, temperature, or atmosphere constraints when known.

Hard constraints should be declarative and versioned where possible.

### 9.3 Similarity features

Score remaining precedents using interpretable components:

- edit topology similarity;
- reactive-handle compatibility;
- atom and attachment context;
- local steric similarity;
- local electronic similarity;
- nearby functional-group similarity by distance;
- spectator/compatibility similarity;
- condition-regime compatibility;
- family agreement weighted by family confidence;
- dataset quality and evidence quality.

Weights belong in versioned recommender definitions. Missing features must be handled explicitly rather than treated as matches.

### 9.4 Recipe aggregation

Group retrieved precedents by canonical `ResolvedConditionRecipe`, not raw source strings. Rank recipes using:

- weighted precedent similarity;
- support count and dataset diversity;
- shrinkage-adjusted yield or success rate;
- evidence and registry confidence;
- robustness across substrate environments;
- penalties for conflicts and missing recipe components.

Return both the aggregated recipe and representative precedents.

### 9.5 Explanation contract

Every recommendation should report:

- retrieval level used;
- whether a named family constrained retrieval;
- matching bond edits and handles;
- important environment matches and mismatches;
- relevant spectator/compatibility warnings;
- support and outcome summary;
- condition identity or role uncertainty;
- precedent reaction IDs.

---

## 10. Implementation Phases

### Phase A: Generic signature foundation

1. Add `ReactionAtomReference`, `ReactionEdit`, `ReactionPartner`,
   `ProductTransformation`, `ReactionTopology`, and `ReactionSignature` models.
2. Normalize the current `BondChange`, `ProductConnection`, family environment, and spectator output into this schema.
3. Generate deterministic L0-L4 signature keys.
4. Add family confidence/evidence without changing current family results.
5. Expose `reaction_signature` from `featurize_reaction()` while retaining current fields temporarily.

**Acceptance:** Existing Suzuki/C-N/C-O/C-S exact cases produce stable signatures; an unknown-family mapped reaction can produce a signature with `named_family=None`.

### Phase B: Family registry and transformation plugins

1. Introduce `ReactionFamilySpec` and a validated family registry.
2. Register existing grammars, operators, environment builders, and product renderers.
3. Remove scattered family conditionals after parity tests pass.
4. Add definition-version hashes to every analysis.

**Acceptance:** Existing reaction labels and verified product connections remain unchanged for the current test corpus.

### Phase C: Unified converter

1. Implement the common input record and source-column adapters.
2. Implement generic admission and reporting.
3. Route the four existing converter CLIs through the common engine.
4. Preserve old CSV columns as review views where useful, but make nested records canonical.
5. Add an auto/mixed dataset mode that does not require a declared family.

**Acceptance:** Pilot counts and labels are explainably reconciled with current outputs; mixed datasets convert in one run; unknown-family verified signatures are admitted.

### Phase D: Resolved condition recipe

1. Add contextual role resolution to `condition_registry`.
2. Add `ResolvedConditionRecipe` and canonical recipe IDs.
3. Reprocess Suzuki, C-N, C-O, C-S, and amide pilots.
4. Report identity coverage, role confidence, ambiguity, and missing substances.

**Acceptance:** Raw identifiers remain traceable; multi-role substances are not collapsed incorrectly; recipes group equivalent source records.

### Phase E: Generic indexing and retrieval

1. Index signature levels, family candidates, environments, spectators, outcomes, and recipe IDs.
2. Generalize retrieval code away from Suzuki-specific `electrophile`/`transfer_partner` assumptions.
3. Implement the fallback ladder and minimum-support rules.
4. Add hard chemistry compatibility filters.
5. Add confidence-weighted family matching.

**Acceptance:** Queries with `named_family=None` return chemically related precedents; named Suzuki queries retain or improve baseline quality.

### Phase F: Evaluation and calibration

1. Split by reaction or substrate scaffold to prevent near-duplicate leakage.
2. Hide conditions from held-out reactions and measure recovery of observed recipes/classes.
3. Compare family-only, generic-only, and hybrid retrieval.
4. Measure coverage, top-k recipe recall, yield-weighted regret, compatibility violations, and explanation fidelity.
5. Calibrate family confidence and retrieval thresholds.

**Acceptance:** Hybrid retrieval improves coverage without increasing hard chemistry violations and does not materially degrade established-family performance.

### Phase G: Reaction-family expansion

Add families in increasing structural complexity:

1. amide, ester, thioester, and sulfonamide formation;
2. generic acyl substitution and nucleophilic substitution;
3. reductions, oxidations, and hydrogenations;
4. reductive amination, elimination, and olefination;
5. cyclizations and multi-edit transformations;
6. rearrangements, cascades, and pericyclic reactions only after multi-edit evidence is robust.

Each family must add a test pack, not a separate conversion architecture.

Implemented high-ROI foundation:

- acyl N/O/S formation and reductive amination use dedicated verified graph
  operators;
- alkyl C–N/O/S substitution reuses one handle-replacement operator and does
  not infer SN1 versus SN2 from connectivity alone;
- terminal-alkene Heck coupling uses endpoint hydrogen counts for deterministic
  regioselection and requires exact product reconstruction before assigning
  `heck`;
- Friedel–Crafts acylation combines an activated acyl electrophile with an
  aromatic C–H site and resolves regioisomers by exact reconstruction;
- unspecified alkene stereochemistry is never promoted to an E/Z claim;
- the 450-row sample evaluation now has 244 exact-product labels, including the
  Friedel–Crafts acylation example, six alkyl substitutions, and four exactly
  reconstructed Heck reactions.

The remaining Wittig, aldol, Minisci, Grignard, cycloaddition, and protection
gaps should be added only after their missing partner-site or multi-edit
contracts are defined. Source reaction names are insufficient evidence for
those grammars.

---

## 11. Proposed File Changes

```text
reactive_taxonomy/
  reaction_models.py                 extend generic models
  reaction_signatures.py             signature construction and keys
  reaction_edits.py                  observed/predicted edit normalization
  families/
    registry.py                      ReactionFamilySpec loader/validation
    plugins.py                       safe Python operator/builder tables
  definitions/
    reaction_families.v1.json
    signature_features.v1.json

condition_registry/
  recipes.py                         ResolvedConditionRecipe
  contextual_roles.py                role resolution with evidence
  definitions/
    role_resolution.v1.json

condition_recommender/
  conversion/
    engine.py
    input_schema.py
    admission.py
    flatten.py
    reports.py
    policies/
  generic_indexing.py
  generic_retrieval.py
  compatibility.py
  similarity.py
  recipe_ranking.py
  definitions/
    generic_retrieval.v1.json
    compatibility.v1.json
```

Names may be adjusted during implementation, but package ownership should not change.

---

## 12. Testing Strategy

### 12.1 Unit tests

- edit normalization from mapped reactions;
- edit normalization from exact reconstruction;
- deterministic signature IDs;
- partner ordering invariance;
- atom-index and component provenance;
- multiple formed/broken bonds;
- family classification confidence;
- spectator distance and retention;
- contextual condition roles;
- retrieval fallback selection;
- missing-feature scoring;
- compatibility exclusions.

### 12.2 Chemistry regression tests

Maintain curated positive, negative, and ambiguous examples for:

- Suzuki C-C;
- sp2 C-N, C-O, and C-S;
- amide and other acyl transfers;
- reactions with competing handles;
- reactions with multiple products;
- mapped and unmapped equivalents;
- unknown-family but valid bond-edit cases;
- invalid atom maps and product mismatches.

### 12.3 Dataset tests

For every pilot dataset, record:

- parsed and signature coverage;
- exact/reconstructed/mapped evidence counts;
- named-family coverage and confidence distribution;
- verified/review/rejected counts and reasons;
- condition identity and role coverage;
- retrieval pool sizes and fallback levels;
- top-k recommendation evaluation.

Dataset snapshot changes must be reviewed as chemistry changes, not merely accepted as code changes.

### 12.4 Performance tests

- lazy SMARTS compilation through `reactive_taxonomy.chemistry.smarts_cache`;
- streaming conversion without loading large datasets entirely into memory;
- cached signature and registry lookup;
- index build and query latency targets defined against the larger private datasets.

Never call `Chem.MolFromSmarts()` directly in the standalone system.

---

## 13. Migration and Compatibility

Backward compatibility with legacy `chemtools` is not a goal. During migration inside the new system:

1. add generic signature fields while keeping current outputs;
2. validate parity on current pilots;
3. route old family CLIs through the generic converter;
4. migrate retrieval to the generic index;
5. remove duplicated family converters and obsolete retrieval fields;
6. bump schema versions and document the final canonical format.

Do not maintain two permanent recommendation paths. Temporary adapters must have removal criteria and tests.

---

## 14. Decisions and Non-Goals

### Decisions

- Named reaction family is optional metadata with confidence.
- `ReactionSignature` is the primary recommendation representation.
- Unknown-family reactions may be verified and admitted.
- Exact product evidence remains the preferred confirmation boundary.
- Generic retrieval uses chemistry constraints before similarity.
- Family-specific logic is implemented as registered overlays/plugins.
- Canonical converted artifacts are nested/versioned; CSV is primarily for review.

### Initial non-goals

- automatically naming every organic reaction;
- inferring mechanisms from conditions alone;
- unrestricted atom mapping or MCS-based product assignment;
- replacing deterministic chemistry with an opaque embedding model;
- supporting cascades and rearrangements before multi-edit validation exists;
- cleaning the entire condition registry as a prerequisite for implementation.

---

## 15. Definition of Done

The first complete type-agnostic release is done when:

1. `featurize_reaction()` returns a versioned `ReactionSignature` for supported exact reactions and usable mapped unknown-family reactions.
2. A mixed dataset can be converted without selecting a family-specific converter.
3. Missing `named_family` does not cause automatic rejection.
4. Conditions are represented by canonical resolved recipes with provenance and uncertainty.
5. Retrieval follows a documented family-to-generic fallback ladder.
6. Compatibility rules run before similarity ranking.
7. Recommendations cite precedents and explain chemistry matches, mismatches, fallback level, and uncertainty.
8. Current Suzuki/C-N/C-O/C-S behavior passes regression tests.
9. Held-out evaluation shows the hybrid system expands coverage without increasing hard chemistry incompatibilities.
10. Family-specific duplicate conversion paths are removed.

---

## 16. Recommended Immediate Work

Start with Phase A only. It provides the common contract needed by every later phase and can be validated without changing recommendation behavior.

The first implementation slice should:

1. add the five generic dataclasses;
2. construct them from the current selected candidate, mapped bond changes, family environment, product connection, and spectators;
3. attach the signature to `ReactionAnalysis`;
4. add deterministic signature-key tests for Suzuki, C-N, C-O, C-S, and one unknown-family mapped reaction;
5. serialize the new signature into pilot review output without yet changing admission or retrieval.

This creates a safe comparison point before converters and retrieval are consolidated.
