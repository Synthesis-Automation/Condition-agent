# Reaction Featurization Workflow

The reaction featurization system is graph-first and label-last. A reaction
name or source label is never the primary classifier.

```text
Reaction SMILES
  -> parse molecules
  -> detect functional groups and reactive sites
  -> collect mapping, reconstruction, and correspondence evidence
  -> reconcile providers and retain distinct edit hypotheses
  -> select verified graph edits when the evidence is unique
  -> minimize mapped edits into atom transitions and remote subgraphs
  -> assign transformation and optional family interpretations
  -> build topology, environments, spectators, and events
  -> generate display labels
  -> generate a deterministic reaction signature
  -> serialize into recommendation records
```

## 1. Parse the reaction

The main entry point is
[`featurize_reaction()`](../reactive_taxonomy/reaction_api.py).

It accepts either:

```text
reactants > agents > products
```

or:

```text
reactants >> products
```

[`parse_reaction_smiles()`](../reactive_taxonomy/reaction_parser.py) splits each
side into components and preserves:

- component index;
- original and canonical SMILES;
- reactant, agent, or product side;
- supplied atom-map numbers;
- parsing errors and warnings.

Agents are analyzed, but they are not used as reaction partners during
structure-evidence collection or grammar candidate enumeration.

Parsing validates syntax and component SMILES; it does not require a globally
balanced chemical equation. Reaction records commonly report only the main
product and omit salts, leaving-group products, and other byproducts. Product
atom provenance is therefore evaluated later by a dedicated
`ReactionCompletenessAssessment` instead of by naive formula equality at this
stage.

## 2. Featurize every molecule

Each component passes through
[`featurize_molecule()`](../reactive_taxonomy/api.py).

For each molecular graph, the system detects and emits:

- functional groups;
- reactive sites and handles;
- atom roles within each site;
- site availability, such as available or blocked;
- aromatic, alkyl, carbonyl, heteroatom, and other local contexts;
- local steric and electronic environments;
- chemist-readable site labels.
- canonical reactive links, bond capacities, and connection endpoints.
- explicit C/N/O/S anion endpoints, polarized C=N capacities, strained
  epoxide/aziridine release links, silyl-ether O–Si links, and selected
  Li/Cu/Al transfer links.

Detection is intentionally conservative. It does not infer anions from neutral
precursors, treat nitro or azide resonance atoms as free anions, classify
cyclopropane as a strained heterocycle, or treat DIBAL alkyl groups as
carbon-transfer partners.

For example, the molecular observations for an aryl bromide and a boronic acid
may include:

```text
aryl bromide:
  site_type = leaving_group
  anchor_context = Ar
  handle_token = Br

boronic acid:
  site_type = transfer_group
  anchor_context = Ar
  handle_token = B(OH)2
```

These are molecular observations, not yet a reaction-family assignment.

## 3. Enumerate reaction candidates

[`enumerate_reaction_candidates()`](../reactive_taxonomy/reaction_candidates.py)
compares canonical connectivity sites with the declarative grammars in
[`reaction_grammars.v2.json`](../reactive_taxonomy/definitions/reaction_grammars.v2.json).

Each grammar specifies:

- required partner roles;
- site-type constraints;
- allowed handle and context combinations;
- same-component or different-component relationships;
- transformation class;
- optional compatible named families.

For an aryl bromide and aryl boronic acid, a candidate may contain:

```text
grammar_id: suzuki_miyaura
transformation_class: c_c_transfer_coupling
roles:
  electrophile -> aryl bromide site
  transfer_partner -> aryl boronic acid site
compatible_named_families:
  suzuki_miyaura
```

At this stage the candidate is only structurally plausible.

## 4. Apply the connectivity rewrite and verify the product

Each candidate is passed to
[`apply_connectivity_rewrite()`](../reactive_taxonomy/connectivity_rewrite.py).

The bounded rewrite:

- removes or changes required bonds;
- forms new bonds;
- records hydrogen changes;
- produces a predicted product;
- returns predicted bond edits.

The rewrite registry uses four reusable shapes: release-and-connect,
split-and-distribute, depart-and-unsaturate, and localized bond-state change.
Split-and-distribute accepts explicit A–B donors and implicit A–H donors
through one normalized contract and may return multiple constitutional
outcomes. Product reconstruction selects a unique orientation; reactants alone
do not force regio- or stereochemistry.

The predicted product is canonicalized without atom maps and compared with the
observed product. Candidate verification becomes one of:

- `exact_product_reconstruction`;
- `product_mismatch`;
- `construction_failed`.

A candidate is selected only when the structural evidence is sufficiently
unique. Symmetry-equivalent candidates can collapse to one interpretation;
chemically different exact candidates remain ambiguous.

Multiple rewrites may also be composed for multi-event reactions.

## 5. Extract and reconcile observed reaction edits

[`resolve_reaction_evidence()`](../reactive_taxonomy/reaction_edits.py)
establishes the actual transformation evidence.

The effective evidence priority is:

1. Valid supplied atom mapping.
2. Exact single-event rewrite reconstruction.
3. Exact multi-event rewrite reconstruction.
4. Conservative unique-scaffold correspondence.
5. Bounded global multi-reactant correspondence.
6. Bounded fragmented-scaffold correspondence for topology-changing,
   single-substrate reactions.
7. Unresolved, ambiguous, or conflicting evidence.

The resolver exposes each attempted source as a typed
`ReactionEvidenceCandidate`. Mapping, exact reconstruction, and correspondence
are providers of structural evidence; they are not separate recommendation
pipelines. The resolver applies one precedence and contradiction policy after
the providers run. `normalize_reaction_edits()` remains a compatibility alias.

An optional RXNMapper provider is available for offline benchmarks and explicit
converter/query-time review mode:

```powershell
python -m pip install -r requirements-mapping.txt
python benchmarks/fischer_rxnmapper_poc.py --progress

python -m condition_recommender.generic_conversion_cli `
  data-processor/reaction_dataset/Fischer_indole_synthesis.csv `
  results/fischer_indole_rxnmapper_conversion `
  --use-rxnmapper

python -m condition_recommender.generic_recommend_cli `
  "<reaction_smiles>" `
  --records "<generic_index.json>" `
  --use-rxnmapper
```

The adapter is implemented in
[`external_atom_mapping.py`](../reactive_taxonomy/external_atom_mapping.py).
It verifies that mapper serialization preserved the input chemistry, projects
only retained-to-unreported attachment boundaries, and reuses the typed mapped
edit normalizer. Generated edits use `external_atom_mapping` evidence and carry
the mapper confidence. They are not silently treated as supplied maps.

The integration first runs ordinary internal analysis. It skips the mapper for
reactions with supplied maps or an existing resolved signature. For an
unresolved or ambiguous reaction it records one of these dispositions:

- `external_mapping_internal_consensus`: exactly one internal edit hypothesis
  matches the external normalized edit profile;
- `external_mapping_only`: the mapper supplies a valid edit profile where no
  internal hypotheses exist;
- `external_mapping_hypothesis_conflict` or
  `external_mapping_ambiguous_hypothesis_match`: retain the original
  hypotheses and route to review;
- `external_mapping_signature_unavailable` or `external_mapping_failed`: retain
  the base analysis and failure provenance.

Converter records serialize this assessment as `external_atom_mapping`.
Mapper-derived signatures are always `review_only` and excluded from the
default index. At query time the first two dispositions may use the normal
signature retrieval ladder against already-verified precedents, with mapper
status, provider, confidence, and mandatory expert-review cautions in the
result.

### Template-free reaction minimization

When normalized edits retain atom correspondence across both sides,
[`build_reaction_core_projection()`](../reactive_taxonomy/reaction_core.py)
builds a `ReactionCoreProjection` directly from the molecular graphs and edits.
It does not load a reaction grammar, template, source label, or reaction name.

The projection:

1. keeps every edit-participating atom as a before/after atom transition;
2. groups connected edits into events;
3. selects a smaller set of primary centers for a concise explanation;
4. removes active atoms from each molecular graph;
5. records each remaining connected component once as a remote subgraph;
6. records every cut connection as a typed attachment port; and
7. matches remote subgraphs across sides as retained, departing, appearing,
   changed, or unresolved.

Classes such as `aryl`, `heteroaryl`, and `alkyl` come from the removed graph:
aromaticity, ring composition, element composition, unsaturation, and
attachment topology. They are not a fixed grammar vocabulary applied before
the cut. Exact fragment SMILES, functional groups, attachment atoms, and bond
orders remain in the result even when the label displays `Ar`, `HetAr`, or
`R`.

For the mapped acetal example, the concise minimized transition is:

```text
C(H)(Ar)(=O) -> C(H)(Ar)(O-R)2
```

The projection emits four keys:

- `RCX2`: exact edit, atom-state, and fragment identity;
- `RCT2`: typed remote-subgraph identity;
- `RSH2`: mapping-robust retrieval shape containing the generic center
  transition, participant handles/sites, retained remote shape, and event
  count;
- `RCS2`: diagnostic center transition only.

Only `RSH2` is retrieval-eligible. `RCS2` is intentionally too broad:
reactions can share a carbon-center transition while requiring different
partner chemistry and conditions. The projection is a sibling observation to
RS3 and never upgrades ambiguous or externally mapped evidence.

The three desktop applications expose a `Use RXNMapper` checkbox that is
checked by default:

```powershell
python -m app.featurizer_gui
python -m app.reaction_converter_gui
python -m app.reaction_recommender_gui
```

The featurizer displays whether mapping was skipped, reached consensus, was the
only edit source, conflicted, or failed. The converter includes the mapper/model
identity in restartable shard identity and uses one worker to avoid loading
multiple model copies. The recommender caches mapper-enabled and internal-only
instances separately. Clearing a checkbox restores the internal-only path.

The global fallback is used only for unmapped, single-main-product reactions
after exact reconstruction and the narrower scaffold fallback fail. It matches
conserved subgraphs from a bounded number of reactant components into
non-overlapping product atoms, requires every product heavy atom to be
accounted for, and accepts only minimum-edit alternatives that imply the same
normalized chemistry. Candidate overflow, product element excess, additional
substantial products, or chemically different best mappings remain unresolved.
Explicit atom and E/Z stereochemistry is compared across each candidate
correspondence; alternatives that imply different stereo observations remain
ambiguous.
This recovers general assemblies such as additions, condensations, and
cycloadditions without assigning a named reaction or inventing a reactant.

When minimum-cost correspondences imply different normalized chemistry, the
resolver does not discard the alternatives. It emits deterministic
`ReactionEditHypothesis` objects containing:

```text
hypothesis_id (REH1)
provider and evidence
confidence
atom-provenanced edits and stereo observations
number of correspondences collapsed into this chemistry alternative
edit cost
topology
warnings
```

These hypotheses are review and query evidence, not observed edits. They do not
create a reaction signature or make a dataset record index-eligible.

Normalized edit types include:

- `formed`;
- `broken`;
- `order_changed`;
- `hydrogen_change`.

Phase 1 also dual-writes an internal evidence-scoped connectivity observation
before those compatibility edit types are consumed. `BondTransition`
distinguishes definite bond/no-bond states from an endpoint absent from the
reported main product or an unresolved state. Aggregated `HydrogenDelta` and
observed formal-charge `AtomStateTransition` objects join those transitions in
a canonical `ConnectivityEditGraph`. Its `CEG1` key is currently shadow
evaluation output only: it does not alter `ReactionSignature.signature_id`,
serialized reaction analyses, admission, retrieval, or recommendation
behavior.

Each `ReactionEdit` retains atom-level provenance:

```text
side
component_index
atom_index
atom_map_number
element
formal_charge
aromaticity
hybridization
local_environment_id
chiral_tag
CIP_code
evidence
confidence
```

Explicit atom and bond stereochemical descriptors are stored separately as
typed `ReactionStereoChange` observations. They distinguish retained, created,
destroyed, and descriptor-changed stereo without claiming a mechanism. If an
operator reconstructs an explicit stereoisomer opposite to the reported
product, the observed correspondence is retained with:

```text
conflicting_stereochemical_evidence
STEREOCHEMICAL_RECONSTRUCTION_CONFLICT
```

This conflict remains review-only and cannot enter recommendation retrieval.

When atom mapping and operator reconstruction agree, the evidence can become:

```text
validated_mapping_and_exact_reconstruction
```

When they disagree, the system retains the mapped edits but marks:

```text
conflicting_edit_evidence
MAPPING_RECONSTRUCTION_CONFLICT
```

It does not silently let a family interpretation override the observed graph.

### Product-atom completeness

After edit normalization, the system evaluates whether every reported product
heavy atom can be accounted for by the supplied reactants. The assessment
records:

- reactant and product heavy-atom counts;
- element counts and element excess on each side;
- reactant and product atom-mapping coverage;
- shared and side-specific atom-map numbers;
- estimated product-heavy-atom coverage;
- suspected missing reactants;
- suspected insufficient reactant multiplicity;
- evidence and warning codes.

Possible statuses are:

- `verified`: exact reconstruction, conservative correspondence, or complete
  product mapping establishes product-atom provenance;
- `incomplete`: the product contains more heavy atoms of one or more elements
  than all supplied reactants can provide;
- `unresolved`: the record is not demonstrably incomplete, but provenance
  cannot be verified.

Reactant atoms absent from the main product are retained in
`reactant_element_excess` rather than treated as an error. This allows normal
omission of byproducts. Product-heavy-atom excess is different: it prevents
reaction-signature generation, and conversion makes the record ineligible for
indexing. The system never creates a missing reactant or duplicates a partner
to force reconstruction.

An incomplete record can still retain a typed, observation-only
`partial_product_transformation` when a single conserved scaffold uniquely
supports replacement of one connected branch. Product-only atoms are grouped
into a rooted `ProductOriginGap`, so both `R-C(=O)-OH -> R-C(=O)-Cl` and
multi-atom cases such as `Ar-Br -> Ar-C#N` retain the installed fragment graph
without inventing its reagent. The gap records product-side atom references,
internal bonds, its scaffold attachment, a deterministic key, and unresolved,
agent-supported, or ambiguous source status. Middle-side source matching is
structural support only and does not create supplied atom correspondence. This
is accompanied by a product-heavy-atom provenance ledger covering both the
conserved scaffold and every product-only fragment atom. This
mechanism-neutral observation does not create a named family, a verified
`ReactionSignature`, or normal index eligibility.

The fallback descriptor projects this observation into two deterministic,
mechanism-neutral contracts:

```text
PTS1 partial-transformation key
  = center + removed branch + installed fragment graph + attachment bond

FSR1 fragment-source requirement
  = installed fragment graph + composition + attachment context
```

These keys are derived from molecular structure, not a reaction label. During
dataset conversion, `condition_recommender` evaluates every `FSR1` requirement
against a versioned registry of curated condition-source capabilities. Only a
precedent whose reported reagent/catalyst recipe supports every requirement may
enter the exact partial-transformation index. A capability match is evidence
that the recipe can supply the fragment; it is not asserted atom
correspondence. Merely containing the same element does not qualify a condition
component.

Completeness and mapping warnings include:

```text
PARTIAL_ATOM_MAPPING
PRODUCT_MAPS_MISSING_FROM_REACTANTS
REACTANT_MAPS_MISSING_FROM_PRODUCTS
UNACCOUNTED_PRODUCT_HEAVY_ATOMS
MISSING_REACTANT_SUSPECTED
INSUFFICIENT_REACTANT_MULTIPLICITY
REACTION_COMPLETENESS_UNRESOLVED
PARTIAL_PRODUCT_CORRESPONDENCE
PRODUCT_FRAGMENT_SOURCE_UNRESOLVED
PRODUCT_FRAGMENT_SOURCE_AGENT_SUPPORTED
AMBIGUOUS_PRODUCT_FRAGMENT_SOURCES
PRODUCT_ATOM_SOURCE_UNRESOLVED:<element>
```

When every enumerated reactant grammar contradicts the reported product, those
candidates remain available as rejected evidence but are vetoed as display
labels with `PRODUCT_CONTRADICTED_GRAMMAR_CANDIDATES`. An unreacted handle such
as Ar-I therefore cannot become the displayed transformation merely because a
compatible reactant-only grammar was enumerable.

Reaction SMILES does not encode experimental quantities. A reported value such
as `0.5 equiv` versus `1.0 equiv` must be preserved as separate condition or
procedure data; it is not inferred from molecular component counts.

## 6. Generate the transformation class and named family

These are separate fields:

- `transformation_class` describes the structural transformation, such as
  `sp2_c_n_substitution`.
- `named_family` is an optional chemistry interpretation, such as
  `suzuki_miyaura`.

The current implementation assigns `named_family` only when the selected
grammar has exactly one compatible named family.

For example:

```text
alkyl halide + alcohol -> ether
transformation_class = sp3_c_o_substitution
named_family = None
```

This is intentional because the graph cannot reliably distinguish SN1, SN2,
protection chemistry, and related mechanistic interpretations.

For a valid edit-backed but grammar-unknown reaction, the system can still
produce a `ReactionSignature` while leaving `named_family` unset. When no
specific edit-pattern class is supported, the signature receives
`generic_graph_transformation` or
`generic_multi_event_graph_transformation`, so structural retrieval does not
depend on a named mechanism.

Source reaction names and `source_declared_family` are not inputs to
`featurize_reaction()`. They are retained later as provenance and do not
determine the structural result.

## 7. Build the other reaction features

After edit normalization, the system derives:

- `spectator_groups`: unchanged functional groups outside the selected event;
- `family_environment`: role-specific steric, electronic, nearby-group,
  coordination, and competing-site features;
- `product_connection`: a compatibility view of a newly formed bond;
- `reaction_topology`: intermolecular or intramolecular scope, tether distance,
  ring formation, and ring-size changes;
- `ReactionEvent` objects: connected groups of edits;
- template-free reaction-core atom transitions, events, remote subgraphs, and
  typed attachment ports when mapped evidence is available;
- retained, created, destroyed, or descriptor-changed stereochemistry;
- event relationships such as `shared_atom`, `shared_component`, or
  `independent_sites`;
- warnings, confidence, and evidence quality.

These observations remain available even if the named-family interpretation is
absent.

## 8. Generate the reaction label

There are two label-generation stages.

First, every grammar candidate receives a prospective grammar label:

```text
Ar1-Br + Ar2-B(OH)2 -> Ar1-Ar2
```

This label is rendered from the assigned reactive sites and a declarative
product-rendering rule. It is not accepted as the final label merely because
the reactants match.

Then
[`build_reaction_display_label()`](../reactive_taxonomy/reaction_display_labels.py)
generates the final evidence-aware label.

| Situation | Final label type |
| --- | --- |
| Mapping conflicts with reconstruction | Conflict-labelled structural edit summary |
| Multiple reaction events | Aggregated event labels |
| Exact reconstructed grammar | Full reactant-to-product grammar label |
| Recognized generic edit pattern | `C-N substitution`, `C=C hydrogenation`, etc. |
| Valid edits without a known pattern | Literal before-to-after edit summary such as `C=O + C-H -> C-C + C-O + O-H` |
| Only one plausible reactant assignment | Reactant-only label ending in an arrow |
| Several plausible assignments | `OR`-joined ambiguous reactant labels |
| No usable evidence | No label |

The structured `display_label` also preserves:

- concise and detailed labels;
- individual edit clauses;
- grammar label;
- generic transformation label;
- structural label;
- contextual reactant and product labels;
- pattern and grammar IDs;
- status, evidence, confidence, and warnings.

The final top-level `reaction_label` is the concise projection of this
structured display.

Typical `reaction_label_status` values include:

```text
exact_product
mapped_edit_summary
observed_edit_summary
mapped_generic_pattern
generic_pattern
multi_event_edit_summary
conflicting_edit_summary
reactant_only
ambiguous_reactants
unavailable
```

Display labels do not participate in reaction identity.

For globally inferred correspondence, the label deliberately describes the
observed graph transition instead of claiming a mechanism. The transformation
class uses a recognized generic edit-pattern class when one is uniquely
supported; otherwise it is `generic_graph_transformation`. The named-family
field remains unset.

## 9. Generate the reaction signature

When normalized edits and topology are available,
[`build_reaction_signature()`](../reactive_taxonomy/reaction_signatures.py)
creates deterministic signature keys:

- L0 `exact_signature_key`: edits, explicit stereochemistry, detailed
  environments, partners, events, and topology;
- L1 `handle_signature_key`: reactive handles, explicit stereochemistry, and
  less-specific edit context;
- L2 `transformation_signature_key`: bond changes, hydrogen changes,
  transformation class, events, and topology;
- L3 `bond_edit_signature_key`: topology-agnostic bond-edit fallback;
- L4 `environment_signature_key`: partner environments and spectators.

The final `signature_id` hashes:

```text
L0-L4 keys
+ signature schema version
+ chemistry definition versions
```

It deliberately excludes:

- reaction display labels;
- source reaction names;
- source row order;
- irrelevant reactant ordering;
- serialization formatting.

## 10. Return and serialize the result

The full output is
[`ReactionAnalysis`](../reactive_taxonomy/reaction_models.py), containing:

```text
parsed components
candidate interpretations
selected candidate and events
evidence candidates from every attempted structural provider
distinct ambiguous edit hypotheses
transformation class
compatible named families
named family
reaction label and structured display label
normalized mapped bond changes
spectators
partner environments
product connection
reaction topology
reaction completeness assessment
reaction signature
reaction core projection
evidence quality
warnings and errors
```

During dataset conversion,
[`convert_record()`](../condition_recommender/conversion/generic.py):

- runs reaction featurization;
- resolves conditions through `condition_registry`;
- evaluates typed product-fragment source requirements against curated,
  versioned condition capabilities;
- attaches source provenance;
- evaluates admission;
- rejects unsupported product-atom-incomplete records from indexing;
- admits a uniquely observed, source-supported partial transformation only to
  the separate review-qualified partial index, never as a verified signature;
- sends unresolved provenance, ambiguous correspondence, and inconsistent
  mapping to review or rejection as appropriate;
- admits a fully resolved ingredient set from an unassigned multistage record
  only at review confidence, with a ranking penalty and explicit caution;
- serializes reaction completeness in canonical JSONL and exposes its status,
  coverage, and warnings in generic review CSV;
- serializes the `PTS1` key, `FSR1` requirements, matched condition components,
  capability IDs, and support status;
- serializes the nested reaction signature;
- serializes the nested reaction-core projection and flattened review keys;
- serializes provider evidence and edit hypotheses for review even when no
  signature can be issued;
- optionally reconciles RXNMapper proposals, serializes their full provenance,
  and keeps every mapper-derived precedent review-only;
- exposes L0-L4 keys in review exports;
- places the record into verified, review, or rejected tiers.

The central design principle is:

```text
molecular graph
  -> observed or reconstructed edits
  -> optional interpretation
  -> human-readable label
  -> recommendation features
```

The label explains the structural analysis; it does not create or control that
analysis.

For an unsigned query carrying an exact `PTS1` key, recommendation checks the
partial-transformation index before ordinary structural fallback:

```text
exact PTS1 identity
  -> all FSR1 requirements supported by each precedent recipe
  -> compatibility and independent-support gates
  -> local-environment ranking within that exact transformation
  -> expert-review result or abstention
```

This makes the behavior general across curated fragment replacements while
preventing a superficially similar but different edit—such as amidation—from
competing with acid-to-acyl-fluoride precedents.
