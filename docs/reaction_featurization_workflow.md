# Reaction Featurization Workflow

The reaction featurization system is graph-first and label-last. A reaction
name or source label is never the primary classifier.

```text
Reaction SMILES
  -> parse molecules
  -> detect functional groups and reactive sites
  -> enumerate grammar-compatible candidates
  -> reconstruct possible products
  -> extract and reconcile graph edits
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

Agents are analyzed, but they are not used as reaction partners during grammar
candidate enumeration.

Parsing validates syntax and component SMILES; it does not require a globally
balanced chemical equation. Reaction records commonly report only the main
product and omit salts, leaving-group products, and other byproducts. Product
atom provenance is therefore evaluated later by a dedicated
`ReactionCompletenessAssessment` instead of by naive formula equality at this
stage.

## 2. Featurize every molecule

Each component passes through
[`featurize_molecule()`](../reactive_taxonomy/api.py).

For each molecular graph, the system detects:

- functional groups;
- reactive sites and handles;
- atom roles within each site;
- site availability, such as available or blocked;
- aromatic, alkyl, carbonyl, heteroatom, and other local contexts;
- local steric and electronic environments;
- chemist-readable site labels.

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
compares the detected sites with the declarative grammars in
[`reaction_grammars.v1.json`](../reactive_taxonomy/definitions/reaction_grammars.v1.json).

Each grammar specifies:

- required partner roles;
- site-type constraints;
- allowed handle and context combinations;
- same-component or different-component relationships;
- transformation class;
- graph operator;
- net edit archetype (`substitution`, `addition`, or `elimination`);
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

## 4. Apply the graph operator and verify the product

Each candidate is passed to
[`apply_operator()`](../reactive_taxonomy/reaction_operators.py).

The operator:

- removes or changes required bonds;
- forms new bonds;
- records hydrogen changes;
- produces a predicted product;
- returns predicted bond edits.

The common operator registry contains `center_replacement`, `pair_addition`,
and `pair_elimination`. Pair addition accepts explicit A–B donors and implicit
A–H donors through one normalized contract and may return multiple
constitutional outcomes. Product reconstruction selects a unique orientation;
reactants alone do not force regio- or stereochemistry.

The predicted product is canonicalized without atom maps and compared with the
observed product. Candidate verification becomes one of:

- `exact_product_reconstruction`;
- `product_mismatch`;
- `construction_failed`.

A candidate is selected only when the structural evidence is sufficiently
unique. Symmetry-equivalent candidates can collapse to one interpretation;
chemically different exact candidates remain ambiguous.

Multiple operators may also be composed for multi-event reactions.

## 5. Extract and reconcile observed reaction edits

[`normalize_reaction_edits()`](../reactive_taxonomy/reaction_edits.py)
establishes the actual transformation evidence.

The effective evidence priority is:

1. Valid supplied atom mapping.
2. Exact single-event operator reconstruction.
3. Exact multi-event operator reconstruction.
4. Conservative unique-scaffold correspondence.
5. Bounded global multi-reactant correspondence.
6. Unresolved or ambiguous evidence.

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

Normalized edit types include:

- `formed`;
- `broken`;
- `order_changed`;
- `hydrogen_change`.

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
supports replacement of one terminal attachment. For example,
`R-C(=O)-OH -> R-C(=O)-Cl` can be reported as an
`acyl_heteroatom_substitution` with `Cl` explicitly marked as having no
reactant source. This mechanism-neutral observation does not create atom
provenance, a named family, a verified `ReactionSignature`, or index
eligibility.

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
evidence quality
warnings and errors
```

During dataset conversion,
[`convert_record()`](../condition_recommender/conversion/generic.py):

- runs reaction featurization;
- resolves conditions through `condition_registry`;
- attaches source provenance;
- evaluates admission;
- rejects product-atom-incomplete records from indexing and sends unresolved
  provenance or inconsistent mapping to review;
- admits a fully resolved ingredient set from an unassigned multistage record
  only at review confidence, with a ranking penalty and explicit caution;
- serializes reaction completeness in canonical JSONL and exposes its status,
  coverage, and warnings in generic review CSV;
- serializes the nested reaction signature;
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
