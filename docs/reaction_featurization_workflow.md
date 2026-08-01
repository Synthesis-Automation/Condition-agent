# Reaction Featurization Workflow

The system is graph-first and label-last. Reaction names and source labels are
optional annotations, never primary structural evidence.

```text
Reaction SMILES
  ↓
Parse molecular graphs
and featurize every component
  ↓
Enumerate grammar-free structural reconstructions
  ↓
Resolve normalized edits
  ↓
Build the structural observation
(alternatives, topology, minimum core, completeness, spectators)
  ↓
Build the generic signature and retrieval keys
  ├─ sufficient evidence → versioned ReactionSignature
  └─ insufficient evidence → structural fallback or abstention
  ↓
Optional grammar and family interpretation
  ↓
Render concise and detailed labels
```

Structural evidence includes detected sites, supplied mapping, product
reconstruction, and graph correspondence. The reconciled evidence is stored in
one grammar-free `ReactionObservation`, which owns the alternatives, topology,
minimum `ReactionCoreProjection`, completeness, and spectators. The generic
signature is built from that observation. Ambiguous or conflicting alternatives
remain explicit review evidence rather than being erased. Optional
interpretation may then add semantic roles, grammar and family meaning, while
rendering produces one structured `reaction_label` with concise and detailed
text.

This is also the strict execution order in `featurize_reaction()`. Structural
reconstruction rules contain no grammar, family, transformation-class, or
display metadata. Grammar annotations are loaded only after the observation,
minimum core, topology, and signature have been built.

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
side into components and calls `featurize_molecule()` while constructing each
component. It preserves:

- component index;
- original and canonical SMILES;
- reactant, agent, or product side;
- supplied atom-map numbers;
- parsing errors and warnings.

Agents are analyzed, but they are not used as reaction partners during
structure-evidence collection or reconstruction-candidate enumeration.

Parsing validates syntax and component SMILES; it does not require a globally
balanced chemical equation. Reaction records commonly report only the main
product and omit salts, leaving-group products, and other byproducts. Product
atom provenance is therefore evaluated later by a dedicated
`ReactionCompletenessAssessment` instead of by naive formula equality at this
stage.

## 2. Featurize every molecular component

Each component passes through
[`featurize_molecule()`](../reactive_taxonomy/api.py).

For each molecular graph, the system detects and emits:

- functional groups;
- reactive sites and handles;
- atom roles within each site;
- site availability, such as available or blocked;
- aromatic, alkyl, carbonyl, heteroatom, and other local contexts;
- local steric and electronic environments;
- chemist-readable site labels;
- canonical reactive links, bond capacities, and connection endpoints; and
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

## 3. Enumerate structural reconstruction candidates

[`enumerate_reconstruction_candidates()`](../reactive_taxonomy/reaction_reconstruction.py)
compares canonical connectivity sites with the declarative structural rules in
[`reaction_reconstruction_rules.v1.json`](../reactive_taxonomy/definitions/reaction_reconstruction_rules.v1.json).

Each rule specifies:

- anonymous structural slots;
- site-type constraints;
- allowed handle and context combinations;
- same-component or different-component relationships; and
- a registered graph operator with explicit slot bindings.

For an aryl bromide and aryl boronic acid, a candidate may contain:

```text
rule_id: boron_transfer_coupling
operator_id: suzuki_release_and_connect
slots:
  slot_1 -> aryl bromide site
  slot_2 -> aryl boronic acid site
```

At this stage the candidate is only structurally plausible. It contains no
chemist-facing role names, named family, transformation class, or display
label.

## 4. Apply the connectivity rewrite and verify the product

Each candidate's registered operator is passed to
[`apply_reaction_operator()`](../reactive_taxonomy/connectivity_rewrite.py).

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

Multiple rewrites may also be composed for multi-event reactions. The
composer reads each registered rewrite's declared discardable attachments and
formed bond, so it supports both one-release substitutions and two-release
transfer couplings. It consumes distinct reactive sites and partner component
instances for every event, then accepts the composition only when the combined
operator reconstructs the observed product exactly. It never duplicates a
missing reactant to balance an equation.

## 5. Resolve edits and build the structural observation

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
the providers run through `resolve_reaction_evidence()`.

After reconciliation,
[`build_reaction_observation()`](../reactive_taxonomy/reaction_observation.py)
assembles the grammar-free `ReactionObservation`. In one place it builds or
retains:

- normalized edits, stereo changes, evidence providers, and edit hypotheses;
- generic reaction topology;
- the minimum reaction-core projection;
- product-atom completeness;
- observed spectator groups; and
- structural reconstruction candidates and their selected single- or
  multi-event reconstruction.

The stored reconstruction candidates use anonymous structural slots and do not
contain semantic roles, transformation classes, named families, or display
labels. An interpretation failure therefore cannot erase or change the
observation.

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

- `external_mapping_internal_consensus`: exactly one internal edit hypothesis,
  or an already resolved reaction signature requested in shadow mode, matches
  the external normalized edit profile;
- `external_mapping_only`: the mapper supplies a valid edit profile where no
  internal hypotheses exist;
- `external_mapping_hypothesis_conflict` or
  `external_mapping_ambiguous_hypothesis_match`: retain the original
  hypotheses and route to review;
- `external_mapping_signature_conflict`: retain an already resolved analysis
  when a shadow mapping proposal contradicts its reaction signature;
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
[`build_reaction_core_projection()`](../reactive_taxonomy/reaction_core/builder.py)
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

Its algorithm-v5 concise label uses a local atom-state transition for a
single-center event. For a multi-center event it instead renders one
normalized edit equation, preventing a shared formed bond from being shown
twice (for example, `Ar–B + Ar–O → Ar–Ar`).
Repeated mapped events are retained in one projection and counted rather than
collapsed; a double Suzuki core is rendered as
`2 × Ar–B + 2 × Ar–Br → 2 × Ar–Ar` with `event_count=2`.
Carbonyl-bearing active neighbors retain their local heteroatom environment in
the single-center label. A decarboxylative connection can therefore appear as
`C(H)2(R)(C(=O)(O-H)) -> C(H)2(R)(ArC)`, instead of reducing the departing
carboxyl group to an indistinguishable `R` substituent.

Conservatively recognized motifs may also expose two graph-derived abstraction
layers. Decarboxylative C-C coupling uses the broad label
`R-C(=O)OH + Ar-H -> R-Ar` and a stable `RCM1` motif key across alkyl and aryl
carboxylic acids. Separate limiter tokens retain distinctions such as
`transfer_center:primary_alkyl`, `transfer_center:secondary_alkyl`,
`transfer_center:aryl`, and `partner_center:heteroaryl`. The GUI therefore
shows a readable limiter such as
`R = R'-CH2- (primary alkyl); Ar = HetAr` alongside the atom-level core. The
motif key is serialized as foundation for a future retrieval tier; current
admission and retrieval behavior is unchanged.

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

The reaction CLI and desktop featurizer display this observation in a
`Reaction minimization` section. The concise view includes the minimized label,
evidence status and confidence, `RSH2`, diagnostic `RCS2`, event/center/active
atom counts, remote-subgraph classes, continuity, exact fragment SMILES,
attachment-port counts, and core warnings. When no core is available, it says
that mapped edit correspondence is required instead of silently omitting the
section.

Reaction batch CSV exports include `reaction_core_available`, the four core
keys, minimized label, evidence status, event and center counts, remote classes
and subgraph summaries, and warnings. Nested JSON remains the complete
lossless representation.

### Graphical minimized reaction

The desktop featurizer can render a `ReactionCoreProjection` as a compact
reaction scheme. The renderer is implemented in
[`visualization/reaction_core_graphic.py`](../visualization/reaction_core_graphic.py);
it is a presentation of the existing molecular observation, not another
reaction classifier.

The renderer:

1. draws every active atom and active-atom bond on each side;
2. replaces only `retained` remote subgraphs with short placeholders;
3. pairs the same retained mapped subgraph across both sides;
4. assigns deterministic labels such as `Ar`, `HetAr`, `R`, or indexed
   variants such as `R1` and `R2`; and
5. retains a visible legend from every placeholder to its exact fragment
   SMILES.

The remote class is read from the graph-derived core projection. Rendering
does not infer `aryl` or `alkyl` from reaction names and does not hide
departing, appearing, changed, or unresolved fragments. Placeholder vocabulary
and abstraction policy are versioned in
[`reaction_core_graphic.v1.json`](../visualization/definitions/reaction_core_graphic.v1.json).
Configured one-atom retained heteroatom fragments remain explicit, preserving
chemically important atoms such as a carbonyl oxygen.
Departing, appearing, changed, and unresolved subgraphs are always drawn
explicitly. Thus, atoms omitted from the reported main product remain visible
on the reactant side of a minimized scheme. The renderer does not invent a
specific side product such as `CO2` or `H2O` when that product was not supplied.

A retained ring or scaffold remainder may have multiple attachment ports. If
all its ports surround one active atom, that atom is rendered as the scaffold
placeholder while its bonds to other active atoms are preserved. This turns a
ring-site transformation into a readable scheme such as
`Ar-CHO + HetAr-H -> Ar-C(=O)-HetAr`. When mapped identity shows that multiple
retained ring paths belong to one scaffold spanning several active atoms, the
renderer reunites the paths and active atoms into one shared placeholder. For
example, a double coupling is drawn conceptually as
`Br-Ar3-Br + Ar1-B + Ar2-B -> Ar1-Ar3-Ar2`, rather than duplicating the central
ring or inventing separate placeholders for its two paths. Other multi-port
fragments retain one boundary node per port.

The featurizer has separate `Full structure` and `Minimized reaction` tabs.
For a resolved but unmapped input, `Map resolved reactions for minimized
graphic` may be enabled to run the optional mapper solely to obtain the
mapped-core evidence needed for this view. This is off by default because it
is slower and because an externally mapped core remains a review proposal.
The mapped edits are reconciled with the existing reaction signature; on
consensus, the original interpretation, detailed label, and signature remain
unchanged and only the external core and provenance are attached. A conflict
retains the original analysis and does not silently replace its chemistry.
Before core construction, validated external map numbers are projected onto
the original component and atom ordering. This keeps core coordinates aligned
with the unmapped structures consumed by the graphic renderer even when the
mapper canonicalizes or reorders its SMILES output.
The tab displays evidence status, confidence, the exact placeholder legend,
and the expert-review warning.

For the complex alkyne/azide example in
[`click_reaction_core_graphic_poc.py`](../benchmarks/click_reaction_core_graphic_poc.py),
the five active atoms render as:

```text
R1-C#CH + R2-N3 -> 1,2,3-triazole(R1,R2)
```

The generated PNG, SVG, and evidence report are
`results/click_minimized_reaction_poc.*`. The current RXNMapper confidence for
this example is about `0.291`, so the chemically useful graphic is explicitly
shown as external review evidence rather than verified reaction identity.

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

Reaction SMILES does not encode experimental quantities. A reported value such
as `0.5 equiv` versus `1.0 equiv` must be preserved as separate condition or
procedure data; it is not inferred from molecular component counts.

## 6. Build the generic reaction signature

When the observation has usable normalized edits, topology, and adequate
product completeness,
[`build_observation_signature()`](../reactive_taxonomy/reaction_signatures.py)
creates a grammar-independent `ReactionSignature`. It delegates the typed key
construction to `build_reaction_signature()`, but its only chemistry input is
the already-built `ReactionObservation`.

The deterministic retrieval levels are:

- L0 `exact_signature_key`: edits, explicit stereochemistry, detailed local
  environments, generic partners, events, and topology;
- L1 `handle_signature_key`: reactive handles, explicit stereochemistry, and
  less-specific edit context;
- L2 `transformation_signature_key`: bond changes, hydrogen changes, generic
  structural transformation class, events, and topology;
- L3 `bond_edit_signature_key`: topology-agnostic bond-edit fallback; and
- L4 `environment_signature_key`: generic partner environments and spectators.

The final `signature_id` hashes:

```text
L0-L4 keys
+ signature schema version
+ identity-bearing chemistry definition versions
```

It deliberately excludes reaction display labels, grammar labels, named
families, source reaction names, source row order, irrelevant reactant order,
and serialization formatting. Removing all optional grammar annotations must
therefore leave the same signature ID and minimum core.

## 7. Evaluate structural fallback evidence

If no exact reconstruction or normalized edit set is available, the system may
derive the conservative `partial_product_transformation` described under
product-atom completeness. This happens after the generic signature attempt and
before optional interpretation. It does not create a verified signature or a
named family.

Later, the serialized analysis also receives a mechanism-neutral
`ReactionFallbackDescriptor`. That descriptor projects verified signatures,
partial transformations, retained edit hypotheses, or structure inventory into
explicit fallback evidence for review and retrieval. It never upgrades weak
evidence into an observed transformation.

## 8. Add optional grammar and family interpretation

These are separate fields:

- `transformation_class` describes the structural transformation, such as
  `sp2_c_n_substitution`.
- `named_family` is an optional chemistry interpretation, such as
  `suzuki_miyaura`.

[`reaction_grammar_annotations.v1.json`](../reactive_taxonomy/definitions/reaction_grammar_annotations.v1.json)
maps a structural reconstruction rule to semantic roles, a transformation
class, rendering metadata, and compatible named families. It cannot change the
observation, minimum core, topology, or signature already produced. The system
assigns `named_family` only when the selected annotation has exactly one
compatible named family.

[`build_reaction_interpretation_candidates()`](../reactive_taxonomy/reaction_interpretation.py)
maps lower-level reconstruction slots onto chemist-facing roles. The resulting
`ReactionInterpretation` owns:

- grammar candidates and the selected interpreted candidate or events;
- semantic reaction partners and role confidence;
- compatible named families and the optional selected family;
- the role-specific family environment;
- the role-labelled product-connection view; and
- interpretation-only warnings and conflicts.

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

## 9. Render the reaction label

Interpretation and rendering remain separate stages.

First, every grammar candidate receives a prospective grammar label:

```text
Ar1-Br + Ar2-B(OH)2 -> Ar1-Ar2
```

This label is rendered from the assigned reactive sites and a declarative
product-rendering rule. It is not accepted as the final label merely because
the reactants match.

Then the sole public renderer,
[`render_reaction()`](../reactive_taxonomy/reaction_rendering.py), generates one
final evidence-aware `RenderedReactionLabel`.

| Situation | Final label type |
| --- | --- |
| Mapping conflicts with reconstruction | Conflict-labelled structural edit summary |
| Multiple reaction events | Aggregated event labels |
| Exact reconstructed grammar | Full reactant-to-product grammar label |
| Recognized generic edit pattern | `C-N substitution`, `C=C hydrogenation`, etc. |
| Valid edits without a known pattern | Literal before-to-after edit summary such as `C=O + C-H -> C-C + C-O + O-H` |
| Only one plausible reactant assignment | Reactant-only label ending in an arrow |
| Several plausible assignments | `OR`-joined ambiguous reactant labels |
| Product contradicts candidates | Reactant-side handles only, with contradiction status |
| No usable evidence | Explicit `Unavailable` label |

The structured `reaction_label` preserves:

- concise and detailed labels;
- individual edit clauses;
- grammar label;
- generic transformation label;
- structural label;
- contextual reactant and product labels;
- pattern and grammar IDs;
- status, evidence, confidence, and warnings.

Its `status` values include:

```text
exact_reconstruction
family_overlay
observed_edits
generic_pattern
multi_event
ring_formation
core_projection
partial_product_correspondence
reactant_only
ambiguous_reactants
product_contradicted_reactants
unavailable
```

There is no parallel string label or standalone label-status field. JSON keeps
the nested object; review CSV projects only its concise and detailed text as
`reaction_display_label` and `reaction_display_label_detailed`. Display labels
do not participate in reaction identity.

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

In compact form, the architectural flow is:

```text
parse
  -> collect structural evidence
  -> resolve normalized edits
       consistent -> minimum core + topology -> signature/retrieval keys
       ambiguous  -> retained alternatives; review or abstain
  -> optional grammar interpretation
  -> rendering
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
