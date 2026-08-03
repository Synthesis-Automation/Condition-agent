# Reaction Featurization Workflow

The system is graph-first, type-agnostic, and label-last. Molecular structures
and normalized graph changes are the source of truth. Molecular reactivity
hypotheses, synthesis patterns, family names, and the reaction label are optional
annotations and cannot create or override structural evidence.

## Canonical workflow

```text
Reaction SMILES
  ↓
Parse molecular graph facts and normalize component inventory
  ↓
Infer atom correspondence
  ↓
Normalize observed edits or retain explicit edit hypotheses
  ↓
Build ReactionObservation
  ├─ topology and completeness
  ├─ minimum ReactionCoreProjection
  │    ├─ edit-provenanced reaction events
  │    ├─ shortest-path relationships between multiple sites
  │    └─ omitted fragments with port-specific R-group profiles
  └─ generic ReactionSignature when evidence is sufficient
  ↓
Add optional annotations
  ├─ molecular motifs and reactive-site hypotheses
  └─ generic transformation and synthesis-pattern matches that cite
     core event IDs, edit indices, and attachment-profile IDs
  ↓
Build one terminal ReactionRenderContext
  ├─ render the reaction label
  └─ render the minimized reaction graphic on demand
  ↓
Serialize for conversion and recommendation
```

This is the execution order in
[`featurize_reaction()`](../reactive_taxonomy/reaction_api.py). Optional molecular
annotations are deliberately attached only after `ReactionObservation` has
been built.

## Layer ownership

| Layer | Owns | Excludes |
| --- | --- | --- |
| Molecular structure | Parsed atoms, bonds, components, maps, charge, aromaticity, stereo | Motifs, reactive sites, reaction families |
| Correspondence and edits | Atom-origin alternatives, formed/broken/order/H changes, confidence, provenance | Named reactions and reaction-site routing |
| `ReactionObservation` | Structure-only components, edits or hypotheses, topology, completeness, minimum core, event relationships, omitted-fragment profiles | Molecular annotations, patterns, family labels, display text |
| `ReactionSignature` | Versioned generic chemistry identity | Display labels, source names, motifs, reactive-site hypotheses, family identity |
| Molecular interpretation | Motif matches, reactive-site hypotheses, local reactivity profiles, connectivity hypotheses | Reaction evidence or identity |
| Reaction interpretation | Patterns supported by existing edits, optional family evidence | Atom correspondence, edits, predicted products |
| Rendering | One chemist-facing reaction label | Chemistry identity or recommendation routing |

## 1. Parse molecular graph facts

[`parse_reaction_smiles()`](../reactive_taxonomy/reaction_parser.py) accepts
`reactants>>products` and `reactants>agents>products`. It calls
[`observe_molecular_structure()`](../reactive_taxonomy/api.py) for each
component and records structure only.

`MolecularStructureObservation` contains parsed atoms, bonds, disconnected
components, canonical SMILES, supplied maps, warnings, and errors.
`MolecularInterpretation` is a separate contract. The composed
`MoleculeAnalysis` is convenient for molecular tools, but only its `structure`
projection may enter a reaction observation.

[`infer_reactant_multiplicity()`](../reactive_taxonomy/reaction_stoichiometry.py)
handles omitted stoichiometric repetition conservatively. It adds a whole
supplied reactant copy only when one bounded structural inventory exactly
accounts for the product-side heavy-element deficit. The inferred copy and its
source component remain explicit, and correspondence must still validate the
expanded inventory. Ambiguous deficits remain unresolved.

## 2. Infer correspondence and normalize edits

[`resolve_structural_evidence()`](../reactive_taxonomy/reaction_edits.py) uses:

1. validated supplied atom mapping;
2. bounded whole-graph, scaffold, and fragmented-scaffold correspondence; and
3. an optional externally mapped proposal after structural validation.

For multi-reactant records, a bounded fallback may expose the largest exact
fragment created by one non-ring attachment cut. This supports tandem records
where a substrate loses a protecting fragment at one locus while reacting at
another. It remains an atom-origin hypothesis and does not assign a reaction
type.

There are no reaction-specific reconstruction rules in this evidence path.
Substitution, elimination, reductive amination, and similar concepts do not
select mappings or generate edits.

When all best correspondences imply the same chemistry, the system emits typed
normalized edits. When atom origins remain chemically distinct, it retains
deterministic `ReactionEditHypothesis` alternatives. Ambiguity is not converted
into a false observation.

## 3. Build the observation, minimum core, and signature

[`build_reaction_observation()`](../reactive_taxonomy/reaction_observation.py)
projects parsed components to `ReactionStructureComponent` and builds:

- provider evidence and correspondence alternatives;
- normalized edits and stereochemical changes;
- generic topology and product completeness;
- the minimum [`ReactionCoreProjection`](../reactive_taxonomy/reaction_core/builder.py);
- a generic [`ReactionSignature`](../reactive_taxonomy/reaction_signatures.py)
  when evidence is sufficient.

For inferred correspondence, deterministic internal atom IDs allow the minimum
core to be generated without mutating or pretending that the input was mapped.
The core keeps active atom transitions, connected events, remote graph shape,
attachment ports, and shortest observed paths between distinct edit events in
the same molecule. Every core event records its normalized edit indices and
reactant/product component provenance. It does not use motifs, reactive-site
hypotheses, source labels, or family names.

Every attachment port carries a versioned `ReactionCoreSubstituentProfile`
derived from the omitted graph rather than its display label. The profile
separates fragment identity from port chemistry and exposes four abstraction
levels:

| Level | Port information |
| --- | --- |
| L0 | alkyl, ring alkyl, aryl, heteroaryl, alkenyl, alkynyl, acyl, or heteroatom |
| L1 | methyl/primary/secondary/tertiary, cyclic, benzylic, allylic, or propargylic |
| L2 | alpha/beta branching, ring sizes, and nearby heteroatoms |
| L3 | deterministic radius-two local-environment key |

`R1`, `R2`, and similar drawing symbols identify continuity. They are not
chemical classifications. The exact fragment SMILES and all atom/port
provenance remain serialized beside the profile.

The signature hashes normalized chemistry and identity-bearing definition
versions. It is invariant to irrelevant component order and serialization. A
missing family never prevents a structurally supported generic signature.

## 4. Add optional annotations

After the observation exists,
[`interpret_parsed_molecules()`](../reactive_taxonomy/reaction_parser.py) may add:

- [`MolecularMotifMatch`](../reactive_taxonomy/models.py) records;
- `ReactiveSiteHypothesis` records;
- local reactivity profiles; and
- optional connectivity hypotheses for downstream reasoning.

These are useful priors for compatibility, explanations, and recommendation,
but disabling their definitions must leave the observation, minimum core, and
signature unchanged.

[`match_reaction_patterns()`](../reactive_taxonomy/reaction_patterns.py) then
matches optional patterns against the completed observation. Each match records
the normalized edit indices, core event IDs, attached substituent-profile IDs,
and fraction of core events that support it. Generic patterns
include net substitution, elimination, addition, coupling, bond cleavage,
bond-order direction, and ring opening or closure. Reusable graph queries live
in `reaction_pattern_predicates.py`; validated v2 metadata lives in
`transformation_patterns.v2.json` and `synthesis_patterns.v2.json`.

Synthesis patterns combine those facts into optional structural annotations.
For example, SNAr, Buchwald–Hartwig C–N, and Ullmann C–N products share one
`sp2_c_n_substitution_like` pattern because their graph-observable change is the
same. The pattern retains all compatible named families and requires condition
evidence to resolve them. Organoboron C–N/O/S coupling is separate and may
support a Chan–Lam candidate, while organoboron C–C coupling supports a
Suzuki-like interpretation.

Pattern definitions contain no operators, structural slots, predicted edits,
or reconstruction instructions. They may rank display interpretations or add
optional family evidence; they cannot modify structural facts.

When two synthesis patterns are supported by non-overlapping edit subsets, the
interpretation retains both as `co_occurring_pattern_ids`. They remain typed
annotations and do not replace graph-derived reaction rendering. Distinct
patterns may occupy non-overlapping edit subsets of one connected core event;
this supports tandem operations on the same atom without inventing a false
event split.

## 5. Build terminal presentation context, render, and serialize

[`build_reaction_render_context()`](../../reactive_taxonomy/reaction_render_context.py)
combines the completed observation, minimum core, signature, annotated
components, topology, partial-product evidence, interpretation, and notation
style. It is the sole typed input shared by the sibling terminal renderers.

```text
Completed reaction analysis
  ↓
ReactionRenderContext
  ├─ render_reaction()                 → RenderedReactionLabel.text
  └─ build_reaction_core_graphic()     → minimized SVG or PNG
```

The minimum `ReactionCoreProjection` remains structural evidence built in
stage 3. Only its graphical presentation occurs here. The graphic and text do
not parse or depend on each other.

[`render_reaction()`](../../reactive_taxonomy/reaction_rendering.py) is the sole
terminal reaction renderer. `RenderedReactionLabel.text` is the only complete
reaction label; `status`, `basis`, `evidence`, `confidence`, and `warnings` are
metadata. The renderer orders and selects normalized edits through the core
events and serializes the contributing core event IDs, R-profile IDs, pattern
IDs, and any unclassified edit indices. It does not read source reaction names
or use optional named-family annotations to override graph evidence.

[`build_reaction_core_graphic()`](../../visualization/reaction_core_graphic.py)
uses the same context, including the same canonical notation style. Its own
definition controls only graphical abstraction policy; it contains no private
`R`/`Alk`/`Ar`/`HetAr` vocabulary.

All molecular-site and reaction rendering resolves shorthand through
`chemist_notation.v1.json`. The canonical fragment vocabulary is:

| Symbol | Definition |
|---|---|
| `R` | Organic substituent whose more specific class is unresolved |
| `Alk` | Non-aromatic substituent attached through an sp3 carbon |
| `c-Alk` | Saturated cyclic substituent attached through an aliphatic ring carbon |
| `Ar` | Carbocyclic aromatic group attached through an aromatic carbon |
| `HetAr` | Heteroaromatic ring system attached through an aromatic carbon |
| `Bn`, `Allyl`, `Propargyl` | Benzyl, allylic, and propargylic substituents at the named attachment atom |
| `X` | F, Cl, Br, or I |
| `Het` | Heteroatom whose element is unresolved |

There are no notation aliases. In particular, `HeteroAr` is invalid; the sole
heteroaryl symbol is `HetAr`.

Review CSV exports use the same label contract:

```text
reaction_label
reaction_label_status
reaction_label_basis
reaction_label_confidence
reaction_label_warnings
```

They also expose optional interpretation without parsing display text:

```text
primary_reaction_pattern
primary_reaction_pattern_count
reaction_pattern_matches
co_occurring_reaction_patterns
identified_reaction_type
compatible_reaction_types
reaction_pattern_confidence
reaction_pattern_requires_condition_evidence
```

For example, Suzuki chemistry may resolve to `suzuki_miyaura`, while an sp2
C–N substitution retains SNAr, Buchwald–Hartwig, and Ullmann as compatible
types until condition evidence distinguishes them.

There are no separate core, graphic, concise, detailed, or contextual reaction
labels. Canonical nested JSON or Parquet remains the lossless converted artifact;
CSV is a review view.

## Invariants

- Parse graph facts before running molecular annotation definitions.
- Correspondence and normalized edits do not depend on reaction patterns.
- The minimum core is structured graph evidence and contains no display label.
- R-group chemistry is derived per attachment port; `R1`/`R2` remain continuity
  labels only.
- Molecular annotations and reaction patterns are optional overlays.
- Patterns consume core-provenanced edits and cite their event/profile evidence;
  they never generate edits.
- Ambiguity, conflicts, confidence, and provenance remain explicit.
- Rendering explains chemistry but does not define chemistry identity.
- Text and minimized graphics consume the same `ReactionRenderContext` and
  canonical notation registry.
