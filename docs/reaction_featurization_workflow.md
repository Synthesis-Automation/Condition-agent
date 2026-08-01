# Reaction Featurization Workflow

The system is graph-first, type-agnostic, and label-last. Molecular structures
and normalized graph changes are the source of truth. Molecular reactivity
hypotheses, synthesis patterns, family names, and display labels are optional
annotations and cannot create or override structural evidence.

## Canonical workflow

```text
Reaction SMILES
  ↓
Parse molecular graph facts
  ↓
Infer atom correspondence
  ↓
Normalize observed edits or retain explicit edit hypotheses
  ↓
Build ReactionObservation
  ├─ topology and completeness
  ├─ minimum ReactionCoreProjection
  └─ generic ReactionSignature when evidence is sufficient
  ↓
Add optional annotations
  ├─ molecular motifs and reactive-site hypotheses
  └─ generic transformation and synthesis-pattern matches
  ↓
Render one concise/detailed reaction label
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
| `ReactionObservation` | Structure-only components, edits or hypotheses, topology, completeness, minimum core | Molecular annotations, patterns, family labels, display text |
| `ReactionSignature` | Versioned generic chemistry identity | Display labels, source names, motifs, reactive-site hypotheses, family identity |
| Molecular interpretation | Motif matches, reactive-site hypotheses, local reactivity profiles, connectivity hypotheses | Reaction evidence or identity |
| Reaction interpretation | Patterns supported by existing edits, optional family evidence | Atom correspondence, edits, predicted products |
| Rendering | One chemist-facing concise and detailed label | Chemistry identity or recommendation routing |

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

## 2. Infer correspondence and normalize edits

[`resolve_structural_evidence()`](../reactive_taxonomy/reaction_edits.py) uses:

1. validated supplied atom mapping;
2. bounded whole-graph, scaffold, and fragmented-scaffold correspondence; and
3. an optional externally mapped proposal after structural validation.

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
and attachment ports. It does not use motifs, reactive-site hypotheses, source
labels, or family names.

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
matches optional patterns against the completed observation. Generic patterns
include net substitution, elimination, addition, coupling, bond-order change,
and ring closure. More specific synthesis patterns include reductive-amination-
like, amide-formation-like, boron-transfer-coupling-like, Heck-like,
cycloaddition-like, and decarboxylative-coupling-like observations.

Pattern definitions contain no operators, structural slots, predicted edits,
or reconstruction instructions. They may rank display interpretations or add
optional family evidence; they cannot modify structural facts.

## 5. Render and serialize

[`render_reaction()`](../reactive_taxonomy/reaction_rendering.py) returns one
`RenderedReactionLabel` with `concise` and `detailed` text plus evidence,
confidence, provenance, and warnings. Rendering may polish a minimum-core label
or add a supported pattern overlay, but display text never participates in core,
signature, conversion admission, or retrieval identity.

Review CSV exports use the same label contract:

```text
reaction_core_label
reaction_display_label
reaction_display_label_detailed
```

The first is the low-level minimum-core label. The latter two are concise and
detailed projections from the single renderer. Canonical nested JSON or Parquet
remains the lossless converted artifact; CSV is a review view.

## Invariants

- Parse graph facts before running molecular annotation definitions.
- Correspondence and normalized edits do not depend on reaction patterns.
- The minimum core is the base generic reaction representation.
- Molecular annotations and reaction patterns are optional overlays.
- Patterns consume edits; they never generate them.
- Ambiguity, conflicts, confidence, and provenance remain explicit.
- Rendering explains chemistry but does not define chemistry identity.
