# Reaction Featurization Workflow

The reaction system is graph-first, type-agnostic, and label-last. Molecular
structure is the source of truth. Reaction names, interpretation annotations, named
families, and display labels cannot create or override structural evidence.

## Canonical workflow

```text
Reaction SMILES
  ↓
Parse and featurize molecular graphs
  ↓
Enumerate interpretation-independent structural reconstructions
  ↓
Resolve mapping, reconstruction, and correspondence into normalized edits
  ↓
Build ReactionObservation
  ├─ evidence and alternatives
  ├─ topology and minimum reaction core
  └─ completeness and spectators
  ↓
Build generic ReactionSignature or retain partial-product evidence
  ↓
Add optional interpretation annotations and named-family evidence
  ↓
Render one concise/detailed reaction label
  ↓
Serialize ReactionAnalysis for conversion and recommendation
```

This is the execution order in
[`featurize_reaction()`](../reactive_taxonomy/reaction_api.py). Interpretation
annotations are loaded only after the observation, minimum core, topology, and
generic signature have been built.

## Contract ownership

| Layer | Owns | Must not own |
| --- | --- | --- |
| Molecular graphs | Functional groups, reactive sites, local environments, connectivity interfaces | Reaction families or reaction labels |
| Structural reconstruction | Anonymous site slots, component relationships, graph operators, predicted products and edits | Semantic partner roles, named families, display text |
| Evidence resolution | Normalized edits, alternatives, conflicts, confidence, provenance | Forced correspondence or family routing |
| `ReactionObservation` | Edits, hypotheses, topology, minimum core, completeness, spectators, reconstruction evidence | Interpretation labels or named-family identity |
| `ReactionSignature` | Versioned generic chemistry and retrieval identity | Display labels, source names, interpretation or family identity |
| `ReactionInterpretation` | Semantic roles, transformation annotation, family candidates, family environment, product-connection view | Changes to the observation, core, or signature |
| Rendering | One chemist-facing concise/detailed label | Chemistry identity or recommendation routing |

## 1. Parse and featurize molecular graphs

[`parse_reaction_smiles()`](../reactive_taxonomy/reaction_parser.py) accepts:

```text
reactants >> products
reactants > agents > products
```

It splits each side into components and calls
[`featurize_molecule()`](../reactive_taxonomy/api.py) for every component. Each
`ReactionComponent` retains its side, component index, input and canonical
SMILES, supplied atom-map status, and molecular analysis.

Molecular featurization detects:

- functional groups and reactive sites;
- atom roles and site availability;
- canonical reactive links, bond capacities, and connection endpoints; and
- local structural, steric, and electronic environments.

These are molecular observations, not reaction-family assignments. Agents are
featurized for provenance and later condition analysis, but they are not used
as reactant partners during structural reconstruction.

## 2. Enumerate structural reconstructions

[`enumerate_reconstruction_candidates()`](../reactive_taxonomy/reaction_reconstruction.py)
matches reactant sites to the anonymous slots in
[`reaction_reconstruction_rules.v1.json`](../reactive_taxonomy/definitions/reaction_reconstruction_rules.v1.json).
Each rule contains only:

- site and availability constraints;
- same- or different-component relationships;
- a registered graph-operator ID; and
- bindings from operator inputs to anonymous structural slots.

The selected operator from
[`connectivity_rewrites.v3.json`](../reactive_taxonomy/definitions/connectivity_rewrites.v3.json)
is executed by `apply_reaction_operator()`. It may break or form bonds, change
bond order, adjust schema-level hydrogen or charge state, and predict a product.

A `ReactionReconstructionCandidate` records the structural rule, operator,
anonymous slot assignments, predicted edits and product, verification status,
and warnings. It deliberately has no interpretation annotation ID, semantic reaction role,
transformation class, named family, or display label.

Exact single-event and balanced multi-event reconstruction are supported.
Multi-event reconstruction consumes distinct site and component instances; it
never duplicates a missing reactant to force product agreement.

## 3. Resolve normalized reaction evidence

[`resolve_reaction_evidence()`](../reactive_taxonomy/reaction_edits.py) reconciles
three structural evidence sources:

1. validated supplied atom mapping;
2. exact single- or multi-event operator reconstruction; and
3. conservative graph correspondence for otherwise unresolved products.

Agreement strengthens confidence. Contradictory mapping and reconstruction are
retained as conflicting evidence; neither is silently discarded. Ambiguous
correspondence is stored as deterministic `ReactionEditHypothesis` alternatives
rather than converted into invented observed edits.

Normalized evidence can include formed, broken, order-changed, hydrogen, charge,
and stereochemical changes. Each atom reference retains component and atom
provenance, element, charge, aromaticity, hybridization, mapping, and local
environment identity.

Optional external atom mapping is a separate review/query enrichment path. It
must pass structural validation and reconciliation and cannot silently replace
the internal analysis. Operational details belong in the primary
[type-agnostic implementation document](new/type_agnostic_reaction_recommendation_implementation.md),
not in this workflow.

## 4. Build the structural observation and minimum core

[`build_reaction_observation()`](../reactive_taxonomy/reaction_observation.py)
assembles one interpretation-independent `ReactionObservation` containing:

- parsed reactants, agents, and products;
- normalized edits and stereochemical changes;
- provider evidence and unresolved edit hypotheses;
- reconstruction candidates and selected reconstruction events;
- generic reaction topology;
- product-atom completeness;
- spectator groups; and
- the minimum `ReactionCoreProjection`, when supported.

The topology describes component scope, participating components, formed-ring
changes, tether distances, and other generic graph facts. Product completeness
is `verified`, `incomplete`, or `unresolved`; missing reported product atoms are
not invented, and omitted byproducts are not treated as missing main products.

[`build_reaction_core_projection()`](../reactive_taxonomy/reaction_core/builder.py)
minimizes the reaction directly from molecular graphs and normalized edits. It
retains active atom transitions, connected edit events, remote subgraphs, and
typed attachment ports. It does not load an interpretation annotation, source label, or
reaction name.

An interpretation failure cannot erase or modify any observation field.

## 5. Build generic identity and partial-product evidence

When normalized edits, topology, and completeness are adequate,
[`build_observation_signature()`](../reactive_taxonomy/reaction_signatures.py)
builds a generic `ReactionSignature` directly from the observation.

Its deterministic retrieval levels are:

- L0: exact edits, detailed environments, events, and topology;
- L1: reactive handles and less-specific edit context;
- L2: generic structural transformation, events, and topology;
- L3: topology-agnostic bond-edit identity; and
- L4: generic partner environments and spectators.

`signature_id` hashes the normalized L0-L4 chemistry, schema version, and
identity-bearing definition versions. It excludes display labels, interpretation
labels, named families, source reaction names, row order, and irrelevant
reactant serialization or ordering.

If verified signature evidence is unavailable, the system retains explicit
structural evidence instead of guessing:

- a unique partial product transformation for conservative branch replacement;
- or the unresolved edit hypotheses already stored in the observation.

These fallbacks can support review or bounded retrieval, but they do not become
verified observations or named reactions.

## 6. Add optional interpretation annotations

Only after generic identity exists does
[`build_reaction_interpretation_candidates()`](../reactive_taxonomy/reaction_interpretation.py)
apply
[`reaction_interpretation_annotations.v1.json`](../reactive_taxonomy/definitions/reaction_interpretation_annotations.v1.json).

An annotation maps anonymous reconstruction slots to chemist-facing roles and
may add:

- a more specific transformation-class interpretation;
- compatible named families and an optional uniquely supported family;
- semantic reaction partners and role confidence;
- a role-specific family environment; and
- a role-labelled product-connection view.

Family identity remains optional. A valid unknown-family reaction can retain
its core and generic signature with `named_family=None`. Interpretation
conflicts are reported and cannot change structural edits, core identity, or
signature identity.

## 7. Render one reaction label

[`render_reaction()`](../reactive_taxonomy/reaction_rendering.py) is the sole
public reaction-label renderer. It combines structural evidence with any valid
optional interpretation and returns one `RenderedReactionLabel` containing:

- `concise` and `detailed` text;
- structured edit clauses;
- structural, contextual, interpretation, and family overlays when supported;
- status, evidence, confidence, provenance IDs, and warnings.

The renderer prefers the strongest non-contradicted evidence. Typical outputs
include exact interpreted labels, multi-event summaries, generic edit-pattern
labels, ring-formation labels, minimum-core labels, partial-product labels, or
an explicit unavailable result.

A reactant-side annotation that contradicts the reported product cannot become
the final product label. Display text never participates in core, signature, or
retrieval identity.

Review CSV exports project the same nested label contract into:

```text
reaction_display_label
reaction_display_label_detailed
```

There is no parallel reaction-label pipeline.

## 8. Serialize for conversion and recommendation

The API returns one
[`ReactionAnalysis`](../reactive_taxonomy/reaction_models.py) containing:

```text
parsed molecular components
structural observation and minimum core
generic signature or explicit fallback evidence
optional interpretation
one rendered reaction label
evidence quality, warnings, and errors
```

Dataset conversion serializes the nested chemistry contracts before condition
normalization, admission, indexing, and retrieval. The review CSV is only a
flattened inspection view; canonical JSON or Parquet remains the lossless
artifact.

Before returning the analysis,
[`build_reaction_fallback_descriptor()`](../reactive_taxonomy/reaction_fallback_descriptors.py)
projects the available signature, partial transformation, interpreted
candidates, retained hypotheses, or structure inventory into one explicit
fallback descriptor. This late projection does not alter any earlier chemistry
contract.

`condition_recommender` may use verified signature tiers, conservative fallback
evidence, chemistry compatibility, and optional family information. It must not
derive chemistry identity from the rendered label or source reaction name.

Recommendation architecture and admission details are documented in the
[type-agnostic reaction recommendation implementation](new/type_agnostic_reaction_recommendation_implementation.md).

## Invariants

- Molecular graphs and normalized edits are the source of truth.
- Structural reconstruction is interpretation-independent.
- The minimum core and generic signature are built before interpretation.
- Interpretation annotations and named-family evidence are optional overlays.
- Ambiguity and conflicts remain typed evidence.
- Rendering explains chemistry but does not define it.
- Conversion and recommendation consume one canonical analysis path.
