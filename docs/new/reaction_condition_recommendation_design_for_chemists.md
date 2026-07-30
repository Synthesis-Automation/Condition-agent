# Reaction-Condition Recommendation: A Chemist's Design Guide

**Status:** Consolidated current design  
**Audience:** Synthetic, medicinal, process, and computational chemists  
**Reviewed against the repository:** 2026-07-30  
**Implementation status:** See
[`type_agnostic_reaction_recommendation_implementation.md`](type_agnostic_reaction_recommendation_implementation.md)

## 1. The short version

The system recommends reaction conditions by asking:

1. What changed between the supplied reactants and reported product?
2. Which atoms and reactive handles participated?
3. What local molecular environments may alter reactivity or compatibility?
4. Which structurally compatible literature precedents used the same or a
   closely related transformation?
5. Which complete condition recipes have independent support?

A reaction name such as *Suzuki–Miyaura* is useful when the molecular evidence
supports it, but it is not required. A mapped, chemically valid reaction can be
represented and retrieved even when its named family is unknown.

The governing idea is:

> A reaction is an evidence-backed molecular-graph change. Reaction names,
> mechanisms, and source labels are interpretations of that change, not
> substitutes for it.

This is a recommendation and evidence-ranking system. It does not prove that a
reaction will work, establish a mechanism, or replace experimental risk
assessment.

## 2. Why the system is organized this way

Reaction datasets are chemically valuable but imperfect:

- reaction names may be absent, inconsistent, or wrong;
- atom maps may be missing or invalid;
- only the main product is often reported;
- salts and byproducts are commonly omitted;
- condition columns do not always describe the actual role of a substance;
- one publication may contribute many highly correlated substrate-scope rows;
- a reported yield is an observation under one procedure, not a universal
  property of a recipe.

Routing every record by its source reaction name would discard useful unnamed
chemistry and allow incorrect labels to override structures. The system instead
separates three layers:

| Layer | Chemical question | Examples |
| --- | --- | --- |
| Observation | What is present and what changed? | C–Br loss, C–N formation, N–H loss, local pyridine ring, retained nitrile |
| Interpretation | What structural class or named family is supported? | substitution; `sp2_c_n_substitution`; possible Buchwald–Hartwig |
| Recommendation | Which compatible precedents and recipes are sufficiently supported? | same edits and handles, compatible spectators, two independent references |

If interpretation fails, valid observations remain available. If the evidence
is inadequate, the recommender abstains or uses a separately labelled,
lower-confidence fallback.

## 3. End-to-end chemical workflow

```text
reaction SMILES
    |
    v
parse reactants, agents, and products
    |
    v
detect molecular features and canonical reactive sites
    |
    v
obtain atom correspondence or exactly reconstruct the product
    |
    v
normalize bond, hydrogen, charge, topology, and stereo observations
    |
    +----> if correspondence is ambiguous, retain every distinct edit
    |      hypothesis without declaring one to be observed
    |
    v
build a versioned reaction signature
    |
    +----> optional transformation and named-family interpretation
    |
    v
resolve condition substances, roles, stages, and recipe identity
    |
    v
admit, quarantine, or reject each precedent
    |
    v
retrieve the narrowest adequately supported compatible precedent pool
    |
    v
apply hard chemistry filters, score similarity, and aggregate recipes
    |
    v
recommendations, precedents, explanations, cautions, or abstention
```

### 3.1 Parse the reported reaction

The normal query is a product-specified reaction:

```text
reactant1.reactant2>>product
```

The three-part form `reactants>agents>products` is also accepted. Original
SMILES, component identity, component order, atom-map numbers, parsing
warnings, and canonical unmapped structures are retained.

Middle-side agents are analyzed but are not automatically treated as atom
donors. Doing so would invent participating reactants. They may provide
explicitly labelled support for an unresolved product-fragment source, but
they do not create supplied atom correspondence.

### 3.2 Describe each molecule without assuming a reaction

Molecular featurization detects functional groups and candidate reactive
features such as:

- C–halogen and other releasable links;
- N–H, O–H, S–H, B–H, Si–H, and related X–H links;
- organoboron and other transfer-carrier links;
- carbonyl, acyl, sulfonyl, and other activated centers;
- polarized C=N bonds in imines, oximes, and hydrazones;
- explicit C, N, O, and S anions without inferred deprotonation;
- epoxide and aziridine ring-opening carbons;
- conditional O–Si release links in carbon-bound silyl ethers;
- bonded Li, Cu, and hydrogen-free Al carbon-transfer links;
- aromatic C–H sites;
- alkenes, alkynes, and other bonds with localized bond-order capacity.

These observations are exposed to reaction grammars through three generic
interfaces:

| Interface | Meaning | Examples |
| --- | --- | --- |
| Reactive link | An explicit A–B bond or schema-level A–H link that may be consumed | C–Br, C–B, N–H, Br–Br |
| Bond capacity | A localized bond whose order may change | C=C, C≡C, C=O |
| Connection endpoint | An atom that may form a new bond under stated valence/release requirements | amine N, aryl C, acyl C |

The presence of a handle is not evidence that it reacted. Product evidence
selects the participating site.

Newly detected sites do not automatically create a reaction grammar. For
example, detecting both an epoxide carbon and an alkoxide only records two
possible reactive features. A verified ring-opening interpretation still
requires product mapping, correspondence, or a separately validated bounded
rewrite.

### 3.3 Establish what changed

The evidence order is:

1. validated supplied atom mapping;
2. exact product reconstruction by one bounded connectivity rewrite;
3. exact composition of multiple rewrites;
4. conservative unique-scaffold correspondence;
5. bounded global correspondence across supplied reactants;
6. bounded fragmented-scaffold correspondence for topology-changing,
   single-substrate reactions;
7. unresolved, ambiguous, or conflicting evidence.

Fragmented-scaffold correspondence is useful when cyclization or atom loss
splits the conserved graph into pieces that an ordinary connected scaffold
match cannot span. It enumerates symmetry-equivalent alternatives and accepts
only one normalized edit interpretation. Close alternatives with different
edits are reported as ambiguous rather than resolved by the reaction name.

The system records:

- formed and cleaved bonds;
- bond-order changes;
- schema-level hydrogen gains and losses;
- observed formal-charge changes;
- retained, inverted, created, destroyed, or changed explicit stereo;
- whether edits form one event or several disconnected events;
- intermolecular, intramolecular, mixed, or unimolecular topology;
- ring formation and inferred ring size where supported.

Hydrogen changes describe a count at an atom; they do not invent the identity
or mechanistic path of a particular hydrogen atom.

For an omitted leaving-group byproduct, disappearance from the reported main
product is a main-product projection. It must not be overstated as an observed
free leaving group or a fully balanced equation.

### 3.4 Reconcile evidence rather than choosing a convenient answer

Valid mapping is the strongest observation, while exact reconstruction is
strong independent support. When they agree, confidence increases. When they
disagree, the mapped observation is retained with explicit conflict evidence
and the record is routed to review. A reaction grammar or family label never
silently overrides contradictory graph evidence.

Candidate rewrites may enumerate several regioisomers or stereochemical
outcomes. Only the reported product may select among them. Reactants alone do
not justify Markovnikov/anti-Markovnikov, syn/anti, E/Z, facial, or mechanistic
claims.

### 3.5 Check product-atom completeness

Each reaction receives one of three product-provenance states:

| State | Meaning | Normal indexing |
| --- | --- | --- |
| Verified | Every reported product heavy atom has supported provenance | May be eligible |
| Incomplete | Product atoms cannot be supplied by the reported reactants | Not eligible for the normal verified index |
| Unresolved | No definite excess, but provenance cannot be established | Review or explicitly labelled fallback only |

Reactant atoms missing from the main product are expected when byproducts are
omitted. Product atoms missing from all supplied reactants are a different
problem: the system never creates the missing reagent or duplicates a partner
to force a match.

For a uniquely supported single-scaffold replacement, an incomplete record may
retain a partial product transformation such as:

```text
Ar–Br -> Ar–C#N  [source of C,N not supplied]
```

This is useful review evidence. It is not promoted to a verified reaction
signature or presented as a fully observed reaction.

## 4. Connectivity grammars: structural, not mechanistic

Reaction grammars are bounded graph-rewrite rules. They declare compatible
sites and the bond, hydrogen, and charge changes needed to construct a possible
product.

Common structural shapes include:

| Structural shape | Net graph description | Chemistry examples |
| --- | --- | --- |
| Release and connect | Release one or two endpoints, then form a new bond | C–N/O/S substitution, Suzuki-type C–C formation, acyl substitution |
| Split and distribute | Lower a multiple-bond order, split A–B or consume A–H, attach A and B across the bond | H₂, HX, X₂, X–H, Si–H, or B–H addition |
| Depart and unsaturate | Lose a substituent and H, increase an adjacent bond order | β-elimination, dehydration |
| Localized bond-state change | Change bond order with supported H/charge changes | hydrogenation, oxidation, reduction |

`substitution`, `addition`, and `elimination` are edit-derived structural
summaries. They do not assert SN1, SN2, radical, ionic, concerted, or
organometallic mechanisms.

Named families are optional overlays. For example:

```text
observed graph change:  sp2 C–Br cleavage + C–C formation + C–B projection
transformation class:   c_c_transfer_coupling
named family:           suzuki_miyaura, if uniquely and exactly supported
```

An alkyl halide plus an alcohol may support `sp3_c_o_substitution` while
leaving `named_family=None`, because connectivity alone does not distinguish
SN1, SN2, Mitsunobu, or protection chemistry.

## 5. Reaction identity and similarity

### 5.1 The reaction signature

Every verified usable reaction is represented by a deterministic, versioned
reaction signature. It contains:

- normalized edits and explicit stereo observations;
- reaction events and topology;
- participating partners and handles;
- local reactivity profiles;
- unchanged spectator groups;
- optional transformation and family interpretations;
- completeness, evidence quality, warnings, and definition versions.

The identifier excludes source reaction names, display wording, row order,
irrelevant reactant ordering, and serialization formatting.

### 5.2 Retrieval levels

The signature provides several views so retrieval can relax in controlled,
explainable steps:

| Level | Main chemical content | Purpose |
| --- | --- | --- |
| L0 exact | Edits, stereo, events, topology, handles, attachment contexts, and detailed profile tokens | Closest structural precedents |
| L1 handles | Edits, stereo, events, topology, and reactive handles | Same reaction core with less local detail |
| L2 transformation | Bond/H changes, events, topology, and transformation class | Same generic transformation |
| L3 bond edits | Bond and H changes without topology | Broad chemistry fallback |
| L4 environment | Reactivity profiles, nearby groups, and spectators | Neighbor selection and reranking |

Topology belongs in L0–L2. L3 can cross an intra/intermolecular boundary only
as a disclosed fallback and must carry an explicit caution.

After exact L3 retrieval, an additional anonymous edit-graph neighbor tier may
relate reactions whose precursor forms produce similar but non-identical net
edits. It requires shared formed-bond chemistry, compatible broken-bond
chemistry, and compatible ring-count direction before similarity is evaluated.
For example, a carbonyl plus arylhydrazine and its preformed hydrazone may enter
an indole-forming transformation through different exact edit sets. This tier
can relate them without using the words “Fischer indole” as a routing key. The
result explicitly states that the exact bond-edit signatures differ.

An ambiguous query does not receive an RS3 signature merely because its
alternatives look plausible. Instead, each distinct edit hypothesis is
converted to an anonymous edit prototype and searched independently against
verified precedents. Recommendations are possible only from the intersection:
the same precedent chemistry must pass the edit and compatibility gates for
every hypothesis. Ranking uses the weakest hypothesis-to-precedent match. This
is useful when ambiguity concerns atom origin or symmetry but all alternatives
still imply the same practical chemistry; otherwise the system abstains.

## 6. Local reactivity profiles

Every selected site has a typed, graph-derived profile with a common outer
organization:

```text
context identity | steric accessibility | electronic activation |
reactive-center state | modifiers and liabilities
```

Context-specific descriptors are used because unrelated chemical quantities
should not be forced onto one axis:

- aromatic systems: ring family, ring sizes, fusion, heteroatom positions,
  ortho occupancy/burden, and graph-derived electron demand;
- alkyl centers: substitution, α/β branching, cyclicity, activation, and β-H
  availability;
- alkenyl and alkynyl centers: endpoint substitution and unsaturated context;
- acyl, sulfonyl, phosphoryl, and related centers: activation and leaving-group
  state;
- N/O/S/P centers: hydrogen count, substitution, conjugation, lone-pair class,
  and availability where meaningful.

These are deterministic 2D molecular-graph descriptors. Electronic demand is
an interpretable chemistry prior, not calculated electron density. Steric
accessibility is a bounded graph proxy, not a conformational Sterimol or buried
volume calculation.

Only stable categorical bins enter reaction identity. Continuous values and
their atom-level contributors remain available for similarity and explanation.
Missing, unresolved, not-applicable, and not-computed states are distinct.

## 7. Condition identity and recipe representation

A source column called `reagent` or `solvent` is only a role hint. Substance
identity is resolved conservatively through the condition registry, and its
role is interpreted in reaction context.

Each component retains:

- the raw source identifier and source field;
- canonical substance identity when resolved;
- contextual role and confidence;
- amount and unit when reported;
- provenance and uncertainty warnings.

Recipes use two related identities:

| Identity | Includes | Chemical use |
| --- | --- | --- |
| `RCORE1` recipe core | Resolved substances and contextual roles | Groups the same underlying condition regime |
| `RCR1` recipe variant | Recipe core plus amounts, temperature, time, concentration, atmosphere, stages, and versions | Preserves reported operating variants |

The recipe vocabulary includes catalysts, ligands, bases, acids, condensation
agents, oxidants, reductants, additives, solvents, and other components.
Unresolved identities remain explicit; they are not guessed or silently
discarded.

Multi-stage procedures remain ordered. If conditions cannot be assigned to the
correct stage or reaction event, the record is quarantined or penalized rather
than treated as a clean single-stage precedent.

## 8. From precedents to recommended recipes

### 8.1 Admission is multidimensional

The converter assesses four independent questions:

- Is the reaction chemistry verified, review-only, or rejected?
- Are condition identities complete, partial, unresolved-but-retained, or
  unusable?
- Is the outcome usable, missing, or invalid?
- Is the record eligible, review-only, or ineligible for indexing?

An unknown named family does not reduce verified chemistry. Missing yield does
not erase a useful chemistry/recipe observation, but that row is excluded from
yield modeling.

### 8.2 Retrieval uses the narrowest adequate pool

The current verified-signature ladder is:

1. exact signature;
2. relaxed handle signature;
3. high-confidence named family, still constrained by compatible edits;
4. generic transformation signature;
5. local-environment neighbors inside the compatible bond-edit pool;
6. broad compatible bond-edit fallback;
7. chemistry-gated anonymous edit-graph neighbors;
8. abstention.

The first tier with enough independent compatible support is selected. Family
evidence never crosses an incompatible-edit gate.

A query without a verified signature does not enter this ladder. If it has
multiple typed edit hypotheses, it may use the separate all-hypothesis
consensus search described above. That search uses only verified precedents,
requires independent support, remains review-required, and does not make the
query itself edit-verified. If consensus is unavailable, a conservative
structure fallback with its own thresholds, trace, and warnings may be tried.

### 8.3 Compatibility precedes similarity

Hard exclusions can include:

- incompatible net edits or event sets;
- impossible valence or incompatible reactive-handle class;
- essential topology conflicts;
- known oxidant/reductant, acid, moisture, or atmosphere conflicts;
- definitely missing mandatory recipe components;
- invalid or conflicting transformation evidence.

Potential coordination, catalyst poisoning, hydrolysis, acid/base, temperature,
or incomplete-recipe risks are generally penalties or cautions when the
evidence is insufficient for a hard exclusion.

Only compatible precedents are scored. Interpretable score components compare
edits, events, topology, handles, local profiles, nearby groups, spectators,
family confidence, and evidence quality. Missing features do not count as
matches.

### 8.4 Recipes are aggregated with correlated evidence controlled

Retrieved rows are grouped by `RCORE1`, with observed `RCR1` variants exposed.
Support is reported separately as:

- raw observations;
- unique canonical reactions;
- reference-local condition series;
- independent publications;
- source datasets.

Many scope examples from one paper show substrate breadth, but they do not
become many independent literature confirmations. Repeated observations within
one publication are deduplicated or receive diminishing support.

Yield is secondary because published reaction corpora are selection-biased.
Expected yield is omitted when there is no usable outcome evidence; it is never
invented.

## 9. What a chemist should expect in the result

A recommendation should state:

- the observed transformation and optional named family;
- the evidence used to establish the reaction;
- the retrieval level and every attempted fallback;
- matching edits, handles, topology, and local environments;
- important mismatches and compatibility cautions;
- the complete resolved recipe and missing fields;
- independent support counts;
- representative precedent reaction and reference IDs;
- score components, definition versions, and uncertainty.

For example, a result may say:

```text
Query:
  HeteroAr–Br + Ar–NH2 -> HeteroAr–NH–Ar
  transformation: sp2 C–N substitution
  family: unresolved

Retrieval:
  L1 handle-compatible precedents; 3 independent references

Important match:
  same C–Br loss, C–N formation, N–H loss, and intermolecular topology

Caution:
  query heteroaryl N may coordinate the catalyst; supporting examples are
  mostly less coordinating carbocyclic aryl bromides
```

The absence of a family name is not an error. A broad fallback, single-paper
support, unresolved condition identity, topology relaxation, or incomplete
recipe should be visible rather than hidden in a score.

## 10. Current scientific limits

The system should abstain or retain review evidence when it cannot support the
claimed chemistry. Important limits include:

- missing reactants or incomplete product-atom provenance;
- ambiguous atom correspondence, regioselectivity, or stereochemistry;
- coordination chemistry and ion pairing not represented by ordinary covalent
  bonds;
- pure proton transfer, isotope exchange, and unsupported radical-state
  changes;
- complex rearrangements, cascades, and pericyclic reactions without validated
  multi-edit contracts;
- unstructured multi-stage procedures;
- incomplete condition-registry coverage;
- limited precedent support for some transformations and environments.

Current calibration and converted corpora are development evidence, not proof
of broad production accuracy. Every proposed recipe still requires chemical
judgment, safety review, and experimental validation.

## 11. Compact glossary

| Term | Meaning |
| --- | --- |
| Reactive handle | A graph-local feature that can participate in a bounded edit, such as C–Br, N–H, C–B, or C=C |
| Reaction grammar | Declarative constraints connecting compatible sites to a bounded graph rewrite |
| Reaction edit | A formed, cleaved, order-changed, or schema-level H change with atom provenance |
| Reaction event | One connected group of edits; a reaction may contain several |
| Topology | Intermolecular/intramolecular scope, tether, and ring-formation information |
| Reaction signature | Versioned, deterministic chemistry identity used by conversion and retrieval |
| Spectator group | An unchanged group retained near or away from the reaction center that may affect compatibility |
| Named family | Optional interpretation such as Suzuki–Miyaura; never the primary structural identity |
| Recipe core | Role-aware condition substances without operating variants |
| Recipe variant | Recipe core plus reported quantities and operating conditions |
| Independent support | Distinct reactions, condition series, references, or corpora rather than raw row count |
| Typed abstention | An explicit result that the available evidence or support is insufficient |
