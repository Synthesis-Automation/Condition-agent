# Connectivity-First Reaction Grammar Design

## Status

This document is a design proposal. It defines a more uniform foundation for
reaction featurization and grammar execution; it does not by itself change the
current public schemas or reaction behavior.

The proposal extends the graph-first principles in
`type_agnostic_reaction_recommendation_implementation.md`. The molecular graph
remains the source of truth, supplied valid atom mapping remains the strongest
evidence, and named reactions remain optional interpretations.

## 1. Purpose

The central representation of a reaction should be the change in atom
connectivity between supplied reactants and observed products:

- which bond units disappear;
- which bond units appear;
- which atoms gain or lose schema-level hydrogen;
- which atoms and components are retained in the reported product; and
- which unchanged local environments provide context around those edits.

Reaction grammars should be small, declarative graph-rewrite rules. They should
not describe a transition state, catalytic cycle, elementary-step sequence, or
named mechanism.

This design has four goals:

1. Make bond breakage and bond formation the primary reaction identity.
2. Express addition, elimination, substitution, coupling, insertion, exchange,
   rearrangement, ring formation, and fragmentation with the same primitives.
3. Use unchanged structure as local context rather than allowing whole-molecule
   similarity to dominate reaction identity.
4. Preserve observed evidence and uncertainty when products are incomplete,
   atom mapping is absent, regioisomers are possible, or byproducts are omitted.

## 2. Non-goals

The base grammar does not attempt to determine:

- concerted versus stepwise mechanism;
- ionic, radical, organometallic, or pericyclic pathway;
- catalyst turnover sequence;
- kinetic or thermodynamic control;
- Markovnikov, anti-Markovnikov, syn, anti, endo, exo, or facial selectivity
  from reactants alone;
- the identity of an omitted byproduct;
- atom correspondence that is not supported by mapping, exact reconstruction,
  or a conservative correspondence method; or
- a named reaction family from source labels alone.

Those properties may be added as evidence-backed interpretation, condition
context, or stereochemical observation, but they are not part of the primitive
connectivity grammar.

## 3. Core principle: use signed bond deltas

Let the mapped reactant graph be \(G_R\) and the observed product graph be
\(G_P\). For each conserved heavy-atom pair \((i,j)\), define:

\[
\Delta b(i,j) = b_P(i,j) - b_R(i,j)
\]

where a localized covalent bond has an integer bond-unit count:

```text
no bond = 0
single  = 1
double  = 2
triple  = 3
```

The sign has one meaning everywhere:

- `delta_units < 0`: bond unit removal;
- `delta_units > 0`: bond unit addition;
- `delta_units = 0`: unchanged connectivity.

Examples:

| Reactant state | Product state | Canonical delta |
|---|---|---:|
| C–Br | C + Br | -1 |
| C + N | C–N | +1 |
| C=C | C–C | -1 |
| C≡C | C=C | -1 |
| C–C | C=C | +1 |
| C–N | C≡N | +2 |

Calling C=C to C–C a "broken double bond" is potentially misleading because
the C–C edge is retained. The precise statement is that one bond unit is
removed from an edge whose order changes from two to one.

### 3.1 Canonical model

A future canonical contract should be equivalent to:

```python
@dataclass(frozen=True)
class BondDelta:
    atom_1: ReactionAtomReference
    atom_2: ReactionAtomReference
    before_order: str | None
    after_order: str | None
    delta_units: int | None
    evidence: str
    confidence: float
```

The contract should store both states as well as the signed delta. Storing both
allows validation to detect inconsistent programs and avoids relying on the
delta alone for aromatic or nonstandard bonds. `delta_units` is required for
localized single, double, and triple bonds and is `None` for a transition whose
bond domain does not support integer arithmetic.

### 3.2 Compatibility with `ReactionEdit`

The existing `ReactionEdit` representation remains useful for serialization
and display. It can be derived from `BondDelta`:

```text
0 -> n  => formed
n -> 0  => broken
n -> m  => order_changed
```

`formed`, `broken`, and `order_changed` should therefore become compatibility
views, not three independent primitive operations. This removes the current
conceptual split in which addition uses an order-change operator while coupling
uses bond-breaking operators.

### 3.3 Hydrogen changes

Most reaction SMILES omit explicit hydrogen atoms. Hydrogen changes must remain
a separate schema-level contract:

```python
@dataclass(frozen=True)
class HydrogenDelta:
    atom: ReactionAtomReference
    before_count: int
    after_count: int
    delta_count: int
    evidence: str
    confidence: float
```

A grammar may use a virtual hydrogen endpoint during product construction, but
that endpoint is bookkeeping, not an invented mapped atom and not a claim that
one particular hydrogen atom was transferred through a particular mechanism.

Examples:

- N–H consumption: `HydrogenDelta(N, 1, 0, -1)`;
- alkene endpoint hydrogen gain: `HydrogenDelta(C, 1, 2, +1)`;
- alcohol oxidation: hydrogen loss at carbon and possibly oxygen;
- carbonyl reduction: hydrogen gain at carbon and oxygen, when supported by the
  observed product.

### 3.4 Orthogonal atom-state and stereo changes

Formal charge, radical state, isotope, and stereochemistry are not covalent
connectivity. They should remain orthogonal observations:

- `AtomStateDelta` may later represent charge or radical changes that are
  present in the supplied structures.
- `ReactionStereoChange` continues to represent created, destroyed, retained,
  or descriptor-changed stereochemistry.

Neither should be silently inferred from a bond grammar. Reactions with no
connectivity change but a real atom-state change should not be forced into a
bond-edit class.

### 3.5 Aromatic bonds

Aromatic bond arithmetic must not depend on an arbitrary Kekulé form.

- Ordinary `delta_units` arithmetic applies directly to localized single,
  double, and triple bonds.
- Aromatic substitution normally changes a substituent bond, not an aromatic
  ring bond.
- Dearomatization, aromatization, and aromatic cycloaddition require a
  Kekulé-invariant normalized representation or an explicit set of localized
  before/after states.
- The executor must not treat `AROMATIC` as a fractional bond order such as
  1.5 and perform ordinary integer arithmetic on it.

Until a validated aromatic normalization contract exists, broad grammars that
modify aromatic ring bonds should require observed mapping or exact product
verification.

## 4. The reaction edit graph

The normalized deltas form a small signed, atom-provenanced edit graph:

- nodes are atoms incident to at least one nonzero edit;
- negative edges are removed bond units;
- positive edges are added bond units;
- a before/after pair records whether an edge disappears, appears, or changes
  multiplicity;
- hydrogen deltas are signed annotations on nodes;
- atom-state and stereo changes are separate annotations.

This edit graph is the primary reaction center. Its canonical ordering must be
invariant to:

- reactant component order;
- atom serialization order;
- grammar role names;
- display labels;
- source reaction names; and
- symmetric endpoint enumeration.

Atom labels at different signature levels can contain different detail:

```text
base edit shape:       unlabeled incidence topology and delta signs
bond-edit identity:    element pair, before order, after order
handle identity:       reactive-link and site families
exact identity:        local environments and atom provenance
```

## 5. Product scope and atom conservation

Reaction datasets often contain the desired product but omit salts, leaving
group products, water, hydrogen halide, spent organometallic carriers, and
other byproducts. The system must separate three concepts:

1. **Observed connectivity delta:** supported by mapped reactant and product
   atoms.
2. **Grammar reconstruction:** a proposed main product generated from supplied
   reactant sites.
3. **Product projection:** which generated connected components are expected to
   correspond to the reported main product.

Breaking C–Br and N–H and forming C–N does not authorize the system to report
H–Br formation unless HBr is an observed product. A grammar may identify Br
and H as discardable endpoints for main-product reconstruction, but it must not
serialize their hypothetical recombination as observed evidence.

### 5.1 Atom accounting rules

- Every heavy atom in a reconstructed main product must originate from a
  supplied participating component.
- Product heavy-atom excess blocks exact reconstruction and signature
  admission.
- Reactant heavy-atom excess is allowed because omitted byproducts are common.
- Missing reactants, reagent equivalents, carbon atoms, heteroatoms, or
  hydrogen sources must never be synthesized to force a grammar match.
- A disappearing reactant fragment is a product-scope observation, not proof of
  its chemical fate.

### 5.2 Product projection in a grammar

A grammar needs declarative product-retention information separate from its
bond deltas:

```text
retained seeds:
  atom roles that must occur in the reconstructed main product

discardable regions:
  fragments disconnected from every retained seed after the rewrite

product component policy:
  one component, several specified components, or exact supplied product set
```

Projection affects forward reconstruction only. It does not participate in the
canonical observed edit identity except through completeness and provenance.

## 6. Reactant featurization around editable links and endpoints

Molecule featurization cannot know which bond actually reacts without product
evidence. It should therefore expose edit-ready observations rather than
declare a reaction.

### 6.1 Normalized site contracts

Existing chemically useful detectors can be normalized behind three small
contracts.

#### `ReactiveLinkSite`

An existing explicit or schema-level link that a grammar may consume:

```text
site_id
endpoint_a: atom endpoint or virtual hydrogen endpoint
endpoint_b: atom endpoint or virtual hydrogen endpoint
before_order
available_units
source_kind: explicit_bond or implicit_hydrogen
source component and carrier provenance
endpoint contexts
endpoint retention hints
availability
symmetry class
```

Examples:

- C–Cl, C–Br, C–I;
- C–B, C–Sn, C–Zn, C–Si;
- N–H, O–H, S–H, B–H, Si–H;
- H–H, X–X, B–B, Si–B; and
- strained C–C or C–O ring bonds when a curated detector supports them.

`retention hints` describe graph construction, not mechanism:

- the carbon of C–B may be retained while the boron carrier is discardable;
- both Br atoms in Br–Br addition may be retained;
- the hydrogen endpoint of Si–H is virtual;
- both atoms of a multiple-bond acceptor remain retained when its order is
  reduced.

A virtual hydrogen endpoint retains its source component and optional carrier
atom for bookkeeping. H–H may therefore expose two hydrogen endpoints even
though neither is a heavy-atom product seed.

#### `BondCapacitySite`

An existing bond whose multiplicity may change:

```text
endpoint_a
endpoint_b
current_order
maximum_decrement
maximum_increment
bond class
ring and aromatic flags
local endpoint contexts
```

Examples include C=C, C≡C, C=O, C=N, N=N, N=O, and selected single bonds that
may be oxidized to multiple bonds. The detector reports capacity and context;
it does not assert that a reaction will occur.

#### `ConnectionEndpointSite`

An atom that can participate in a new connection under a bounded grammar:

```text
atom
current valence and charge
hydrogen count
availability
context
required accompanying release, hydrogen loss, or charge change
```

This contract should not enumerate every atom in every molecule as generically
reactive. Curated site definitions and grammar constraints must bound candidate
generation.

### 6.2 X–H and A–B unification

`pronucleophile_XH` and `addition_donor` are chemically useful detector labels,
but the grammar executor should see a common link interface.

```text
RNH–H       => atom endpoint N + virtual hydrogen endpoint
RO–H        => atom endpoint O + virtual hydrogen endpoint
R3Si–H      => atom endpoint Si + virtual hydrogen endpoint
B–H         => atom endpoint B + virtual hydrogen endpoint
Br–Br       => atom endpoint Br + atom endpoint Br
R2B–BR2     => atom endpoint B + atom endpoint B
R3Si–BR2    => atom endpoint Si + atom endpoint B
```

The same connectivity template can then handle explicit A–B and implicit A–H
links. Chemistry-specific definitions still determine which links are
available and whether their endpoints are transferred or discarded.

### 6.3 Reaction-centered feature layers

Features for recommendation should be organized by distance from the edit
graph:

```text
Layer 0: edit core
  edited atom pairs, before/after orders, hydrogen deltas

Layer 1: reactive endpoint state
  element, charge, hybridization, H count, aromaticity, ring status

Layer 2: local environment
  immediate substituents, conjugation, steric and electronic descriptors

Layer 3: nearby functional context
  functional groups within a bounded graph radius

Layer 4: unchanged scaffold and spectators
  broader context used for ranking and cautions
```

Retrieval and scoring should compare these layers in this order. A large
unchanged scaffold must not outweigh an incompatible edit graph.

## 7. Declarative connectivity grammar

A grammar is a constrained graph-rewrite rule with five independent parts:

1. **Roles:** site and endpoint predicates.
2. **Relationships:** same/different component, tether, and distinct-site
   constraints.
3. **Rewrite:** signed bond and hydrogen deltas.
4. **Enumeration:** allowed endpoint permutations and symmetry collapse.
5. **Projection and verification:** retained product seeds, discardable
   fragments, valence checks, and evidence requirements.

Named-family compatibility and display labels are optional overlays.

### 7.1 One generic executor

The executable registry should converge on one operator:

```text
apply_connectivity_rewrite
```

It should:

1. resolve grammar endpoint selectors to atom indices;
2. validate every before-state;
3. apply all negative deltas;
4. apply all positive deltas;
5. apply explicit schema-level hydrogen changes;
6. project the permitted product component or components;
7. sanitize and validate valence;
8. canonicalize product and stereochemistry without inventing stereo;
9. emit normalized predicted edits; and
10. return all distinct constitutional outcomes.

Convenience names such as `pair_addition` or `center_replacement` may remain as
declarative rewrite-template IDs, but they should expand to delta programs.
They should not require separate product-building implementations.

### 7.2 Proposed grammar shape

An illustrative definition is:

```json
{
  "id": "generic_pair_addition",
  "roles": {
    "acceptor": {
      "site_type": "bond_capacity",
      "minimum_order": "DOUBLE",
      "decrement_units": 1
    },
    "donor": {
      "site_type": "reactive_link",
      "available_units": 1,
      "transfer_endpoint_count": 2
    }
  },
  "role_relationships": [
    {
      "roles": ["acceptor", "donor"],
      "component_relation": "same_or_different"
    }
  ],
  "rewrite": {
    "bond_deltas": [
      {
        "endpoints": ["acceptor.a", "acceptor.b"],
        "delta_units": -1
      },
      {
        "endpoints": ["donor.a", "donor.b"],
        "delta_units": -1
      },
      {
        "endpoints": ["acceptor.a", "donor.a"],
        "delta_units": 1
      },
      {
        "endpoints": ["acceptor.b", "donor.b"],
        "delta_units": 1
      }
    ],
    "permutations": [
      ["acceptor.a", "acceptor.b"]
    ],
    "retained_seeds": [
      "acceptor.a",
      "acceptor.b",
      "donor.real_atom_endpoints"
    ]
  },
  "requires_product_verification": true
}
```

When `donor.b` is virtual hydrogen, compilation converts the two deltas
incident to it into hydrogen loss and gain. It does not create a heavy-atom
bond to hydrogen in the public edit schema, and a virtual endpoint is never
used as a heavy-atom product-retention seed.

### 7.3 Grammar matching modes

The same grammar should operate in two directions.

#### Observed-edit matching

For a valid mapped reaction:

- derive the edit graph directly;
- match grammar rewrite shapes against observed deltas;
- use matching grammars only as interpretations;
- retain the observed signature even if no grammar matches; and
- report conflicts when a source label or candidate grammar contradicts the
  observed edit graph.

#### Forward reconstruction

For an unmapped reaction:

- enumerate grammar-compatible sites;
- apply the declarative rewrite;
- enumerate symmetry-distinct outcomes;
- compare reconstructed main products with supplied products; and
- select only a unique, exact, chemistry-valid outcome.

Reactant-only plausibility may produce candidates, but not observed reaction
facts.

## 8. Fundamental rewrite templates

The primitive layer needs only signed deltas. A small macro vocabulary is still
useful for authoring and validation.

### 8.1 `connect`

```text
negative deltas: none
positive deltas: A-B
```

Examples:

- association or direct bond formation from two open-valence/charged endpoints;
- intramolecular ring closure when no explicit leaving link is represented;
- recombination in a fully specified mapped record.

Pure connection grammars must be conservative because reactant-only valence
availability is often insufficient to identify a unique product.

### 8.2 `cleave`

```text
negative deltas: A-B
positive deltas: none
```

Examples:

- deprotection with a departing fragment;
- hydrogenolytic or reductive cleavage when hydrogen changes cap the endpoints;
- ring opening;
- fragmentation;
- dealkylation.

The absence of a formed heavy bond does not make the reaction uninformative.
Hydrogen changes and product projection are often essential.

### 8.3 `release_and_connect`

Release one or both retained endpoints from carrier links and connect them:

```text
A-L  -> A + L
B-M  -> B + M       optional when B is already an open endpoint
A + B -> A-B
```

This is the common base logic for:

- alkyl and aryl C–N/C–O/C–S replacement;
- C–C transfer coupling;
- amide, ester, thioester, sulfonamide, and related condensation products;
- oxidative C–H coupling when hydrogen loss supplies one or both released
  endpoints;
- many intramolecular cyclizations; and
- protection reactions that replace H with a protecting-group connection.

At the connectivity level, "substitution" and "coupling" are contextual
interpretations of this shared release-and-connect shape. The base grammar
must not claim SN1, SN2, SNAr, reductive elimination, or another mechanism.

### 8.4 `split_and_distribute`

Consume one bond unit from a retained multiple bond, split a donor link, and
connect the two donor endpoints across the acceptor endpoints:

```text
U=V   -> U-V       remove one unit
A-B   -> A + B     remove one unit
U-A               add one unit
V-B               add one unit
```

This is the uniform addition template.

### 8.5 `depart_and_unsaturate`

Remove endpoint substituent links and add a bond unit between already adjacent
atoms:

```text
U-A   -> U + A
V-B   -> V + B
U-V   -> U=V
```

This covers beta elimination, dehydration, dehydrohalogenation,
dehydrogenation, and related net eliminations. If A–B is not an observed
product, its formation is not invented.

### 8.6 `insert`

Replace one link by two links through a retained atom or fragment:

```text
A-B -> A-I-B

negative: A-B
positive: A-I, I-B
```

This topology supports atom or fragment insertion, including selected oxygen,
nitrogen, carbene, silylene, sulfur, or carbonyl insertions when all product
atoms have supplied provenance.

### 8.7 `extrude`

The inverse connectivity shape removes an intervening atom or fragment and
connects its former neighbors:

```text
A-I-B -> A-B

negative: A-I, I-B
positive: A-B
```

Examples include decarboxylative or denitrogenative net transformations when
the observed atom accounting supports the lost fragment. The grammar describes
only the net graph; it does not claim how extrusion occurred.

### 8.8 `exchange`

Two or more existing links are replaced by crossed links:

```text
A-B + C-D -> A-C + B-D
```

This can represent:

- sigma-bond exchange;
- olefin or alkyne metathesis when bond multiplicities are included;
- disulfide exchange;
- transesterification or transamidation as a more specific
  release-and-connect/exchange case;
- group-transfer reactions; and
- olefination-like product rewiring when product projection discards the
  heteroatom/carrier fragment.

Endpoint retention and product scope distinguish exchange from addition.

### 8.9 `migrate`

A retained atom or fragment changes attachment:

```text
A-M -> B-M

negative: A-M
positive: B-M
```

Additional bond-order shifts may accompany the migration. This covers the
connectivity of many rearrangements without assigning a cationic, radical,
anionic, or concerted pathway.

### 8.10 `annulate`

Several bond deltas jointly create or remove a ring system. Annulation is not a
new primitive; it is a connected multi-edge edit graph plus
`ReactionTopology`:

- two or more formed intramolecular/intermolecular connections;
- zero or more consumed bond units;
- positive cycle-rank delta or observed ring-size change.

Cycloaddition, electrocyclization, and cascade ring formation should use
explicit multi-delta programs rather than mechanism-named operators.

## 9. Deriving structural archetypes

`edit_archetype` should be derived from the normalized edit graph. A grammar
may declare an expected structural shape for validation, but declared metadata
must not override observed deltas.

Counts alone are insufficient. For example, explicit A–B addition and a
two-bond exchange may both contain two negative and two positive edges. Their
edit-incidence topology differs.

The classifier should examine:

- number and sign of deltas;
- whether negative and positive edges share atoms;
- whether a changed atom pair remains adjacent after a negative delta;
- whether a positive delta strengthens an already adjacent pair;
- which edited atoms are retained in the reported product;
- whether a new atom or fragment lies between formerly adjacent atoms;
- ring cycle-rank change;
- connected components of the edit graph; and
- hydrogen annotations.

A proposed base edit-shape vocabulary is:

```text
connection
cleavage
release_and_connect
split_and_distribute
depart_and_unsaturate
insertion
extrusion
exchange
migration
annulation
bond_multiplicity_change
composite
unresolved
```

These are structural summaries, not primary identity. `substitution`,
`coupling`, `condensation`, `addition`, and `elimination` are useful
chemist-facing interpretations derived from a base shape plus endpoint and
retention context. More than one interpretation may be compatible with the
same base shape. For example, C–N bond formation after loss of C–Br and N–H may
be called substitution or coupling depending on conditions and convention; its
canonical edit graph and its `release_and_connect` base shape are unchanged.

## 10. Expansion across diverse reaction types

The following sections test whether the primitive logic is sufficiently
general.

### 10.1 Substitution and transfer coupling

Examples:

```text
C-X + N-H -> C-N
C-X + O-H -> C-O
C-X + S-H -> C-S
C-X + C-B -> C-C
C-X + C-Zn -> C-C
acyl-Y + N-H -> acyl-N
sulfonyl-Y + N-H -> sulfonyl-N
```

Typical deltas:

- remove the substrate-carrier bond;
- optionally remove a partner-carrier bond or partner hydrogen;
- form the connection between retained endpoints.

The same grammar shape covers inter- and intramolecular cases. Site constraints
distinguish sp3 carbon, aryl, heteroaryl, alkenyl, acyl, sulfonyl, phosphorus,
silicon, and other centers.

### 10.2 Addition to C=C and C≡C

Examples:

```text
C=C + H-H
C=C + H-X
C=C + X-X
C=C + O-H / N-H / S-H
C=C + Si-H / B-H
C≡C + A-B
```

For one equivalent of addition:

- decrement the C–C bond order by one;
- consume the A–B or A–H link;
- connect one addend to each endpoint;
- enumerate both orientations when the endpoints or addends are
  constitutionally different; and
- require the supplied product to resolve regioselectivity.

Repeated addition to an alkyne is represented as two events or a two-unit
program only when sufficient reactant equivalents and product evidence exist.

### 10.3 Addition to polarized multiple bonds

The same topology applies to:

```text
C=O + H-H
C=O + C-M
C=O + N-H / O-H / S-H
C=N + H-H
C=N + C-M
N=N + H-H
```

The acceptor is a general `BondCapacitySite`, not only an alkene. Element,
charge, valence, and product-state constraints determine whether an outcome is
valid.

For organometallic addition to C=O:

- decrement C=O by one unit;
- remove C–M from the transferred carbon when represented;
- form carbonyl-C to transferred-C;
- record oxygen charge or hydrogen state only when observed.

The workup proton source must not be invented. A neutral alcohol product
without a supplied/provenanced hydrogen source may still have exact heavy-atom
connectivity but should carry explicit hydrogen-provenance uncertainty.

### 10.4 Reduction and oxidation

Reduction and oxidation are not primitive reaction classes. They are patterns
of bond-unit and hydrogen deltas.

Examples:

```text
C=C -> C-C       bond delta -1, endpoint H gains
C=O -> C-O       bond delta -1, H/state changes
N=N -> N-N       bond delta -1
C-O -> C=O       bond delta +1, H losses
C-C -> C=C       bond delta +1, H losses
C=N -> C#N       bond delta +1, H/state changes
```

Condition interpretation may call these hydrogenation, reduction, oxidation,
dehydrogenation, or aromatization, but the signature remains based on observed
deltas.

### 10.5 Elimination

Examples:

```text
C(X)-C(H) -> C=C
C(OH)-C(H) -> C=C
C(H)-C(H) -> C=C
C(X)=C(H) -> C#C
```

Typical deltas:

- remove a leaving-group connection;
- remove one or more schema-level hydrogens;
- increment the bond order between adjacent retained atoms.

Regioisomeric beta sites are separate outcomes. Stereochemistry is accepted
only from observed product evidence or a validated stereospecific input
contract.

### 10.6 Condensation and olefination

Examples include imine formation, hydrazone formation, oxime formation, and
carbonyl olefination.

Their observed main-product deltas may contain:

- complete loss of a carbonyl C–O edge;
- formation of C=N or C=C with one or more bond units;
- loss of N–H, O–H, or carrier connections;
- disappearance of oxygen or a carrier fragment from the reported main
  product.

These should be explicit exchange or release-and-connect programs with product
projection. They should not be approximated as a single C–N or C–C formation
when the observed bond multiplicity and lost carbonyl oxygen are known.

### 10.7 Insertion and extrusion

Examples:

```text
A-B -> A-O-B
A-B -> A-N-B
A-B -> A-C(=O)-B
A-B -> A-CR2-B
A-I-B -> A-B
```

The inserted atom or fragment must have supplied product provenance. If an atom
appears only in the product, the reaction remains incomplete rather than
assuming a reagent.

Organometallic M–C connectivity requires an explicit bond-domain policy.
Ordinary covalent, dative, ionic, and coordination bonds must not be mixed
without normalization.

### 10.8 Rearrangement and migration

Examples include attachment migration, ring expansion/contraction, and
skeletal rearrangement.

Base patterns include:

- one broken and one formed bond involving the migrated atom/group;
- accompanying bond-order shifts;
- insertion/extrusion-like changes;
- ring cycle-rank or ring-size changes.

The edit graph can represent Wagner–Meerwein-like, allylic, sigmatropic, or
other migration connectivity without claiming the mechanism or the identity of
an unobserved intermediate.

### 10.9 Ring closure and ring opening

Ring closure is not a separate family grammar. It is a rewrite with:

- one or more intramolecular formed connections;
- a positive graph cycle-rank delta; and
- a derived formed-ring size.

Ring opening has a negative cycle-rank delta and usually contains:

- cleavage of a ring bond;
- capping or substitution at one or both endpoints; or
- bond-order redistribution.

The same rewrite template should support intermolecular connection and
intramolecular closure; `ReactionTopology` supplies the distinction.

### 10.10 Cycloaddition and electrocyclic transformations

These are connected multi-bond rewrite programs.

For a simple [2+2] connectivity pattern:

- decrement one bond unit on each alkene;
- form two cross-connections;
- record the new ring topology.

For a Diels–Alder connectivity pattern:

- decrement the two diene double bonds;
- decrement the dienophile double bond;
- increment the internal diene bond;
- form two new terminal cross-connections; and
- record the new six-membered ring.

This says nothing about concertedness, orbital symmetry, endo preference, or
facial selectivity. Those are separate annotations or product observations.

Electrocyclization similarly combines a new ring connection with several
bond-order shifts. It should be represented as one connected edit event, not as
several independent reactions.

### 10.11 Metathesis and bond exchange

Olefin metathesis can be represented as:

- removal of two C=C edges or the required bond units;
- formation of two new C=C edges with crossed endpoint pairing; and
- product projection when only one product is reported.

The grammar must enumerate symmetry-distinct endpoint pairings and require
product verification. It must not infer E/Z geometry.

Sigma-bond exchange and disulfide exchange use the same crossed-link topology
with different site and bond-order constraints.

### 10.12 Fragmentation, deprotection, and decarboxylation

These may contain more negative than positive bond deltas in the reported main
product:

```text
protected-X-R -> X-H
R-C(=O)O group loss -> R-H or R-new_connection
strained ring -> fragments
N-N or O-O cleavage -> capped products
```

The system must permit valid cleavage-centered signatures. It must not require
every accepted reaction to form a heavy-atom bond.

When a fragment is absent from the product side, the completeness assessment
records reactant excess. A label such as "decarboxylation" requires structural
evidence for loss of the corresponding carbon-containing fragment, not merely
a source reaction name.

### 10.13 Protection and deprotection

Protection is generally release-and-connect:

- lose X–H;
- lose protecting-group carrier–handle;
- form X–protecting-group.

Deprotection is generally cleavage plus X–H gain. The base edit graph is
mechanism-neutral; "Boc protection", "benzyl deprotection", and similar names
are contextual annotations based on the retained and departing fragments.

### 10.14 Polymerization and repeated edits

A finite observed oligomer or polymerization step can be represented as
multiple connection/addition events when atom correspondence and stoichiometry
are supplied.

The system should not infer an arbitrary repeat count, duplicate monomer
components, or create a polymer from a single monomer instance. Bulk
polymerization conditions may eventually require a separate repeat-unit
contract layered over the same local edit grammar.

### 10.15 Multi-event and cascade records

Disconnected signed edit graphs are separate `ReactionEvent` objects.
Connected multi-edge rewrites such as cycloaddition or rearrangement remain one
event.

A record with several events does not establish their temporal order.
Sequential, tandem, cascade, and one-pot labels require external evidence. The
base schema stores:

- the canonical event multiset;
- shared atoms and components;
- topological relations; and
- no invented sequence.

## 11. Difficult and excluded cases

### 11.1 Tautomerism and resonance

Different tautomer or resonance representations can create apparent
bond-order and hydrogen edits that are not the intended synthetic
transformation. Standardization must be explicit and versioned. The system
should:

- distinguish input normalization from reaction edits;
- avoid collapsing chemically distinct tautomers without policy;
- recognize representation-only resonance changes where possible; and
- retain uncertainty when normalization changes the apparent center.

### 11.2 Coordination and ionic association

Metal coordination, ion pairing, and salts are not ordinary covalent
connectivity. They require typed bond domains or condition/speciation models.
The covalent grammar should not turn every metal-ligand association into a
formed single bond.

### 11.3 Reactions with missing participating reagents

If the product contains an atom not present among supplied participating
components, a grammar cannot reconstruct it exactly. Agents may be promoted to
participants only through an explicit, typed policy; they must not
automatically be treated as atom donors.

### 11.4 Pure proton transfer or charge-state change

These can have no heavy-atom connectivity edit. They should use hydrogen or
atom-state deltas and remain outside heavy-bond retrieval unless the
recommendation system explicitly supports acid-base or salt-forming processes.

### 11.5 Isotopic exchange

Isotope changes require atom identity and isotope-state deltas. They must not be
collapsed into an ordinary hydrogen-count change when isotope provenance
matters.

## 12. Signature and retrieval consequences

The canonical signed edit graph should become the chemistry core of the
reaction signature.

A proposed hierarchy is:

```text
L0 exact:
  canonical signed edit graph
  atom-local environments
  reactive-link identities
  topology
  stereo observations

L1 handle:
  canonical signed edit graph
  normalized link/site families
  endpoint contexts
  topology

L2 structural rewrite:
  element-labeled edit graph
  derived rewrite shape
  topology
  hydrogen deltas

L3 bond delta:
  before/after element-pair bond states
  hydrogen delta types
  no family requirement

L4 environment:
  local and spectator context for reranking, not reaction identity
```

The exact ordering can preserve the current public L0-L4 names during
migration, but the bond-edit key should ultimately hash canonical
before/after bond deltas rather than separately assembled formed, broken, and
order-change lists.

### 12.1 Context weighting

Recommendation should apply:

1. hard edit compatibility;
2. hard valence, handle, and condition compatibility;
3. rewrite-shape and topology compatibility;
4. reactive endpoint environment similarity;
5. nearby functional-group and spectator compatibility; and
6. broader scaffold similarity.

This order prevents a close scaffold match with the wrong reaction center from
outranking a less similar scaffold with the correct bond rewrite.

## 13. Evidence and conflict policy

Evidence priority remains:

1. validated supplied atom mapping and observed deltas;
2. exact product reconstruction from a declarative rewrite;
3. exact composition of multiple rewrites;
4. conservative unique atom correspondence;
5. unresolved or conflicting candidates.

Grammar metadata never overrides observed edits.

Conflicts should be typed, for example:

```text
GRAMMAR_DELTA_CONFLICT
SOURCE_FAMILY_EDIT_CONFLICT
PRODUCT_PROJECTION_CONFLICT
HYDROGEN_PROVENANCE_UNRESOLVED
AROMATIC_NORMALIZATION_UNRESOLVED
MULTIPLE_EXACT_REWRITE_OUTCOMES
```

Broad grammars should require product verification. Reactants alone generally
cannot resolve regioselectivity, stereochemistry, chemoselectivity, or which of
several compatible sites participated.

## 14. Validation requirements

Definition loading should validate:

- every endpoint selector resolves to a declared role and atom role;
- every negative delta has a compatible reactant before-state;
- every positive delta has a valid product target state;
- no pair receives contradictory deltas;
- hydrogen deltas reference valid carrier atoms;
- product-retained seeds cannot be discarded;
- permutation groups contain compatible endpoint types;
- declared expected shape agrees with the compiled delta program;
- aromatic edits use an approved normalization mode;
- all generated products pass RDKit sanitization and explicit valence checks;
  and
- executable behavior is selected only from a fixed Python registry.

No arbitrary Python callable may be named or imported from a JSON definition.

## 15. Migration plan

### Phase 1: canonical delta view

- Add typed `BondDelta` and `HydrogenDelta`.
- Convert normalized `ReactionEdit` objects to and from the delta view.
- Prove deterministic ordering and partner-order invariance.
- Add the delta representation without changing current signature identity.

### Phase 2: generic executor

- Add `apply_connectivity_rewrite`.
- Compile an existing addition grammar to a declarative delta program.
- Compile an existing C–C coupling grammar to the same executor.
- Compare predicted products and edits with current operators.

### Phase 3: normalized site interfaces

- Introduce `ReactiveLinkSite`, `BondCapacitySite`, and
  `ConnectionEndpointSite` interfaces.
- Adapt existing leaving-group, transfer-group, X–H, addition-donor, and
  unsaturated-bond detectors.
- Preserve existing detector labels as chemistry annotations.

### Phase 4: operator migration

- Migrate `center_replacement`.
- Migrate `join_two_anchors`.
- Migrate `pair_addition` and `pair_elimination`.
- Migrate simple bond-order-change grammars.
- Migrate bespoke operators only after their complete delta and projection
  semantics are explicit.
- Remove obsolete executable paths after parity tests pass.

### Phase 5: derived edit shapes

- Replace grammar-declared archetype authority with edit-graph inference.
- Keep a declared expected shape only as grammar validation metadata.
- Add coupling, insertion, extrusion, exchange, migration, annulation, and
  cleavage classifications where the topology is unambiguous.
- Preserve `unresolved` rather than forcing a class.

### Phase 6: signature migration

- Add canonical edit-graph tokens to L0-L3.
- Measure collision and retrieval changes on existing datasets.
- Bump schema and definition versions.
- Update converted artifacts only after chemistry and retrieval gates pass.

## 16. Required test matrix

Every migrated template requires positive, negative, ambiguous, conflict, and
partner-order-invariance tests.

Minimum chemistry coverage:

### Release and connect

- Suzuki C–C;
- C–N, C–O, and C–S;
- alkyl and aryl contexts;
- acyl N/O/S formation;
- intramolecular ring closure;
- missing carrier atom in reactants;
- two equally valid reactive sites.

### Addition

- H2, HX, X2, N–H, O–H, S–H, Si–H, and B–H;
- alkene and one-equivalent alkyne addition;
- polarized C=O and C=N examples;
- symmetric and asymmetric donor;
- symmetric and asymmetric acceptor;
- both regioisomeric outcomes;
- explicit stereoproduct conflict;
- no available bond unit.

### Elimination and multiplicity changes

- beta-halo elimination;
- dehydration;
- dehydrogenation;
- carbonyl reduction and alcohol oxidation;
- alkyne formation;
- no beta hydrogen;
- multiple beta sites;
- aromatic edit rejection or approved normalization.

### Complex rewrites

- insertion and extrusion;
- attachment migration;
- ring opening and closure;
- [2+2] and Diels–Alder connectivity;
- olefin exchange/metathesis;
- condensation or olefination with omitted byproduct;
- disconnected multi-event reaction;
- connected multi-edit rearrangement.

### Evidence and serialization

- valid mapped unknown-family reaction;
- invalid and partial atom maps;
- exact mapped/operator agreement;
- mapped/operator conflict;
- omitted byproducts;
- product atom excess;
- deterministic IDs;
- reactant component-order invariance;
- symmetric outcome collapse;
- signature schema round trip.

The complete repository test suite must pass before any operator is removed.

## 17. Acceptance criteria

The design is successfully implemented when:

1. Addition to C=C and C≡C is represented as removal of one acceptor bond unit,
   consumption of one donor link, and formation of two endpoint connections.
2. X–H and explicit A–B donors use one normalized link interface.
3. C–C, C–N, C–O, and C–S coupling use the same generic rewrite executor.
4. Existing `formed`, `broken`, and `order_changed` outputs are deterministically
   derived from canonical deltas.
5. Mapped unknown reactions receive valid signatures without a grammar or named
   family.
6. Broad grammars do not force regioselectivity, stereochemistry, byproducts,
   or atom correspondence.
7. Reaction retrieval is keyed first by the signed edit graph and only then by
   unchanged molecular context.
8. Specialized executable operators are reduced to a small, validated generic
   registry.
9. Aromatic, coordination, missing-reagent, and hydrogen-provenance limitations
   remain explicit.
10. Chemistry parity and full-suite regression tests pass.

## 18. Final architectural rule

The system should answer these questions in order:

```text
1. Which atom pairs changed bond state?
2. Which atoms changed schema-level hydrogen state?
3. What connected edit graph do those changes form?
4. Which reactant links and endpoints can express that graph?
5. What unchanged local environments modify compatibility and conditions?
6. Which optional structural class or named family is supported?
```

The governing rule is:

> A reaction is an evidence-backed graph rewrite. Grammars constrain and
> execute rewiring of supplied atoms. Structural archetypes, named reactions,
> mechanisms, and unchanged scaffold features are layered interpretations or
> context, never substitutes for the observed connectivity change.
