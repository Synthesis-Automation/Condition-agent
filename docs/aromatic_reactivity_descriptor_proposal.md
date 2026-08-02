# Aromatic Reactivity Descriptor Proposal

The phased implementation and migration sequence is defined in
[`context_aware_reactivity_descriptors_implementation_plan.md`](context_aware_reactivity_descriptors_implementation_plan.md).

## 1. Purpose

The current partner summary mixes reactive-center identity, steric measurements,
electronic heuristics, and low-value implementation details:

```text
steric ortho hindered, 2 ortho substituent(s), 8 local heavy atom(s), r<=2;
electronic neutral (local score +0)
```

This is difficult to compare across partners. In particular:

- `local_heavy_atoms_r2` is a debugging measurement, not an immediately useful
  reactivity explanation;
- `class`, `ortho_substituent_count`, and `local_heavy_atoms_r2` partly describe
  the same steric environment at different abstraction levels;
- the current electronic score is a distance-weighted functional-group tag sum,
  not a physical electron-density calculation;
- `Ar` and `HetAr` omit ring identity, ring size, fusion, and the positions and
  electronic types of ring heteroatoms;
- an aryl electrophile and an aromatic N-H partner are affected by different
  electronic quantities. Ring electron demand is useful for the former, while
  lone-pair availability, resonance delocalization, and acidity are often more
  useful for the latter.

The proposed design creates one uniform, mechanism-neutral observation schema
and a separate role-aware interpretation. The default display shows only the
highest-value factors; complete structured evidence remains available for
retrieval, rules, review, and debugging.

## 2. Design rules

1. **Site relative:** every descriptor is computed relative to the atom that
   participates in the reaction edit.
2. **Graph derived:** ring identity, ring positions, substituents, and local
   environments come from the molecular graph, not the reaction name or source
   label.
3. **Observation before interpretation:** record ring topology and substituent
   effects first; interpret them as electrophile activation or nucleophile
   availability only when the reaction context supports that role.
4. **Uniform structure, role-specific meaning:** every partner has the same
   top-level descriptor groups, but optional role-aware assessments may differ.
5. **No false precision:** use `electronic_demand`, not `electron_density`, for a
   graph-derived heuristic. A future quantum-calculated value must be a separate
   descriptor with its method and provenance.
6. **Categorical identity, continuous support:** serialize both stable categories
   and normalized scores. Signature keys use categories or bins, not raw
   floating-point values.
7. **Concise by default:** raw shell counts, contributing atoms, and intermediate
   scores belong in an expanded evidence view, not the primary reaction summary.
8. **Versioned and deterministic:** vocabulary, algorithms, bin boundaries, and
   rendering rules are declarative and versioned.

## 3. Proposed partner contract

Replace the open-ended meaning of the current `steric` and `electronic`
dictionaries with typed nested descriptors:

```python
@dataclass(frozen=True)
class AromaticHeteroatom:
    element: str
    formal_charge: int
    aromatic_role: Literal["pyridine_like", "pyrrole_like", "other"]
    ring_distance_from_anchor: int
    positional_relation: Literal[
        "anchor", "ortho", "meta", "para", "remote", "fused_other_ring"
    ]


@dataclass(frozen=True)
class AromaticSystemDescriptor:
    system_class: Literal[
        "carbocyclic", "heteroaromatic", "mixed_fused", "not_aromatic"
    ]
    ring_family: str
    ring_sizes: tuple[int, ...]
    aromatic_ring_count: int
    fused: bool
    heteroatoms: tuple[AromaticHeteroatom, ...]
    anchor_in_ring: bool
    confidence: float
    method: str


@dataclass(frozen=True)
class StericProfile:
    accessibility_class: Literal["open", "moderate", "hindered", "severe"]
    accessibility_score: float
    ortho_occupancy_count: int | None
    ortho_capacity: int | None
    ortho_burden_class: Literal["none", "low", "medium", "high"] | None
    ortho_burden_score: float | None
    center_substitution_class: str | None
    attached_group_burdens: tuple[dict[str, Any], ...]
    method: str


@dataclass(frozen=True)
class ElectronicContribution:
    source_id: str
    effect: Literal["withdrawing", "donating", "mixed"]
    pathway: Literal["inductive", "resonance", "charge", "aromatic_intrinsic"]
    positional_relation: str
    contribution: float


@dataclass(frozen=True)
class ElectronicProfile:
    demand_class: Literal[
        "electron_rich", "slightly_rich", "balanced",
        "slightly_poor", "electron_poor"
    ]
    demand_score: float
    contributions: tuple[ElectronicContribution, ...]
    confidence: float
    method: str


@dataclass(frozen=True)
class ReactiveCenterProfile:
    element: str
    hybridization: str
    formal_charge: int
    hydrogen_count: int
    heavy_atom_attachment_count: int
    conjugation_class: str | None
    lone_pair_class: str | None
    lone_pair_availability: Literal["high", "medium", "low", "not_applicable"]
    acidity_class: str | None


@dataclass(frozen=True)
class PartnerReactivityDescriptor:
    aromatic_system: AromaticSystemDescriptor
    steric: StericProfile
    electronic: ElectronicProfile
    reactive_center: ReactiveCenterProfile
    flags: tuple[str, ...]
    definition_version: str
```

The public contracts should use immutable typed dataclasses. Small dictionaries
inside contribution records are acceptable temporarily during migration, but
their keys must have a validated vocabulary and a removal criterion.

## 4. Aromatic-system representation

### 4.1 Ring family

`ring_family` is a graph-derived chemist-facing class such as:

- `benzene`;
- `pyridine`, `pyridazine`, `pyrimidine`, or `pyrazine`;
- `pyrrole`, `furan`, `thiophene`, `imidazole`, or related five-membered rings;
- `naphthalene`, `quinoline`, `indole`, `carbazole`, or another validated fused
  system;
- `other_carbocyclic_aromatic`, `other_heteroaromatic`, or
  `ambiguous_fused_system` when a specific validated class is unavailable.

Specific ring names are annotations derived from graph topology. Generic fields
such as ring sizes, heteroatom composition, and relative positions remain the
machine-readable source of truth.

### 4.2 Heteroatom positions

Do not depend on arbitrary molecule atom numbering. Store heteroatoms relative
to the reactive anchor:

- cyclic distance from the anchor;
- chemist-facing relation (`ortho`, `meta`, `para`, or `remote`) when meaningful;
- whether the heteroatom is in the anchor ring or another fused ring;
- pyridine-like versus pyrrole-like electronic role;
- element and formal charge.

For symmetric rings, normalize the two possible directions and retain the
lexicographically smallest position sequence. For fused systems, use shortest
ring-path distance plus same-ring/fused-ring membership rather than forcing
benzene-style positions.

### 4.3 Ring size and fusion

Store the sorted sizes of aromatic rings in the connected aromatic system, the
number of aromatic rings, and whether the system is fused. This distinguishes,
for example, benzene from pyridine, pyridine from pyrimidine, and pyrrole from
indole without relying only on `Ar` versus `HetAr`.

## 5. Steric descriptor

### 5.1 What should affect the primary summary

For an aromatic reactive site, the primary steric factors are:

1. occupancy of positions immediately adjacent to the reactive anchor;
2. the size and branching of those substituents;
3. ring fusion adjacent to the anchor;
4. for a reactive heteroatom, shielding by directly attached groups;
5. an overall graph-derived accessibility class.

`ortho_occupancy_count` alone is not a total steric descriptor. Two methyl
groups and two tert-butyl groups both have occupancy two but impose very
different burdens. The proposal therefore reports both occupancy and a
size-weighted burden.

### 5.2 Deterministic graph proxy

Use a versioned graph proxy rather than an unexplained atom count:

- identify each approach-side branch adjacent to the reactive center;
- accumulate heavy atoms with distance attenuation over a bounded radius;
- add explicit branching and adjacent-fusion contributions;
- normalize the result to `[0, 1]`;
- bin it into `open`, `moderate`, `hindered`, or `severe`.

Keep the individual branch contributions as evidence. Do not show the raw
radius-shell heavy-atom count in the default UI.

For an N-H partner, separately retain:

- substitution at N (`primary`, `secondary`, and so on);
- burden of each directly attached group;
- total ortho burden on attached aryl systems.

This avoids the misleading combination `secondary center, attached: Ar, Ar`
while still preserving both facts.

## 6. Electronic descriptor

### 6.1 Mechanism-neutral observation

The base electronic descriptor is **electron demand at the reactive anchor**.
Use a signed normalized convention:

```text
-1.0 = strongly electron donating / electron rich
 0.0 = balanced
+1.0 = strongly electron withdrawing / electron poor
```

The score should combine explicit, auditable contributions from:

- formal charge;
- intrinsic aromatic heteroatom effects;
- recognized substituent inductive effects;
- recognized substituent resonance effects at meaningful ring positions;
- attenuation by graph distance;
- protonation or deprotonation state.

Every contribution records its source, pathway, relative position, and signed
value. A score of zero with no recognized contributors means `balanced with
limited evidence`, not proven neutral electron density.

### 6.2 Role-aware interpretation

The same electronic observation has different reactivity implications:

- **aryl electrophile:** report local aryl electron demand; optionally derive an
  `electrophile_activation` assessment only when the transformation or condition
  rule defines how that demand matters;
- **N-H partner:** prioritize lone-pair class, resonance delocalization, formal
  charge, N substitution, and acidity. Derive `lone_pair_availability` separately
  from ring electron demand;
- **other partners:** add a role interpretation only through a registered,
  versioned rule.

For an unassigned family, the UI should avoid claims such as “electron-poor is
more reactive.” It should report the observation and any broadly valid caution,
while leaving mechanism-specific effects unresolved.

## 7. Recommended concise display

Rename the section from `Steric/electronic analysis` to `Reactivity profile`.
Use one line per partner and a fixed field order:

```text
Reactivity profile:
  Electrophile — Ar-Br | ring: benzene (6-membered) |
    ortho sterics: high (2/2 occupied) | ring electronics: balanced
  N partner — AromN-H | ring: <graph-derived family> |
    N environment: secondary N-H, lone pair delocalized |
    sterics: open | lone-pair availability: low
```

For a single-line UI:

```text
Electrophile Ar-Br — benzene; ortho burden high (2/2); electron demand balanced
N partner AromN-H — <ring family>; access open; N lone pair delocalized/low
```

Default display rules:

- show ring family and ring size for every `Ar` or `HetAr` anchor;
- show heteroatom positions only when present;
- show ortho occupancy only for aromatic anchors;
- show N substitution and lone-pair availability for an N reactive center;
- omit zero-valued contributor lists, method names, raw shell atom counts, and
  exact continuous scores;
- expose scores, contributors, method, confidence, and raw graph evidence in an
  expandable details view or canonical JSON.

`N partner` is preferable to `nucleophile` when the role interpretation is not
confident. The structured partner role and its confidence remain unchanged.

## 8. Signature and retrieval use

The current L0 and L4 identities include the full partner environment.
Introducing the new descriptors therefore cannot be treated as a display-only
schema edit.

Use stable, binned tokens in signatures:

```text
aromatic_family:benzene
aromatic_ring_size:6
aromatic_heteroatom:N@ortho:pyridine_like
steric_access:hindered
ortho_occupancy:2_of_2
ortho_burden:high
electronic_demand:balanced
lone_pair:pyrrole_like
lone_pair_availability:low
```

Do not place raw floating-point scores, display labels, atom indices, or source
reaction names in signature identity. During scoring, continuous values may be
compared separately with normalized distances. Partner features must be aligned
by participating edit atoms and optional roles, not by input component order.

Aromatic family is an environment feature, not a mandatory routing key. Hard
filters should still be based on graph transformations, compatible handles,
valence, and explicit chemistry rules.

## 9. Proposed implementation sequence

### Stage 1: definitions and typed observations

1. Add a versioned aromatic-family vocabulary and descriptor-bin definitions.
2. Add the typed descriptor contracts in `reactive_taxonomy`.
3. Implement ring-system perception, invariant relative heteroatom positions,
   graph steric burden, and auditable electronic contributions.
4. Emit the new profile in shadow fields without changing the existing
   `ReactionSignature`, admission, retrieval, or recommendation ranking.

### Stage 2: display migration

1. Add a single renderer owned by the descriptor contract rather than separate
   GUI, CLI, and CSV interpretations.
2. Switch the GUI and concise review output to the clean summary.
3. Keep an expanded evidence view for raw values and provenance.
4. Add snapshot tests for display ordering, omitted low-value fields, and
   unclassified or low-confidence cases.

### Stage 3: signature and scoring migration

1. Evaluate new environment features against the existing corpus in shadow
   mode.
2. Define categorical signature tokens and continuous similarity components.
3. Intentionally bump descriptor, signature-feature, and reaction-signature
   versions when activating the new environment identity.
4. Rebuild indexes and explain coverage, retrieval-level, and ranking changes.
5. Calibrate thresholds from reaction outcomes where sufficient data exist;
   otherwise label them chemistry priors.

### Stage 4: remove compatibility fields

Remove legacy `class`, `local_heavy_atoms_r2`, and `qualitative_sum` paths only
after:

- new typed descriptors serialize in all canonical records;
- GUI, CLI, CSV, condition rules, signatures, and similarity no longer consume
  the old keys;
- parity and migration tests pass;
- old index/schema versions are explicitly rejected or migrated.

## 10. Minimum regression set

The descriptor implementation should cover:

- unsubstituted, mono-ortho, and di-ortho aryl halides;
- equal occupancy with different substituent burden, such as methyl versus
  tert-butyl;
- 2-, 3-, and 4-halopyridines;
- diazines with distinct ring-N arrangements;
- five-membered pyridine-like and pyrrole-like heteroaromatics;
- fused benzene, indole, quinoline, and carbazole systems;
- aniline, diarylamine, and aromatic N-H partners with different lone-pair
  delocalization;
- formal charge and protonation changes;
- symmetric ring numbering and alternative valid SMILES serialization;
- irrelevant reactant-order changes;
- unknown-family mapped reactions;
- conflicting source labels versus graph-derived ring identity;
- deterministic descriptor IDs and signature IDs.

Tests should assert positive, negative, ambiguous, and low-confidence cases.

## 11. Decision summary

Adopt four consistent descriptor blocks for every reaction partner:

```text
aromatic system | steric accessibility | electronic demand | reactive center
```

For aromatic systems, always capture ring family, ring sizes, fusion, and
heteroatom types and positions. Use ortho occupancy plus a size-weighted burden
for sterics. Replace the ambiguous “local electronic score” with auditable
electronic-demand contributions, and treat N lone-pair availability as a
separate reactive-center property. Keep the primary display qualitative and
compact; retain numeric evidence and provenance in structured output.
