# Connectivity-First Reaction Grammar V2

## Status

Implemented on 2026-07-28.

V2 is a clean chemistry-contract cutover. It does not retain the v1 reaction
operator path, shadow execution, parity routing, or authority flags. The
molecular graph and the emitted connectivity rewrite are the source of truth.

V2 has since been superseded for environment identity by the context-aware
profile activation. The connectivity model remains current, while signatures
now use schema `3.0`, the `RS3` namespace, `signature_features.v3.json`, and
`taxonomy_manifest.v3.json`.

## Core model

The base reaction description is:

```text
molecular observations
  -> canonical connectivity sites
  -> grammar-compatible role assignment
  -> bounded connectivity rewrite
  -> normalized bond/H/charge edits
  -> edit-derived structural archetype
  -> optional family annotation
  -> versioned reaction signature
```

Unchanged atoms and bonds are context. They may affect compatibility,
conditions, ranking, and explanations, but they do not define the reaction.

## Definition bundle

The identity-bearing v2 definitions are listed by
`taxonomy_manifest.v2.json`:

- `site_patterns.v2.json`: molecular annotation patterns and atom roles;
- `context_facets.v2.json`: typed context facets, semantic IDs, and display
  tokens;
- `site_interfaces.v2.json`: declarative adapters from annotations to
  connectivity sites;
- `reaction_grammars.v2.json`: compatible site-role constraints only;
- `connectivity_rewrites.v2.json`: bounded graph rewrite programs;
- `descriptor_rules.v1.json`: local environment construction;
- `signature_features.v2.json`: L0-L4 signature projections.

Reaction signature identity contains the schema version and a content hash for
every identity-bearing definition. A definition edit therefore cannot silently
reuse an old signature identity.

## Canonical molecular connectivity sites

`featurize_molecule()` now emits `connectivity_sites` directly. Detector
annotations remain useful provenance, but reaction execution consumes these
canonical interfaces:

- `ReactiveLinkSite`: an explicit A-B bond or atom-H link that can be
  released;
- `BondCapacitySite`: a localized bond whose order can increase or decrease;
- `ConnectionEndpointSite`: an atom that can form a new connection subject to
  declared release, H, charge, or capacity requirements.

Each endpoint retains component, atom index, element, charge, aromaticity,
hybridization, valence, hydrogen count, ring state, local environment, and
context tokens. Reaction parsing rebases molecule-local IDs into reaction
component scope and rejects assignments that contradict the molecular
observation.

## General A-B addition

A-B addition to a localized multiple bond is represented uniformly:

```text
A-B              broken, or A-H consumed at schema level
E1=E2            order decreased by one
E1-A             formed
E2-B             formed
```

The same rewrite accepts:

- H-H and X-H;
- N-H, O-H, and S-H;
- Si-H and B-H;
- X-X, including Br-Br;
- interhalogens such as I-Cl;
- Si-B and other explicitly detected A-B donor bonds.

The grammar does not encode a mechanism, catalyst, regioselectivity rule, or
named reaction. Symmetric outcomes collapse deterministically. Distinguishable
endpoint assignments are enumerated and only product evidence may select one.

## Substitution, coupling, and elimination

Substitution and most two-partner couplings share release-and-connect:

```text
C-LG     broken
P-X      broken or P-H/charge consumed when required
C-P      formed
```

The partner may be carbon, nitrogen, oxygen, sulfur, or another registered
endpoint. “Coupling” is an optional interpretation of the same edit topology,
not a separate base executor.

Elimination is the inverse connectivity shape:

```text
C-LG     broken
adjacent C-H consumed
C-C      order increased
```

Simple reduction and oxidation are localized bond-capacity changes plus
explicit schema-level hydrogen changes.

## Derived reaction archetypes

V2 removes `edit_archetype` from reaction grammar definitions. The runtime
derives it from emitted or observed changes:

- bond-order decrease plus new connection/H gain: `addition`;
- bond-order increase plus bond/H loss: `elimination`;
- bond break plus bond formation: `substitution`;
- isolated multiplicity change: `bond_order_change`;
- mixed unmatched topology: `composite`;
- insufficient evidence: `unresolved`.

These are structural summaries, never mechanism assignments.

## Removed v1 execution contracts

The cutover removes:

- `reaction_operators.py`;
- connectivity rewrite parity and shadow comparison;
- per-grammar authority flags;
- construction-absence fallback to legacy operators;
- grammar `operator` objects;
- grammar-declared `edit_archetype`;
- `OperatorOutcome` and `operator_outcome_id`.

The public rewrite result is `RewriteOutcome`, and candidates store
`rewrite_outcome_id`.

## Signature and converter versions

- historical V2 reaction signature schema: `2.0`;
- historical identifier namespace: `RS2`;
- current reaction signature schema and namespace: `3.0` / `RS3`;
- current recommendation record schema: `3.1`;
- current generic converter definition: `generic_conversion.v2.1`.

Old artifacts are intentionally rejected by the v2 index validation rather
than silently mixed with new chemistry identities. They should be reconverted
from source data.

## Explicit current limits

V2 does not claim that every possible reaction is expressible. New chemistry
should first be reduced to definite bond, H, charge, and localized bond-order
changes. Additional bounded instructions are appropriate only when the
existing instruction set cannot represent those edits without inventing atom
correspondence or omitted byproducts.

High-value future edit shapes include cleavage, insertion/extrusion,
exchange/metathesis, migration, annulation/cycloaddition, and condensation
with explicit product projection. They should extend the same connectivity
site and rewrite contracts, not add named-reaction executors.

## Validation

`validate_taxonomy()` verifies:

- the manifest and all required definition files;
- context facet identities;
- site-interface adapters;
- SMARTS and atom-role integrity;
- one and only one registered rewrite for every grammar;
- absence of legacy grammar execution fields;
- signature schema agreement.

Regression coverage includes direct v2 execution across all 32 grammars,
partner-order invariance, multi-event reconstruction, invalid provenance,
mapped evidence, deterministic IDs, and explicit Br-Br, I-Cl, Si-H, and Si-B
observations.
