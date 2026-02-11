# Taxonomy Motif Backlog (10x100 Pilot)

This backlog is derived from:

- `results/reaction_coverage_discovery.10x100.json`
- `results/taxonomy_expansion_todos.10x100.json`

Scope: 10 reaction datasets, up to 100 rows per dataset (`920` reactions total).

## Applied Taxonomy Edits

The following reaction-type updates were applied in `chemtools/taxonomy/data/reaction_types.v4.0.json`:

1. Added `Alcohol_to_Alkyl_Halide`

- Captures alcohol-to-halide conversion motifs (Appel/SOCl2/PBr3/deoxyfluorination-like patterns).
- Added aliases: `Appel_halogenation`, `Appel halogenation`, `Chlorination_SOCl2_oxalyl_chloride`, `Chlorination SOCl2 oxalyl chloride`, `Deoxy_fluorination`.

1. Added `Aliphatic_Halide_Exchange`

- Captures halide-exchange patterns (Finkelstein-like) such as `RCH2-Cl -> RCH2-I`.
- Added aliases: `Aliphatic_Halide_Exchange`, `Aliphatic Halide Exchange`, `Finkelstein`.

1. Alias cleanup

- `BalzSchiemann` and `Sandmeyer_reactions` now map to `Sandmeyer`.
- Removed unrelated alias overloads from `Halogenation_aromatic` to reduce misrouting.

## Top Unresolved Clusters (Non-`none -> none`)

1. `Alkyl-Si*|RCH2-NR2|RCH2-OR -> R2CH-CO2H|RCH2-NH2|R_acidic-H`

- Count: 27
- Main source label: `Aldol (classic & Mukaiyama)`
- Reason profile: mostly `unknown_reaction_type` with frequent `mapping_warning`.
- Backlog action: likely multi-transform/noisy record handling, not immediate motif addition.

1. `none -> Alkenyl-Br`

- Count: 22
- Source labels: `Addition of Halogens to Double or Triple Bonds`, `Allylic, Benzylic and Vinylic Halogenations`
- Backlog action: improve reacted-motif extraction when only product-side halide motif is observed.

1. `R2CH-OH -> Alkenyl-Alkenyl|RCH2-Br`

- Count: 19
- Source label: `Addition of Halogens to Double or Triple Bonds`
- Backlog action: likely mixed/multi-event records; add event-gating before reaction-family assignment.

1. `Ar-C=N|Ar-NR2|Bn-OH -> Ar-Ar|Ar-AromN|Bn-Br`

- Count: 9
- Source labels mixed; likely multi-transform.
- Backlog action: route to multistep bucket before taxonomy assignment.

1. `Alkyl-H|Ar-COR|Ar-H -> Ar-C=N`

- Count: 7
- Candidate overlap low with existing families.
- Backlog action: inspect whether this is a real missing family or extraction artifact.

## Priority Motif Backlog

Top motifs currently outside reaction-slot taxonomy (weighted by unresolved counts):

1. `Alkyl-N(R)CO2R` (72)

- Often present in reductive amination/alkylation outputs.
- Action: verify whether it should be promoted into additional reaction product slots.

1. `Ar-C=N` (35)

- Frequent in unresolved clusters.
- Action: evaluate if dedicated imine/intermediate families are needed or if this is transient/noisy.

1. `Alkyl-H` (16)

- Likely broad/reactivity-context motif; often ambiguous as reaction-center signal.
- Action: keep as low-priority for direct slot usage; prefer contextual/event handling.

1. `Alkyl-CF3` (9)

- Action: inspect fluorination/trifluoromethylation families for missing product coverage.

1. `Inorganic` (37)

- Not a target organic motif for reaction-family slots.
- Action: keep in reagent/spectator handling, avoid using as reaction center slot motif.

## Recommended Next Work Cycle

1. Re-run discovery on the same 10x100 sample after the applied taxonomy edits.
2. Compare deltas in:

- unknown reaction-type rate,
- unresolved rate,
- top cluster composition.

3. If `none -> Alkenyl-Br` remains dominant, prioritize extraction-layer fix over motif expansion.
2. Only add new motifs/families when a cluster remains high-count after extraction + alias corrections.
