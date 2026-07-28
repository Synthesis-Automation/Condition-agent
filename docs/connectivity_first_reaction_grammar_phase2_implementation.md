# Connectivity-First Reaction Grammar: Phase 2 Implementation

## Status

Implemented on 2026-07-28.

Phase 2 adds a bounded, separately versioned connectivity-rewrite compiler and
executor. It deliberately does not replace the current reaction operators.
The legacy operator result remains authoritative while the new executor is
available through an exact shadow-parity comparison.

## 1. Delivered contracts

The implementation adds:

- `connectivity_rewrites.v1.json`, with its own schema and instruction-set
  versions;
- `compile_connectivity_rewrite_definitions`, which rejects unknown
  instructions, arbitrary callables, unsupported bond domains, duplicate bond
  edits, invalid selectors, and projection without a declared broken
  attachment;
- `apply_connectivity_rewrite`, which validates before-states, applies negative
  then positive bond deltas, applies schema-level H and observed charge
  changes, performs explicitly authorized projection, validates product seeds
  and valence, preserves supported existing alkene stereo, and returns
  canonical symmetry-collapsed outcomes; and
- `compare_connectivity_rewrite`, which compares products, normalized
  compatibility edits, warnings, outcome IDs, and ordering without changing
  production selection.

The executable instruction vocabulary is limited to:

```text
change_localized_bond_state
change_schema_hydrogen_count
change_observed_formal_charge
declare_product_seed
declare_projection_discardable_attachment
enumerate_endpoint_permutation
```

Definitions cannot name or import executable Python.

## 2. Pilot rewrites

Five existing grammar IDs are compiled through four shared rewrite programs:

| Rewrite | Template | Grammar coverage |
|---|---|---|
| `suzuki_release_and_connect` | `release_and_connect` | `boron_transfer_coupling` |
| `c_n_release_and_connect` | `release_and_connect` | `sp2_c_n_substitution`, `sp3_c_n_substitution` |
| `alkene_split_and_distribute` | `split_and_distribute` | `addend_pair_addition_to_alkene` |
| `beta_depart_and_unsaturate` | `depart_and_unsaturate` | `beta_halo_elimination` |

The addition rewrite uses one heavy-atom delta program for explicit A-B links
such as Br-Br and one virtual-H lowering of the same logic for links such as
Si-H. Endpoint permutations are definition data, and symmetry-equivalent
products are collapsed deterministically.

The release-and-connect rewrite supports both intermolecular coupling and
intramolecular ring closure. Product projection is allowed only for the
discarded side of a bond explicitly changed to `NONE`.

## 3. Compatibility boundary

This phase does not:

- alter `reaction_grammars.v1.json` or its version;
- add connectivity-rewrite versions to `ReactionSignature` identity;
- change `ReactionAnalysis` serialization;
- change candidate generation, selection, admission, retrieval, or scoring;
- switch production reconstruction from legacy operators;
- infer omitted byproducts; or
- approximate aromatic, dative, coordination, radical, or isotope changes.

Shared low-level bond, H, and alkene-stereo editing helpers were extracted so
the old and new executors use the same deterministic RDKit primitives without
the generic executor importing legacy operator implementations.

## 4. Parity gates

The regression matrix includes:

- Suzuki C-C coupling and alkenyl stereo retention;
- aryl and alkyl C-N formation;
- intramolecular C-N ring closure;
- Br2 and R3SiH alkene addition;
- asymmetric-alkene endpoint enumeration and stable outcome ordering;
- beta-halo elimination;
- invalid declared before-state rejection;
- unsupported grammar rejection; and
- compiler rejection of executable or unsupported bond-domain instructions.

The focused Phase 2 gate passes 19 tests, and the complete repository suite
passes 512 tests.

Run the focused gate with:

```powershell
pytest -q tests/reactive_taxonomy/test_connectivity_rewrite_executor.py
```

The repository-wide gate remains:

```powershell
pytest -q
```

## 5. Handoff to Phase 3

Phase 3 should add normalized site interfaces:

- `ReactiveLinkSite`;
- `BondCapacitySite`; and
- `ConnectionEndpointSite`.

Existing leaving-group, transfer-group, X-H, addition-donor, and
unsaturated-bond detectors should adapt to those interfaces while retaining
their current chemical annotations and site IDs. Phase 3 should not yet make
the new executor authoritative. That switch belongs after corpus-level parity
and the remaining Phase 4 operator migrations.
