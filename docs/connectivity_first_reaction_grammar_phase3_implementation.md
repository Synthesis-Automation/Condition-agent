# Connectivity-First Reaction Grammar: Phase 3 Implementation

## Status

Implemented on 2026-07-28.

Phase 3 introduces immutable normalized reactive-site interfaces over the
existing chemistry-specific detectors. Detector results remain authoritative
observations and retain their site IDs, site types, signatures, chemist labels,
availability, and context. The normalized interfaces are internal adapter
views and do not change molecule or reaction serialization.

The generic connectivity-rewrite executor remains a shadow implementation.
Current operators still determine production reconstruction.

## 1. Delivered contracts

`reaction_site_interfaces.py` defines:

- `ConnectivityEndpoint`, representing either a real atom or one
  schema-level virtual hydrogen with carrier provenance;
- `ReactiveLinkSite`, representing a consumable explicit A-B or implicit A-H
  link;
- `BondCapacitySite`, representing bounded localized bond-order decrement and
  increment capacity;
- `ConnectionEndpointSite`, representing a curated atom endpoint plus any
  required link release, hydrogen loss, formal-charge change, or bond-capacity
  consumption; and
- `NormalizedSiteInterfaces`, retaining all normalized views derived from one
  detector site.

All contracts are frozen dataclasses with schema version `1.0`, typed endpoint
kinds, deterministic component-qualified IDs, explicit invariants, and
JSON-serializable `to_dict` views.

No interface enumerates arbitrary atoms. A normalized connection endpoint
exists only when a curated detector has already supplied the site.

## 2. Detector adapters

The adapters currently normalize:

| Existing detector | Normalized view |
|---|---|
| `leaving_group` | explicit reactive link plus two bounded connection endpoints |
| `transfer_group` | explicit reactive link plus two bounded connection endpoints |
| `pronucleophile_XH` | atom/virtual-H link plus an endpoint requiring H loss |
| `aromatic_CH` | atom/virtual-H link plus an endpoint requiring H loss |
| `addition_donor` explicit A-B | explicit link with endpoint symmetry |
| `addition_donor` implicit A-H | atom/virtual-H link with carrier provenance |
| `unsaturated_bond` | localized bond capacity plus two connection endpoints requiring one consumed bond unit |

Examples normalized through the same contracts include C-Br, C-O-sulfonate,
C-B, C-Si, N-H, O-H, S-H, B-H, Si-H, Br-Br, Cl-Br, B-B, Si-B, C=C, C#C,
and C#N.

Chemistry-specific detector labels are retained in `source_site_type`,
`source_signature`, `source_chemist_label`, and `annotation_tokens`. They are
annotations, not alternate routing keys.

## 3. Rewrite integration

Connectivity rewrite schema and instruction-set versions are now `1.1`.
Definitions declare the required site-interface schema version.

The bounded selector vocabulary now supports:

```text
role.reactive_link.endpoint_a
role.reactive_link.endpoint_b
role.reactive_link.carrier
role.bond_capacity.endpoint_a
role.bond_capacity.endpoint_b
role.connection_endpoint.atom
```

Interface predicates may inspect only a fixed whitelist of typed fields.
Unknown interface names, fields, or selector members are rejected during
static compilation.

Suzuki coupling, C-N release-and-connect, and alkene pair addition now execute
through normalized selectors. The beta-elimination pilot retains its existing
atom-role selectors until the remaining site/operator migration phase.

`xh_addition_to_alkene` now shares the same normalized
split-and-distribute rewrite as explicit A-B, Si-H, and B-H addition. A virtual
hydrogen is compiled to carrier H loss and acceptor-endpoint H gain; it is
never treated as a heavy atom or product seed.

Compatibility bond-change labels are still emitted from each normalized
endpoint's source atom role, preserving exact legacy edit parity.

## 4. Compatibility boundary

Phase 3 does not:

- change existing detector output or `ReactiveSite` serialization;
- attach normalized interfaces to public molecule or reaction payloads;
- change `reaction_grammars.v1.json` or reaction-signature identity;
- change candidate generation, selection, admission, retrieval, or scoring;
- make the generic executor authoritative;
- normalize every current site family;
- infer generic reactivity for otherwise inert atoms; or
- support aromatic bond arithmetic, coordination, radicals, or isotopes.

`normalize_compound_sites` and `normalize_reaction_assignment` are explicit
inspection APIs. Calling them does not mutate their source analyses.

## 5. Validation matrix

The Phase 3 regressions cover:

- halide and sulfonate leaving links;
- carbon-boron and carbon-silicon transfer links;
- distinct transfer and Si-H addition modes on one molecule;
- N-H, O-H, S-H, B-H, and Si-H virtual endpoints;
- symmetric Br-Br and asymmetric Cl-Br links;
- alkene, alkyne, and nitrile bond capacity;
- required link-release, H-loss, and bond-capacity metadata;
- no promotion of inert alkane atoms;
- unique component-qualified interface IDs;
- component-order-invariant interface chemistry;
- unchanged public molecule serialization;
- normalized X-H addition product/edit parity; and
- compiler rejection of unknown normalized selectors.

The focused interface and rewrite gates pass 44 tests, and the complete
repository suite passes 537 tests.

Focused gates:

```powershell
pytest -q tests/reactive_taxonomy/test_reaction_site_interfaces.py
pytest -q tests/reactive_taxonomy/test_connectivity_rewrite_executor.py
```

Repository-wide gate:

```powershell
pytest -q
```

## 6. Handoff to Phase 4

Phase 4 may migrate the remaining operator definitions to normalized
connectivity programs:

1. adapt eliminable pairs and simple order-change sites;
2. add alkyne pair addition variants;
3. cover C-O and C-S release-and-connect;
4. cover activated acyl/sulfonyl replacement;
5. dual-run every migrated grammar across the regression corpus; and
6. remove a legacy operator only after exact product, edit, warning,
   ambiguity, stereochemistry, and ordering parity.

The generic executor must remain non-authoritative until those gates pass.
