# Observed Route-Action Labels and Replay: Design and Status

**Status date:** 2026-08-14

## Decision

Route-learning evidence and executable-operator admission are separate
contracts. A reported reaction can reliably identify the chosen product site,
retained transformation, synthon partition, and exact precursor realization
even when a leaving group disappears and prevents every reactant-side bond edit
from being checked on the product graph.

`ObservedRouteActionLabel` therefore exposes independent supervision facets.
It does not relax generic operator compilation: executable templates still need
the existing strict core-quality and source-round-trip checks.

## Evidence facets

The versioned `observed_route_action_label.v1@1.1` policy requires:

- verified product completeness;
- matching route and reaction product identity;
- complete active-atom mapping;
- a changed primary center and resolved retained continuity;
- no blocking core reason; and
- for review-grade cores, only `not_all_edits_graph_checked` as a review reason.

The last case admits route supervision only when unresolved edits concern atoms
absent from the product. OP1 remains the existing handle-independent signature
of product-retained edits.

Every label reports these booleans independently:

- `product_site_verified`;
- `retained_edits_verified`;
- `synthon_partition_verified`;
- `exact_precursors_verified`;
- `strategy_verified`;
- `realization_verified`;
- `operator_roundtrip_verified`; and
- `search_eligible`.

It also preserves core quality, completeness, retained/departing edit counts,
departing edit descriptors, operator rejection stage/reason, and all limitations.
Departing bond descriptors order their endpoints by atomic number, preventing
orientation aliases such as `C-O` versus `O-C` from splitting audit strata.

Schema `2.0` nests this label inside every `RouteActionStepEvaluation`.
Candidate replay remains optional. `--labels-only` extracts supervision without
running the expensive single-step search.

## Candidate supervision

When replay is enabled, every returned candidate remains hard graph- and
signature-validated. Its relationship to the observed route is one of:

- `observed_exact`;
- `strategy_equivalent`;
- `same_site_operator`;
- `same_site`; or
- `unchosen_alternative`.

An unchosen alternative is not asserted to be chemically wrong. It is weak
preference evidence because the patent records one concrete route rather than a
comparison against every feasible route.

Exact precursor rank is over concrete candidates. SITE1, OP1, SYN1, and STRAT1
ranks are over first occurrences of distinct identities. Each recall metric has
its own evidence-appropriate denominator.

## Supervision audits

### Random 50-route regression

The complete labels-only audit produced:

| Measure | Result |
| --- | ---: |
| Routes | 50 |
| Observed reaction occurrences | 188 |
| Retained-edit labels | 184 (97.9%) |
| Synthon-partition labels | 184 (97.9%) |
| Exact-precursor labels | 184 (97.9%) |
| Product-site labels | 139 (73.9%) |
| STRAT1 labels | 139 (73.9%) |
| Realization labels | 139 (73.9%) |
| Executable operator round trips | 18 (9.6%) |
| Blocked observations | 4 |
| Product-site identity unavailable | 45 |
| Route conversion failures | 0 |

The original 9.6% result was valid for executable operators but was the wrong
denominator for route learning. The new contract raises STRAT1 supervision to
73.9% and retained-transformation supervision to 97.9% without admitting one
additional operator into the library.

Artifact:

```text
results/core_retrosynthesis/route_action_evaluation/
  routes.poc.random50.action_labels.v3.jsonl.gz
  routes.poc.random50.action_labels.v3.jsonl.gz.manifest.json
```

The JSONL SHA-256 is
`4a09e5e9201fc3e7ecbc0ef65e1fee234ae0dcfeaa136b24ca5b270ec37bd958`.

Converter `3.1` avoids repeating graph analysis for review-grade steps and
avoids loading the search library in labels-only workers. Its shard manifest
also binds the label-policy definition, so `--resume` rejects stale labels.

### Complete 5,000-route POC

The full deterministic audit completed with no rejected routes:

| Measure | Result |
| --- | ---: |
| Routes | 5,000 |
| Observed reaction occurrences | 18,647 |
| Retained-edit labels | 18,280 (98.0%) |
| Synthon-partition labels | 18,078 (96.9%) |
| Exact-precursor labels | 18,307 (98.2%) |
| Product-site labels | 13,170 (70.6%) |
| STRAT1 labels | 13,099 (70.2%) |
| Realization labels | 13,099 (70.2%) |
| Executable operator round trips | 1,201 (6.4%) |
| Blocked observations | 326 |
| Route conversion failures | 0 |

The split-level fractions are stable: STRAT1 coverage is 70.9% on test, 70.2%
on train, and 69.9% on validation. This makes a split-specific extraction bug
unlikely. The merged artifact is:

```text
results/core_retrosynthesis/route_action_evaluation/
  routes.poc.full5000.action_labels.v3.jsonl.gz
  routes.poc.full5000.action_labels.v3.jsonl.gz.manifest.json
```

Its JSONL SHA-256 is
`e266b8262e2efe86ad3cc8b5bd903ddd80de29ba36226bdbbe3194bd2550d286`.

The most common product-absent bond evidence is `C-O` single (7,403), `C-Cl`
single (2,434), `C-N` single (2,265), `C-Br` single (1,809), and `C-O` double
(1,505). These groups should dominate the first chemist audit because they
cover the largest promoted-label populations. The manifest contains the full
descriptor and optional annotation distributions rather than only this summary.

## Command

```powershell
python -m core_retrosynthesis evaluate-route-actions `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.random50.tree.v2.jsonl.gz `
  results/operator_retrosynthesis_poc/full_scale_v3/compact/operator_library_v3.json.gz `
  results/core_retrosynthesis/route_action_evaluation/routes.poc.random50.action_labels.v3.jsonl.gz `
  --labels-only --workers 8 --overwrite
```

Large jobs use stable SHA-256 route-ID partitions. A shard is safely reused only
when `--resume` validates its source, operator library, search configuration,
selection, converter contract, output checksum, and shard identity. The merge
requires every index exactly once, verifies every checksum and route assignment,
rejects duplicate route IDs, and writes routes in deterministic identity order.

```powershell
python -m core_retrosynthesis evaluate-route-actions SOURCE LIBRARY PART `
  --labels-only --workers 8 --shard-count 16 --shard-index 0 --resume

python -m core_retrosynthesis merge-route-action-shards MERGED PARTS...
```

Remove `--labels-only` to run candidate replay with the same shard/resume
contract. Forward validation, rather than label extraction, dominates replay
runtime.

## Chemist review

There are 17,079 promoted steps: retained-edit evidence is verified, but the
strict executable operator is not admitted. The self-contained review samples
120 examples round-robin across core quality, optional named annotation,
departing-edit descriptor, and STRAT1-versus-retained-only status:

```text
results/core_retrosynthesis/route_reviews/
  higher_level_routes_promoted_action_labels_v3.html
```

Every card shows the observed reaction, evidence facets, OP1/SITE1/SYN1/STRAT1,
departing edits, strict rejection, and limitations. Accept/reject/uncertain
status and notes persist locally and export as JSON. Static verification found
120 cards, 120 reaction SVGs, 120 status/note controls, and no render errors.
Interactive browser verification was unavailable in the generation session.

## Bounded replay smoke test

A five-route, 17-step smoke test froze `top_k=5`, 100 applied templates, 20
validation attempts, and lazy validation. All 17 searchable steps returned five
validated candidates in 61 seconds. Exact-precursor recall was 23.5% at five;
SITE1 recall was 61.5%; STRAT1 recall was 38.5%. Four candidates matched the
observed exact precursor and three additional candidates were
strategy-equivalent.

```text
results/core_retrosynthesis/route_action_evaluation/
  routes.poc.random5.action_replay.v3.jsonl.gz
  routes.poc.random5.action_replay.v3.jsonl.gz.manifest.json
```

Its JSONL SHA-256 is
`d03a5f7978ed8be99d8bc3cea9fc9c470b3ac8b49080589aef91e30f54c23bc7`.
This sample is a pipeline/runtime check only. It is too small for accuracy
claims, and the compact library has not been rebuilt to exclude evaluation
patents and close chemistry analogues.

## Interpretation and next gate

The first route-policy POC is implemented. A bounded 50-route replay yielded 81
listwise choice sets and exercised deterministic training, serialization, CLI
loading, and multistep ranking. Only five examples belonged to validation, so
the versioned safety gate selected a zero residual and the artifact has zero
planner influence. See `docs/new/multistep_route_action_policy_poc.md`.

The next gate is to:

1. complete chemist decisions in the generated promoted-label review;
2. rebuild an evaluation library excluding route patents and close analogues;
3. replay larger leakage-controlled shards with a frozen search budget;
4. require a sufficiently populated held-out validation split before activation;
5. evaluate next-action recovery and whole-route search outcomes on untouched
   test routes; and
6. keep deterministic candidate generation and forward validation authoritative.

The 45 observations without product-site identity can still supervise retained
OP1/SYN1 and exact-precursor tasks. A future unary/transformation-site identity
may cover reductions and deprotections, but it must be versioned and matched on
both observed and candidate reactions before it contributes to SITE1 or STRAT1
metrics.
