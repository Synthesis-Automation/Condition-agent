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

The versioned `observed_route_action_label.v1@1.0` policy requires:

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

## 50-route supervision audit

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
  routes.poc.random50.action_labels.v2.jsonl.gz
  routes.poc.random50.action_labels.v2.jsonl.gz.manifest.json
```

The JSONL SHA-256 is
`37be39eb3ca4b2b79ca3abf74fb1abd6086649dd39399eee8b14851eee178253`.

## Command

```powershell
python -m core_retrosynthesis evaluate-route-actions `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.random50.tree.v2.jsonl.gz `
  results/operator_retrosynthesis_poc/full_scale_v3/compact/operator_library_v3.json.gz `
  results/core_retrosynthesis/route_action_evaluation/routes.poc.random50.action_labels.v2.jsonl.gz `
  --labels-only --workers 8 --overwrite
```

Remove `--labels-only` to run candidate replay. Large replay jobs should be
partitioned into deterministic route shards because forward validation, not
label extraction, dominates runtime.

## Interpretation and next gate

The supervision contract is now broad enough for a route-policy POC. Before
training:

1. conduct chemist review of promoted review-grade labels, stratified by
   departing edit type;
2. add deterministic resumable replay shards and merge validation;
3. rebuild an evaluation library excluding route patents and close analogues;
4. replay the full 5,000-route corpus with frozen search configuration;
5. train a confidence-weighted or positive-unlabeled action reranker; and
6. keep deterministic candidate generation and forward validation authoritative.

The 45 observations without product-site identity can still supervise retained
OP1/SYN1 and exact-precursor tasks. A future unary/transformation-site identity
may cover reductions and deprotections, but it must be versioned and matched on
both observed and candidate reactions before it contributes to SITE1 or STRAT1
metrics.
