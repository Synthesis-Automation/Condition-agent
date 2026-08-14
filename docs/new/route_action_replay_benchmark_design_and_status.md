# Route-Action Replay Benchmark: Design and Status

**Status date:** 2026-08-14

## Decision

Route-action replay is the first use of the observed multistep corpus for model
development. At every observed target or intermediate, it runs the existing
validated single-step operator ladder and measures whether the recorded action
is recovered. It does not learn chemistry, bypass forward validation, or turn
route motifs into executable rules.

Only an observed step that passes the existing data-driven compiler, source
round trip, reaction-core quality, product completeness, and identity checks is
an eligible positive. Every other step remains in the artifact with its failure
stage and reason.

## Contracts

Schema `1.0` has three immutable records:

- `RouteActionEvaluation`: route provenance, split, search configuration, and
  every reaction occurrence;
- `RouteActionStepEvaluation`: the observed exact precursor, SITE1, OP1, SYN1,
  STRAT1 identities, recovery ranks, outcome, and search diagnostics; and
- `RouteActionCandidate`: a compact validated action with ranking features,
  evidence, matching flags, and visible-precedent leakage evidence.

Candidate labels are explicit:

- `observed_exact`;
- `strategy_equivalent`;
- `same_site_operator`;
- `same_site`; and
- `hard_negative`.

Exact precursor rank is over concrete returned candidates. SITE1, OP1, SYN1,
and STRAT1 ranks are over first occurrences of distinct identities, so repeated
handle realizations do not consume multiple strategy ranks.

The converter writes deterministic gzip JSONL atomically. The manifest records
the source-tree and operator-library hashes, complete search configuration,
aggregate and split metrics, eligibility reasons, search work, and output hash.
Elapsed time is manifest metadata and is excluded from record identity.

## POC audit

The strict top-25 replay over the existing random 50-route review sample used
the full-scale v3 library with 31,791 templates. It produced:

| Measure | Result |
| --- | ---: |
| Routes | 50 |
| Observed reaction occurrences | 188 |
| Eligible positive steps | 18 (9.6%) |
| Core-not-verified steps | 169 |
| Source-round-trip failures | 1 |
| Eligible steps with candidates | 18/18 |
| Returned validated candidates | 450 |
| Exact observed candidates | 12 |
| Same-site candidates with different operators | 16 |
| Hard negatives | 422 |
| Exact precursor recall at 1 | 33.3% |
| Exact precursor recall at 5 | 55.6% |
| Exact precursor recall at 10/25 | 66.7% |
| SITE1 recall at 1 | 44.4% |
| SITE1 recall at 5 | 66.7% |
| SITE1 recall at 10/25 | 72.2% |
| STRAT1 recall at 1 | 33.3% |
| STRAT1 recall at 5 | 61.1% |
| STRAT1 recall at 10/25 | 66.7% |
| Validation attempts | 1,800 |
| Visible same-patent precedent overlaps | 0 |
| Route conversion failures | 0 |

Artifact:

```text
results/core_retrosynthesis/route_action_evaluation/
  routes.poc.random50.action_replay.v1.jsonl.gz
  routes.poc.random50.action_replay.v1.jsonl.gz.manifest.json
```

The JSONL SHA-256 is
`ccf9ef10f4ef8614a33bc3c721afd247bf33f0b2c96d186c312eb03f25b926d5`.

These are diagnostic POC results, not release accuracy. The sample contains
only three eligible test-split steps, and the operator library has not yet been
rebuilt with all route patents and near-analogues excluded. The patent-overlap
flag covers retained candidate precedents only and cannot prove absence of
training leakage.

## Command

```powershell
python -m core_retrosynthesis evaluate-route-actions `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.random50.tree.v2.jsonl.gz `
  results/operator_retrosynthesis_poc/full_scale_v3/compact/operator_library_v3.json.gz `
  results/core_retrosynthesis/route_action_evaluation/routes.poc.random50.action_replay.v1.jsonl.gz `
  --workers 8 --overwrite
```

## Interpretation and next gate

The benchmark infrastructure works, candidate coverage is complete for eligible
steps, and the current ranking leaves useful headroom. The immediate blocker is
positive coverage: 90.4% of observed steps do not yet satisfy the strict generic
operator admission contract.

Before training a route-context reranker:

1. stratify the 169 core-quality failures by chemistry and evidence defect;
2. improve generic single-step support only where graph evidence is sufficient;
3. build an operator library that excludes evaluation patents and close chemical
   analogues;
4. replay the full 5,000-route corpus with the frozen library and configuration;
5. require materially broader eligibility and a sufficiently large untouched
   test set; and
6. only then train a simple pairwise action reranker over the retained candidate
   features and route context.

The system must not weaken core quality thresholds merely to increase the number
of positives.
