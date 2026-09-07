# Bounded search exploration: implementation and effectiveness review

Reviewed: 2026-09-07. This is a development experiment, not a release benchmark.

The first implementation addresses lost downstream choices and incorrect
operator-order guidance within the existing `plan_multistep_routes` path.
Partition contracts, chemistry admission, the operator ladder, and terminal
predicates retain their existing responsibilities.

## Implemented changes

1. **Dependency-aware ordering.** `LiteratureRouteOrderingGuidance` now follows
   reaction parent/child occurrence IDs. It never treats adjacent expansions on
   sibling branches as a chemical sequence. Missing or conflicting producer links
   receive no ordering evidence. Scores average log support over the actual
   dependencies and are invariant to traversal order. This is a structural
   parent/child prior, not proof of a same-site chemical enabling relationship.
   The ordering algorithm is `dependency_route_ordering.v2`; the mined catalogue
   remains `route_state_learning.v1` because its observations did not change.
2. **Bounded guidance.** Blocked-state status precedes optional guidance. Guidance
   can reorder states only within the same accumulated-cost band, currently one
   cost unit. It cannot outweigh an arbitrarily more expensive route. The band
   and widening interval are validated in `search_exploration.v1.json`.
3. **Deferred widening.** With `widening_factor > 1`, the planner requests a larger
   bounded pool from the same validated one-step expander. It initially releases
   `per_step_top_k` actions and stores continuations for the remaining batches.
   A continuation receives a turn after four expansions, or when the ordinary
   queue empties. Revisited pools are cached: releasing another batch does not
   rerun one-step chemistry. Each release consumes an expansion, stays subject
   to depth/cycle/validation/exclusion gates, and is recorded in diagnostics.
   The ordinary beam bounds ordinary states; deferred work is separately bounded
   by the action-pool and total expansion limits and is counted in frontier size.
   Pending root batches are not returned as zero-step partial routes. Partition
   reports also retain the exploration-policy version; historical reports with
   no recorded version remain unspecified.
4. **Audit fields.** Multistep results use schema `1.9` and expose exploration and
   guidance versions, widening factor, revisit count, and pending deferred states.
   Widening defaults to `1` (disabled). Ordinary unguided behavior is preserved.

This is progressive **release of a bounded validated pool**, not resumable
template enumeration. A larger initial pool can cost more validation work, and
its top choices can differ from a smaller pool after diversity ranking. Neither
baseline preservation nor improved completion is guaranteed. Exhausting this
pool is also not proof that no additional chemistry exists.

## Reproduce the comparison

```powershell
python -m core_retrosynthesis.search_effectiveness `
  results/core_retrosynthesis/synthetic_partition/dataset10_kway_review.v1.json `
  results/core_retrosynthesis/route_step_operator_library/v1/operators_validated_departures/operator_library_v3.json.gz `
  results/literature_molecule_index.sqlite `
  results/core_retrosynthesis/search_effectiveness/v1_parallel `
  --catalog results/core_retrosynthesis/route_state_learning/route_state_learning.v1.json.gz `
  --allow-untyped-literature-terminals --workers 4
```

Frozen defaults for this experiment:

| Control | Value |
| --- | ---: |
| Maximum route depth | 6 |
| Expansion budget per target | 30 |
| Beam width | 8 |
| Returned routes | 3 |
| Initial actions per expansion | 3 |
| Wide pool | 9 |
| Applied templates per ladder level | 100 |
| Validation attempts per ladder level | 25 |
| Molecular-weight terminal threshold | 150 Da |

The same targets, library, stock index, depth, expansion budget, and beam are used
for four modes:

- `baseline`: ordinary top-three search;
- `fixed_wide`: ordinary top-nine search, controlling for additional breadth;
- `deferred_wide`: top-three batches from a pool of nine;
- `deferred_ordered`: deferred widening plus dependency-aware ordering, with no
  latent-state reservation selector, so the ordering effect is isolated.

Depth six removes the known depth ceiling for all ten reference routes. This
baseline is a new controlled run and must not be compared as if it used the older
depth-three/15-expansion configuration.

`comparison.json` records input SHA-256 hashes, configurations, metrics, actual
validation counts, runtimes, and per-target positive/negative deltas. Each mode
is checkpointed to JSON and HTML before the next mode starts.
`comparison.html` links to paired graphical routes with review notes and export.
`--workers 4` runs the modes in independent processes with separate index
connections used only for lookups. The default `--workers 1` runs serially. Both retain the same
per-target deterministic search controls; the worker count is recorded.

The expansion ceilings are equal; validation work and runtime are **measured,
not equalized**. Runtime is descriptive: serial modes can benefit from warm
chemistry caches; concurrent modes compete for machine resources. A gain over narrow search alone
cannot establish scheduling efficiency: compare against `fixed_wide` as well.

## Interpret the results

The panel contains nine train routes and one previously inspected test route.
It is not untouched and cannot establish generalization. No parameters should be
tuned against the remaining untouched evaluation set.

- Report heuristic-complete and supplier-complete targets separately. A literature
  occurrence or a molecule below 150 Da is not commercial availability evidence.
- Report observed root recovery and maximum reference action matches separately
  from route completion. Matching the action multiset does not independently
  verify reference-route topology or experimental usefulness.
- Inspect per-case losses as well as aggregate gains. A new partial route is not
  automatically a useful alternative.
- Compare validation counts and expansions alongside recovery. Widening releases
  consume expansion budget even when the pool was already validated.
- Graph and signature validation remain prerequisites. Condition feasibility,
  experimental selectivity, and protection necessity are not established by this
  ablation.

## Measured outcome and decision

The completed run is in
[`v1_parallel/comparison.html`](../../results/core_retrosynthesis/search_effectiveness/v1_parallel/comparison.html),
with machine-readable controls, hashes, and per-case deltas in
[`comparison.json`](../../results/core_retrosynthesis/search_effectiveness/v1_parallel/comparison.json).

| Mode | Heuristic-complete targets | Known actions retained | Known roots retained | Validation attempts |
| --- | ---: | ---: | ---: | ---: |
| Baseline, top three | 6/10 | 13/40 | 8/10 | 4,171 |
| Immediate top nine | 5/10 | 9/40 | 6/10 | 6,899 |
| Deferred widening | 4/10 | 3/40 | 3/10 | 6,075 |
| Deferred widening + corrected ordering | 4/10 | 3/40 | 3/10 | 6,250 |

All modes have **zero supplier-complete targets** with the supplied literature
index. Every wider mode regressed on both completion and reference-action
retention while spending more validation work. Corrected ordering added no
recovery gain over deferred widening on this panel.

Immediate top nine loses completion for case 3 and loses known actions for cases
1 and 3. Both deferred modes lose completion for cases 3 and 8 and lose known
actions for cases 1, 3, 4, 8, and 10. There are no completion or known-action gains
to offset these losses on this panel.

Among displayed routes for partial-only targets, median reaction count is five
for the baseline, four for immediate top nine, and one for both deferred modes.
This exposes a continuation/selection problem: the broader searches consume the
budget while retaining much shallower partial prefixes. It does not prove that
every short alternative is chemically inferior. A larger initial reservoir can
also change the one-step diversity ordering; this run does not isolate every
contribution of generation, scheduling, and final partial-route selection.

**Decision: do not promote widening into the default planner.** Keep the
dependency/ambiguity correctness fix and bounded guidance, and retain widening
only as an explicit experimental control. Increasing breadth alone is not the
next justified improvement. The next experiment should measure continuation
value and partial-route progress, with a trace distinguishing generated,
admitted, explored, and displayed actions. Any subsequent policy should preserve
the existing useful candidate prefix and be compared against the frozen baseline.
These are follow-up hypotheses, not improvements demonstrated by this run.

The serial baseline in `v1/baseline.json` and the independently executed baseline
in `v1_parallel/baseline.json` are byte-identical. All evaluated portfolios have
zero zero-step partial routes. Runtime is reported but is not used to claim an
efficiency win, because processes share machine resources.

Verification: **1,268 tests passed** in the final full `pytest -q` invocation.
This includes the dependency, widening, low-budget reporting, policy-provenance,
CLI, process-isolation, and report tests. A real-library process smoke run also
completed. HTML checks verified molecular SVGs, paired metrics, review controls,
and export markup. Browser visual/interaction QA was unavailable because the
browser runtime reported no connected browser.

## Use on a new target

```powershell
python -m core_retrosynthesis plan-routes LIBRARY STOCK_INDEX TARGET `
  --max-depth 6 --per-step-top-k 3 --widening-factor 3
```

The existing `evaluate-partition-review-routes` command also accepts
`--widening-factor`. Optional `--route-state-catalog` plus
`--route-state-ordering` on `plan-routes` includes the existing state reservation
as well as corrected ordering; the comparison runner isolates ordering alone.

## Verification and promotion gate

Fast behavior tests:

```powershell
pytest -q tests/core_retrosynthesis_tests/test_search_exploration.py tests/core_retrosynthesis_tests/test_search_effectiveness.py
```

They cover late-action recovery, reuse without revalidation, deterministic
results, expansion exhaustion, invalid later actions, sibling traversal
invariance, conflicting/missing dependency links, guidance bounds, CLI forwarding,
and report pairing/stock-evidence distinctions. The full `pytest -q` suite is the
handoff gate.

Before enabling widening by default, freeze a separate development/validation
panel, compare at measured compute budgets with blinded chemistry review, and
then run the untouched set once. Continue with explicit precursor-state
requirements only if this experiment supports the added search complexity.
