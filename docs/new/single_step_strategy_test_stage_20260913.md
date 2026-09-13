# Single-step strategy search: first testable stage

## Scope

The Workbench now requests distinct verified strategies, with up to three
concrete precursor realizations per strategy. This is a search and application
integration change. Reaction observations, operator chemistry, STRAT1 identity,
shared-core definitions, source admission, and condition matching are unchanged.
Existing operator libraries can be used without rebuilding the dataset.

The canonical grouped entry point is
`core_retrosynthesis.disconnect_strategies_detailed`. The existing tuple-returning
`disconnect_strategies` delegates to it. The flat operator search remains in use
for existing planner/evaluation contracts; it is no longer the Workbench's
single-step result path.

## Implemented behavior

- Applicable templates and generated proposals are interleaved across operator
  IDs before budgets are applied. Similar variants of one operator cannot occupy
  the entire prefix when other applicable operators are available.
- Identical mapped proposals for one operator share a validation attempt.
  Distinct supplied correspondences remain separate. Final deduplication in
  this path preserves validated strategy identity as well as precursor identity.
- Each retained realization still passes forward analysis and operator-signature
  agreement. Scheduling buckets do not grant chemical admission.
- Search proceeds L2 to L1 to L0 until the requested number of distinct verified
  strategies is found, or the permitted tiers are exhausted. Disabling L0 is
  honored. Limits are **per attempted tier**, not global request limits.
- Existing structural, compatibility, hierarchical and optional precursor-realism
  ranking are reused within each tier. Group representatives retain the first
  ranked realization; correlated source support is not summed.
- Stage counters report applicable templates, generated proposals, validation
  attempts, accepted candidates, duplicate proposals, validation failures and
  budget exclusions. A budget-limited shortfall is not labeled dataset absence.
- Proposals with a verified signature but incomplete strategy identity are
  preserved as `unresolved_candidates`. They do not enter ranked strategies.
  The Workbench exposes their reactions in a separate review section.
- Each strategy has a precursor-choice selector. Conditions load for all
  representatives first, followed by retained alternatives. Each request uses
  that realization's own mapped reaction and precedent IDs. Selection survives
  progressive condition updates.

## API and application impact

`POST /api/v1/retrosynthesis` returns schema **2.0**. `top_k` now counts strategies.
The former flat `candidates` / `candidate_count` fields are replaced by
`strategies`, `strategy_count`, `requested_strategy_count`, and
`returned_realization_count`. Each strategy retains the owning package's
representative/alternate contract, enriched with application references,
condition-loading state, and optional forward audits. `search_diagnostics` and
`unresolved_candidates` are explicit. Custom clients must migrate to this shape.

The condition lookup endpoint is unchanged. No second Workbench route or
parallel legacy response is introduced. Browser JSON export contains the grouped
result. The existing additional forward audit remains advisory and distinct from
the mandatory signature validation.

## Frozen development comparison

Before changing search, four target queries and the Compact library SHA-256 were
frozen in `results/single_step_strategy_stage/baseline.json`, together with the
original search/grouping code hashes and original candidates. Both executions
use top five, 100 templates and 30 validations per attempted tier.

| Target | Original complete strategy IDs in five flat hits | New strategies | New retained precursor choices |
| --- | ---: | ---: | ---: |
| Indazole nitrile, `N#CC1=NNC2=CC=CC=C12` | 2 | 5 | 6 |
| Ethylamine, `CCN` | 2 | 5 | 8 |
| Ethanol, `CCO` | 2 | 5 | 8 |
| Biphenyl, `c1ccc(-c2ccccc2)cc1` | 1 | 5 | 9 |

The frozen baseline's raw `distinct_strategies` count includes the empty string
when incomplete IDs occur. The table deliberately counts only nonempty IDs;
the frozen source artifact is preserved. Original incomplete-hit counts are
2, 3, 2 and 0 respectively.

These queries are a consumed **development smoke panel**, not a held-out accuracy
benchmark or independent chemist review. Returning more identities does not
establish more useful syntheses. All four new searches reached five strategies
at L2 and reported budget exclusions. Fallback behavior is exercised separately
by deterministic regression tests.

The new searches preserve 6, 13, 15 and 7 incomplete-identity proposals for review
respectively. This exposes an existing identity limitation rather than repairing
it by inventing synthons. New query times were approximately 7.3, 4.9, 5.4 and
12.3 seconds versus 5.6, 4.0, 5.9 and 11.9 seconds originally. These are single
local runs, not a controlled latency benchmark; no speed improvement is claimed.
Detailed results and reproducible local scripts are in the same result folder.

## Try it

1. Restart the application: `python -m app.web_api --workbench`.
2. Select **Single-step retrosynthesis**, **Compact**, and **Top strategies: 5**.
3. Enter one of the targets above and click **Plan one step**.
4. Select a strategy, switch its **Precursor choice**, and inspect its own
   conditions, forward audit and supporting precedents.
5. Open **Search coverage** and, when present, **Unresolved strategy identities**.

No dataset rebuild is needed. Frontend assets must be rebuilt after pulling
source changes (`--build` on the launcher or `npm run build` in the web folder).

In the workspace used for this validation, the configured **Compact** operator
library is present (3,710 operators / 31,980 templates), but the **Full** operator
library is absent. A direct Full request correctly returns HTTP 503 with
`retrosynthesis operator library is unavailable for full mode`. The Full
condition index is a separate artifact and does not supply executable retro
operators. Use Compact for this test stage; Full needs its operator artifact
built or configured before it can be tested.

## Verification

- Focused strategy tests cover distinct-result fallback, no-L0 behavior,
  operator scheduling, identity ambiguity/conflicts, mapped reconstruction,
  deterministic results, and preservation of invalid/unresolved/conflicting
  validation gates.
- Browser integration passed on Compact: five distinct strategies, selection of
  an alternate precursor, its own condition request, stable selection throughout
  progressive updates, search counters, and zero JavaScript errors.
- The production frontend build passes. The browser screenshot is
  `results/single_step_strategy_stage/workbench.png`.
- Complete Python suite: **1,473 passed in 375.98 seconds**. Ruff Python
  error checks and `git diff --check` also pass.

## Remaining foundation work

This stage balances work across **operators**, not full provisional target-site
and synthon groups. Different sites within one operator can still compete for
budget. It also stops after enough strategies are found; it does not promise
broader-tier exploration once that quota is met.

The old flat ladder's final scaffold-candidate reservation is not applied by
this grouped path. Structural-complexity evidence and per-tier ranking remain,
but strategy-level reservation and complex-target parity need evaluation before
claiming improved strategic synthesis quality.

Next, audit incomplete identities (including retained-skeleton sanitization),
add bounded site-aware scheduling, and freeze a chemistry-stratified untouched
evaluation with strategy/site recovery, invalid rate, useful alternatives,
latency and chemist review. Shared observation-to-operator compilation and
required-versus-ranking context remain a subsequent chemistry-contract change.
Those changes may require operator rebuilding; this stage does not.
