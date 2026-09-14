# Single-step retrosynthesis foundation follow-up

This follows the accepted [strategy test stage](single_step_strategy_test_stage_20260913.md).
The Workbench continues to use the grouped single-step API. Search now completes
previously missing identities, shares validation work across target sites, and
reserves qualifying scaffold-simplifying strategies.

## Chemistry and search contracts

- `reactive_taxonomy.reaction_identity_graphs` describes edit sites using a
  canonical labeled graph of the complete product and anchored edit events.
  Cleavage and hydrogen changes can locate a site even when no product bond was
  formed. Precursor-only endpoints are external ports, not a list of leaving
  group names. Supplied mapping must pass the existing correspondence checks.
- Retained precursor fragments can receive a graph identity without being
  sanitized as standalone molecules or capped with invented hydrogens. Original
  atom states, bonds and component membership remain observations. An identity
  descriptor is never used as an executable precursor molecule.
- Existing nonempty `SITE1`/`SYN1` identities remain stable. Versioned
  `SITE2`/`SYN2` descriptors complete missing runtime identities. The compiler's
  existing operator and template identities have not changed; existing operator
  artifacts do not need rebuilding for this search change.
- Templates share work across operators. Generated proposals then share the
  validation budget across provisional product sites within each operator.
  These provisional sites schedule work only: forward signature validation and
  operator agreement still gate every returned proposal.
- Scaffold reservation operates on distinct strategy representatives using the
  existing versioned ranking policy, specificity tier and score-band guards.
  It does not promote an arbitrary low-scoring proposal to fill a category.
- `single_step_strategy_search.v2` records the scheduling change. Search
  coverage exposes the number of provisional operator/site groups considered.
  Limits remain per tier; the ladder stops when it has enough strategies.
  Workbench ranking controls use chemistry labels (reaction site, synthon and
  precursor choice) instead of hard-coded identity namespace names.

The identity compatibility rule is temporary. Remove the old runtime descriptor
path when observation-to-operator compilation adopts the same graph descriptors,
the operator/strategy identity migration is versioned together, source replay
and partner-order/symmetry regressions pass, and dependent completion indexes
are regenerated. Until that gate, mixed identity namespaces are explicit and
must not be interpreted as proof that every chemically equivalent realization
has been collapsed into one strategy. This follow-up does not claim a complete
replacement of operator compilation or of required-versus-ranking context.

## Development regression

The four previously consumed development targets are indazole nitrile,
ethylamine, ethanol and biphenyl. All **41 previously incomplete proposals**
(6, 13, 15 and 7 respectively) now receive complete strategy identities while
preserving their operator signatures. New searches have no unresolved strategy
identities for these four examples. These are regression results, not an
independent estimate of chemical usefulness.

Evidence: `results/single_step_foundation_20260914/development_audit.json`.

## Frozen evaluation protocol

The separate panel samples 1,500 Full condition-index observations and freezes
60 queries across structural transformation annotations. Sampling excludes
direct observation, publication and reaction-equivalence overlaps with the
listed earlier shared-core panels. The existing grouped splitter keeps related
publication/reaction groups out of both training and test partitions. This is
a held-out comparison within this experiment; it is not a claim of exclusion
from every historical experiment in the repository.

The split contains 1,203 training and 297 held-out observations. Publication,
canonical-reaction and mapping-equivalence overlap counts are all zero. A
separate product audit also finds zero exact canonical query products in the
training inputs. Scaffold overlap is allowed; this is not a chronological
evaluation or a claim of complete scaffold novelty.

Only the training partition supplies executable templates. Flat and grouped
engines receive the same library, top-five request, 100-template and
30-validation limits per attempted tier. The flat comparator includes the same
identity-completion fix, isolating the search/grouping comparison. Grouped
exact-precursor recovery considers all retained precursor choices and therefore
has a larger realization allowance than five flat hits. Every valid single
product is searched even when its source reaction cannot be compiled. Coverage
includes all selected queries; recovery uses the explicitly reported number of
available reference identities. Source failures remain visible separately.
The comparison excludes optional stock-based precursor realism and the separate
forward-product competition audit. Those Workbench options can change displayed
ranking or add evidence. Mandatory graph/signature checks remain enabled in
both evaluation engines.

Inputs, partitions, protocol and execution code hashes are recorded before
search outcomes. The results are not used to tune this implementation.
Signature-status counts are automated structural checks, not measured chemical
precision. The HTML packet is explicitly an **unblinded diagnostic review**;
human usefulness scores start blank and are preserved when regenerating it.

The initial v1 evaluation incorrectly skipped product search when source
compilation failed. Its inputs and outputs remain in `evaluation/` for audit;
its coverage numbers are superseded. `evaluation_v2/` copies the exact frozen
panel and corrects that measurement error. Search rules and ranking were not
tuned between runs. This panel is now consumed and must not be called a new
untouched set in subsequent work.

The two runs' training-library JSON payloads are identical. Execution manifests
confirm that the search, ranking and identity modules did not change. The
initial trial of build-cache compatibility was removed after the shard audit
showed changed source files; that build-only module is not used by either
evaluation search. Timing is diagnostic only: engines run in fixed order with
shared caches, while a concurrent Full build contends for CPU. No speed
improvement is claimed from these measurements.

A later build-lifecycle fix releases cached whole reaction analyses between
shards. It does not change the compilation, identity or ranking functions. The
audit verifies that removing the newly added cache-cleanup helpers reproduces
the frozen compiler and operator-application source hashes; regression tests also compare compiled JSON
before and after cleanup. This change therefore does not require retuning or
reinterpreting the frozen search comparison.

Forward-library preparation and web artifact freshness were subsequently fixed
and tested separately. They are outside this evaluation's search path, which
does not enable the optional forward-product competition audit.

### Corrected evaluation result

The training-only library admitted 929 of 1,203 source observations. Of the 60
held-out queries, 33 supplied usable compiled reference reactions; the other
27 remain visible with their compilation failures. One disconnected salt-form
product is outside the current single-component target contract. Its component
was not selected or removed implicitly.

| Metric | Flat top five | Grouped top five |
| --- | ---: | ---: |
| Queries returning proposals | 59/60 | 59/60 |
| Total distinct strategies across queries | 230 | 295 |
| Mean distinct strategies per valid single target | 3.90 | 5.00 |
| Observed strategy recovered | 14/33 | 14/33 |
| Observed edit site recovered | 17/33 | 17/33 |
| Exact precursors recovered | 11/33 | 10/33 |
| Returned candidates violating signature-status contract | 0 | 0 |

The demonstrated benefit is more distinct alternatives with unchanged observed
strategy/site recovery. Exact precursor recovery loses one case, even with
grouped alternatives included. This is not an accuracy improvement claim.
These results use a small training-only library, not the production Full
library. Chemical usefulness and feasibility remain unscored pending review.

The corrected report is
`results/single_step_foundation_20260914/evaluation_v2/report.json`.
The [review packet](../../results/single_step_foundation_20260914/evaluation_v2/chemist_review/review_packet.html)
contains all 60 cases and 590 representative cards, with alternate precursors
expandable. Browser inspection found no drawing errors or horizontal overflow.
The blank review sheet is beside the packet. The panel remains consumed after
this evaluation; improvements motivated by its disagreements require a new
prospective evaluation.

```powershell
python -m app.retrosynthesis_evaluation freeze datasets/literature/full/generic_index.sqlite results/NEW_PANEL --size 1500 --queries 60 --exclude results/shared_core_followup/frozen_panel.sqlite --exclude results/shared_core_validation/generic_index.sqlite
python -m app.retrosynthesis_evaluation evaluate results/NEW_PANEL/panel.json
python -m app.retrosynthesis_evaluation review results/NEW_PANEL/report.json
```

Use a new directory for a new experiment. Frozen panels and started executions
cannot be silently overwritten. A review refresh preserves an existing
`chemist_review/review.csv` containing human judgments.

## Full operator artifact

The Full condition index and Full retrosynthesis library are separate artifacts.
This work resumes the missing Full operator build from its source shards.
Reuse requires the existing exact build configuration and source path, size
and timestamp checks. Many saved older shards refer to source files that have
since changed; those must be recompiled. Their manifests are not rewritten to
pretend that stale outputs describe the current source. The source compiler and
its strict `pass_only` admission policy are unchanged.

The long build exposed a memory issue: each worker retained up to 20,000 whole
reaction analyses across source shards. Twelve workers reached about 12 GB
combined and continued growing. Both serial and parallel batch paths now clear
that compilation cache after each shard, including failed jobs, and collect
released objects. Focused tests verify unchanged compiled payloads, cache
cleanup on success/reuse/failure, and resumable merge behavior. Initial checks
after restarting the build showed about 2.6-3.8 GB across twelve workers. A
later measurement at 683 completed shards was 5.57 GB combined; the build then
finished without the earlier memory growth. These are sampled working-set
measurements, not a peak-memory benchmark. Evidence is in `memory_progress.jsonl`.

The resumable build and its progress evidence live in
`results/single_step_foundation_20260914/`. `full_progress.json` reports actual
completed work. After the merge, the local verification helper automatically
checks artifact counts, replays 40 fixed-seed source precedents, and calls the
Full API for indazole nitrile, ethylamine and biphenyl. Results go to
`full_verification.json`, `full_api_smoke.json` and `full_completion.md`.
A `waiting_for_build` status does not mean Full is ready. A process failure is
reported explicitly; it is not turned into a successful release status.

All 725 source shards have compiled and merged. The final artifact contains
660,190 source rows, 458,466 unique admitted observations (69.44%), 10,770
operators, 28,954 realizations and 95,362 templates. Merging removes 3,715
duplicate admitted observations. Rejection reasons remain in the build report;
the largest groups are unverified materialized cores (99,949) and unavailable
atom mappings (76,627). These are compiler admission results, not a comparison
against the separate condition index's coverage.

The matching forward artifact admits 78,081 of 95,362 templates. The remaining
17,281 fail the unchanged independent source forward/reverse replay gate; they
are not silently included in the forward-admitted library. Both source artifacts
are prepared, and **Full is ready for local Workbench use**. The
[build report](../../results/operator_retrosynthesis_poc/full_scale_v3/full/build_report.json)
and [verification report](../../results/single_step_foundation_20260914/full_verification.json)
record the final counts and artifact hash.

Full source replay passes for all 40 fixed-seed samples: 36 exact and four under
the compiler's existing stereo-relaxation policy. The three API targets each
return five distinct strategies, with 19 precursor choices in total. All retain
verified-signature status. The separate advisory forward audit reports three
structurally supported, eight inconclusive and eight out-of-scope choices.
Targeted replay reproduces 18 choices; one generating operator is absent from
the independently forward-admitted library. Blind search has its own coverage
and budget limits, so these API checks do not establish feasibility or complete
forward support for every displayed strategy.

## Try the updated Workbench

The prepared forward library must be at least as recent as its paired retro
artifact. The API previously loaded an older prepared file solely because it
existed; it now derives and caches from the newer retro source instead. Release
preparation builds both artifacts together so the default forward audit does
not require this expensive preparation on the first request. This timestamp
guard detects the local stale-artifact case; it is not a content-fingerprint
provenance guarantee for arbitrarily copied files.

Offline forward preparation supports `python -m forward_synthesis build-library
SOURCE OUTPUT --workers 12`. Serial and parallel builds retain the same source
forward/reverse replay gate, rejection reasons, ordering, support and precedent
aggregation. A parity regression covers passing, failing, malformed and
duplicate source operators. Worker compilation caches are released between
checks to bound retained memory.

```powershell
python -m app.web_api --workbench --build
```

Select **Single-step retrosynthesis**, request five strategies, and inspect
**Search coverage** and each strategy's **Precursor choice**. The frontend build
must be refreshed and an already running Python server restarted to load the
updated code. The separate source condition dataset does not need conversion
again for this change.

## Verification

- Complete Python suite after the cache-lifecycle, forward-build and freshness
  fixes: **1,501 passed in 624.33 seconds**.
- The final production frontend build passes. Existing bundle-size warnings
  remain; no asset-loading failure was observed.
- The Workbench browser check passes with its default forward audit enabled:
  five distinct Compact strategies,
  selection of an alternate precursor, its own condition request, persistent
  selection while evidence loads, search coverage and zero JavaScript errors.
- The same browser workflow passes on Full using biphenyl, including its
  alternate precursor and resulting conditions. The initial Full test assumed
  ethylamine would have alternate choices, as Compact does; its five Full
  strategies each had one choice. The browser fixture was changed to a target
  with alternates; no production ranking or admission rule was changed to make
  that assertion pass. Both screenshots were inspected after the passing runs.
- The post-build verification workflow was exercised on Compact: all 40 sampled
  precedents replayed under the existing compiler acceptance policy (37 exact,
  three with the compiler's existing stereo relaxation). Its three API queries
  each returned five distinct strategies with verified-signature status and
  complete forward-audit status counts. The matching prepared forward artifact
  admits 26,325 of 31,980 source templates; 5,655 fail the unchanged independent
  source forward/reverse replay gate. The previous stale artifact had 31,791
  source templates and 26,156 admitted operators.
  Full verification is recorded separately above.
- Ruff Python error checks pass. Generated libraries, evaluation outputs and
  screenshots remain local under `results/` and are not source-code additions.

Useful artifacts include `workbench_compact_default.png`, `workbench_full.png`,
`pytest_release_final.log`, `browser_compact_default.log`,
`browser_full_default.log`, `frontend_build_final.log` and
`compact_verifier_smoke/full_verification.json` in
`results/single_step_foundation_20260914/`.

## Remaining chemistry research

This release completes the implementation and local artifact-verification
stage. A chemist usefulness review and a new prospective evaluation remain
necessary before claiming stronger synthesis quality. Priorities are shared
observation-to-operator identity compilation, an explicit distinction between
required chemistry and ranking context, and better blind-forward search
coverage. Source observations also need step-granularity evidence: a valid net
graph transformation can summarize a sequence of operations. None of these
uncertainties is resolved merely by returning five graph-validated proposals.
