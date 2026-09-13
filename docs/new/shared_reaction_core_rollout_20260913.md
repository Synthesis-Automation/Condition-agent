# Shared-core larger-corpus validation — 2026-09-13

Historical v1 report. See the [v2 general edit-graph beta](shared_reaction_core_v2_20260913.md)
for the current implementation and rebuilt Full/Compact artifacts. The counts
and hashes below describe the earlier v1 artifacts.

The user's feedback, “the review look ok, do the next steps”, was recorded as
qualitative acceptance of the development packet. The per-case review CSV is
preserved; no individual decisions or formal adjudication were invented.

## Decision

Keep shared-core retrieval opt-in. The fresh engineering holdout shows substantial
coverage and historical recipe-recovery regressions. This is a measured reason
to postpone default replacement, even though larger artifacts are useful for
testing cyanation and other supported joins.

| Frozen split | Queries | Existing coverage | Shared-core coverage | Existing top-5 recipe recovery | Shared-core top-5 recipe recovery |
| --- | ---: | ---: | ---: | ---: | ---: |
| Publication/reaction grouped | 201 | 73.63% | 30.85% | 6.97% | 1.00% |
| Scaffold disjoint | 154 | 74.03% | 20.78% | 6.49% | 3.25% |

Both paths have zero reported hard compatibility violations and complete
explanations for covered queries under the existing automated checks. Neither
finding establishes chemical transferability or experimental success.

The cohort contains 1,000 records, with development observations, publications,
canonical reactions and mapping-equivalent reactions excluded. Input identities,
split membership, seed, code/definition hashes and metrics were frozen before
evaluation. Both paths use identical train/test partitions. Shared projections
contain training records only. Publication/reaction overlap is zero in both
splits; scaffold overlap is also zero in the scaffold-disjoint split.

Limitations: these are 799/846-record training indices, not Full-corpus retrieval
evaluations. The sampling script retained the first 1,000 eligible positions from
a sorted 1,200-position random draw, which biases selection toward earlier source
positions. The two splits reuse one cohort. These results are engineering evidence,
not a formal untouched release panel. This cohort must not become a tuning target
while still being described as untouched.

Reproducible inputs, scripts and reports are under
`results/shared_core_rollout/`: `frozen_validation_manifest.json`,
`freeze_validation.py`, `holdout_comparison.json`, `summarize_holdout.py`, and
the four split-specific evaluation directories.

## What the coverage gap means

In the grouped test split, 62 queries qualify for broad L1/L2 projections,
136 retain observed-local L0 only, and three lack a usable projection. The
scaffold-disjoint split has 33 broad-eligible, 120 L0-only and one unavailable
query. The implementation's broad contract currently models a qualified single
intermolecular join; it is not yet a general abstraction of all graph changes.

The largest grouped-query reasons for withholding broad projections are
unsupported edit patterns (67), protected state changes (24), multiple events
(22), and missing qualified ports at both endpoints (16). These protections
must not simply be disabled to improve coverage.

The next chemistry contract should:

1. Represent the connected reaction-center edit graph for formed, broken,
   order-changed and hydrogen/state edits, without mandatory partner roles.
2. Preserve edit topology, atom types and chemically relevant state changes
   while relaxing surrounding molecular context in explicit levels.
3. Treat departure/source-port generalization as one independently qualified
   graph operator, rather than the prerequisite for every broad projection.
4. Use the same comparison gate for direct, product-side and retro-seeded
   candidates; extra retrieval channels must not bypass graph compatibility.
5. Develop positive, negative, ambiguous and conflicting cases outside this
   consumed cohort, then freeze a new independent panel before default cutover.

## Engineering changes

Completed derived artifacts:

| Library | Source rows | L0 eligible | L1/L2 eligible | Size |
| --- | ---: | ---: | ---: | ---: |
| Full | 571,157 | 561,283 | 178,347 | 2,058,027,008 bytes |
| Compact | 94,643 | 92,597 | 28,171 | 339,681,280 bytes |

Full completed in 830.95 seconds using six workers; Compact in 290.27 seconds
using three workers, with both builds overlapping on this machine. These are
local build timings, not a controlled benchmark. The chemistry definition hash
remains `SCD1:06028c3a041943b28783e3413ab6bc068da1e4f821befe5fd5db044575e81686`.
Both artifacts pass SQLite integrity and foreign-key checks. Their source IDs
match the canonical indices; 128 sampled projections per artifact also pass
row-binding and payload-hash verification. The audit scripts and full manifests
are in `full_artifact_audit.json`, `compact_artifact_audit.json` and
`audit_artifacts.py` under the rollout results directory.

For `Brc1c2c(cccc2)[nH]n1.N#C>>N#CC1=NNC2=CC=CC=C12` and the corresponding
iodide/HCN query, Full now returns **five recommendations** at top-k five from
20 qualifying observations and 14 independent evidence units. Automatic search
examines 20 candidates; broad search examines 102 and retains the same 20.
Compact retains four observations from three independent evidence units and
returns four recommendations. These are analogue precedents with explicit
source/input differences; five is an output cap, not a requirement to fill.

The full-index headless Edge/Playwright test passed: actual Full/Compact counts,
five returned cyanation recommendations, L1/L2 display, and no page JavaScript
errors. Screenshot: `results/shared_core_rollout/workbench_full.png`.

The projection builder supports bounded worker processes, deterministic
source-order writes, batch checkpoints and resumption. A completed artifact is
published atomically. Checkpoints are bound to source identity and chemistry
definitions and cannot be loaded as completed artifacts. Serial/parallel equality,
cancel/resume behavior and mismatched-source rejection have regression coverage.

```powershell
python -m condition_recommender.shared_core_cli build `
  datasets/literature/full/generic_index.sqlite `
  datasets/literature/full/generic_index.shared_core.sqlite `
  --workers 6 --resume `
  --progress-file results/shared_core_rollout/full_progress.jsonl
```

Offline evaluation accepts `--experimental-shared-core`; its projection artifact
is built exclusively from the training partition. The workbench capability
response exposes the configured index's stored identity and record count.
The selector and status display the actual count; a custom development index
is labeled “Custom index” instead of “Full”. These are additive API fields;
recommendation request contracts and chemistry definitions are unchanged.

## Running the larger experimental workbench

After both derived artifacts exist, stop the old server in its terminal and run:

```powershell
$env:CONDITION_SHARED_CORE_EXPERIMENTAL = '1'
python -m app.web_api --workbench
```

Use the default dataset paths; an old `--index results/shared_core_validation/...`
argument explicitly selects the small development sample. To return to the
default retrieval path, unset `CONDITION_SHARED_CORE_EXPERIMENTAL` and restart.
The flag is process-wide, so use a terminal dedicated to this experiment.

Full source conversion, admission and canonical recipes have not changed. The
new SQLite files are derived projections of existing admitted observations.
A raw-data rebuild is not required to use these artifacts.

## Final validation

- Complete `pytest -q`: **1,437 passed in 462.18 seconds**.
- Frontend `npm run build`: passed in 30.81 seconds; existing bundle-size warning.
- Full-index Playwright: **1 passed in 12.8 seconds**, with no page JavaScript errors.
- Modified Python modules pass Ruff undefined/unused-name checks; `git diff --check` passes.
- Full and Compact completed artifact integrity, source binding and sampled row
  verification pass. Four frozen evaluation runs completed without unhandled failures.

Logs, hashes and final status are recorded in
`results/shared_core_rollout/rollout_manifest.json`. The temporary port-8026
browser-test server was stopped. The user's existing server was not restarted;
the commands above activate the larger experimental artifacts after a restart.
