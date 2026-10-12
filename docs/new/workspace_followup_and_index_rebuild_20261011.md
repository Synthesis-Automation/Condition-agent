# Workspace follow-up fixes and index rebuild

The completed Normal run `8373df5f19c04ca59572f2fc003b99c1` exposed two remaining
integration errors: optional corpus loading aborted `assess_retro_validity`, and
route authoring omitted saved inspections for reported steps although publication
required them. Its final answer was repaired manually and completed successfully.

## Fixes

- `assess_retro_validity` contract 4 retains structural/template assessment when
  loading optional condition/shared-core data raises a missing-file or validation
  error. `corpus_availability` preserves the error; artifact warnings and both
  summary views disclose the gap. Unexpected errors still fail. A structural
  check does not establish experimental feasibility.
- Route authoring reuses exact matching saved inspections for both proposed and
  reported steps. Reported steps still require explicitly supplied source evidence;
  automatic inspection attachment cannot create a reported claim. Conditions retain
  their explicit source attribution, without borrowing the automatic inspection.
- Prompt/help guidance gives the supported fragment limits and topology values,
  explains target graph mismatches and parent call references, and identifies the
  actual artifact-inspection parameters. The environment prompt provides a portable
  UTC clock alternative to the unsupported PowerShell 5.1 `Get-Date -AsUTC` option.

The saved 17-step input now passes answer/evidence validation without the manual
reference patch. The original failed validity candidate completes its structural
assessment with the stale corpus explicitly unavailable. Historical evidence is
unchanged; the current-code scientific call is in a fresh pinned investigation:
[`verification.json`](../../results/ai_native/workspace_followup_fix_20261011/verification.json).

Validation completed on October 11:

- Targeted workspace/authoring/prompt regressions: **97 passed**.
- Full `python -m pytest -q`: **2,999 passed, 4 skipped**, in 1,200.18 seconds.
- A separate probe rebuilt the first completed 500-record shard, loaded compatible
  condition/shared-core indexes, and completed recorded precedent retrieval and
  condition recommendation. Admission, chemistry, condition and family counts
  matched the old shard. This is sample evidence only; full-corpus comparison
  remains pending. See
  [`shard_probe/report.json`](../../results/ai_native/workspace_followup_fix_20261011/shard_probe/report.json).

## Full rebuild

The user authorized the rebuild on October 11. The existing unified builder is
regenerating all **2,054,437** physical observations from
`datasets/intermediate_datasets`, with eight workers and 500-row shards. It builds
the derived condition/shared-core indexes, evidence and fragment stores, and
operator artifacts as one coherent release. It never changes old index version
labels to pass compatibility checks.

- Previous active release: `bd7a78d653af3156ba57d2d2`.
- New build directory: `datasets/processed_datasets/releases/05fd3b97e436c3ba79147e26`.
- Driver: `results/dataset_builds/20261011_workspace_index_refresh/run_build.py`.
- Durable progress: `progress.jsonl` beside the driver.
- Build result: `result.json`, written only after the builder completes and publishes.
- Post-build checks: `verify_when_ready.py`, with `verification_status.json` and
  `verification_report.json` beside it.

This document records initiation, **not completion**. Inspect those status files
and the active manifest for the current outcome. The old release stays active until
the unified builder validates complete source coverage, artifact hashes and index
dependencies and atomically publishes the new manifest. Completed shards and stages
are resumable with the same driver; do not run duplicate builders concurrently.

The automatic follow-up pins the published release in a new scientific workspace,
checks metadata compatibility, executes source-record retrieval, condition
recommendation and the target's validity assessment, and records before/after
coverage, admission, rejection and family counts. It reports any failure explicitly.
A development-only first-shard probe is separate from full-corpus validation.

Restart the scientific-chat server and start a new Normal investigation to use the
code fixes. After publication, restart consumers that pinned the previous release;
existing investigations retain their frozen data. No HTTP endpoint changed.
Independent chemistry review, untouched evaluation and experimental route validation
remain pending; a technically complete build does not satisfy those gates.

## Resumed after editor closure (October 12)

The original builder and verification watcher stopped when VS Code was closed.
The fragment scan's last logged count was 1,158,000 observations. Conversion,
route/evidence catalogs, condition index and shared-core index had completed.
The fragment builder cannot resume an incomplete observation scan, so that stage
must restart; completed stages are reused after their hashes are checked.

At 23:13 Shanghai time, the same release identity was confirmed and the driver
was relaunched with `--finish-only`, skipping conversion. Both the driver and
verification watcher now run through hidden Windows `Start-Process` launches,
independently of the VS Code terminal. Current PIDs are recorded in `jobs.json`;
restart output is in `resume.stdout.log`, `resume.stderr.log`, and the corresponding
`verify_resume` logs beside the driver. Publication and verification remain pending.
