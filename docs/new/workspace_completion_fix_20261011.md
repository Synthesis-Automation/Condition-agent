# Workspace completion fix — 2026-10-11

The Normal run `2b7e036859cf494b94f814face855488` reached its 1,800-second
deadline without a final receipt. Its saved findings-only draft was readable,
but repeated step inspections failed because optional condition enrichment tried
to load an incompatible recommendation index. The earlier molecular-motif update
changed its fallback descriptor definition identity. Presence-only readiness did
not disclose this at startup.

## Changes

- `capabilities().index_compatibility` checks SQLite condition-index metadata using
  the canonical lazy loader. It reads metadata without scanning corpus rows. Missing,
  incompatible and merely present artifacts remain distinct. Non-SQLite artifacts
  are explicitly not checked. Compatibility does not certify corpus content.
- Non-pure workspace prompts disclose the condition-index incompatibility before
  research, explain which evidence remains usable and discourage unchanged retries.
- `inspect_step_precedents` contract **5** preserves independently loaded template,
  reference and procedure evidence when optional condition enrichment raises a
  missing-file or validation error. It records the error, empty observations and
  `condition_transfer: not_assessed`. Compact and detailed summaries preserve this
  gap. Required source/selection errors and unexpected runtime errors still fail.
- Structured runtime prompts give a deadline and reserve up to two minutes for
  finalization. Free-text modes retain their existing prompt contract.
- After the runtime stops a timed-out structured attempt, the service may recover
  only that attempt's `answer-draft.json`. Normal file-integrity, baseline, schema
  and evidence validation must pass. The answer remains `partial: true`, unreviewed,
  and the turn remains `timed_out`, with the original error. No completion receipt
  is fabricated. Cancellation, missing/invalid drafts and older repair-attempt
  drafts are not recovered. The browser displays an incomplete-result notice.

No API endpoint or chemistry admission rule changed. Historical investigations and
their outputs are not rewritten. Recovery applies to newly executed attempts, not
retroactively to every historical timeout. Schema/evidence checks do not establish
chemical feasibility or semantic support for an agent's claims.

## Use and remaining data migration

Restart the scientific-chat server, refresh the browser and create a new **Normal**
investigation. Do not change a pinned investigation's baseline to continue with new
scientific code. Source-led planning can continue with explicit condition-data gaps.

The configured `bd7a78d653af3156ba57d2d2` release still carries older fallback
descriptor definitions. This patch does **not** refresh its 1.7 million reaction
rows or make condition recommendations available from it. A versioned canonical
corpus refresh/rebuild and normal validation are required; changing metadata alone
would misrepresent stale descriptors and is not a migration.

## Validation

Regression tests cover stale/invalid index preflight, preservation of actual source
evidence and answer binding, unexpected failures, valid/invalid/missing timeout
drafts, repair-attempt isolation, cancellation, a real child-process timeout,
free-text mode isolation, and the browser's incomplete banner.

The local production-data regression is recorded in
[`workspace_completion_fix_20261011/review.html`](../../results/ai_native/workspace_completion_fix_20261011/review.html).
It re-runs the original target disconnection and formerly failing selected-step
inspection in a fresh pinned workspace, then validates a proposed-step answer.
It is an engineering regression, not a new model trial or chemistry release gate.
