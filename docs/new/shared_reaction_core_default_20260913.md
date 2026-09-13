# Shared reaction core v2: default promotion

Date: 2026-09-13

V2 is the operational default following the user's instruction, “the new system
seem good, make it as the default.” This is a selection of the working foundation,
not an assertion that the frozen numerical gates or independent review have passed.
The earlier [v2 validation report](shared_reaction_core_v2_20260913.md) and its
generated manifests remain the evidence for chemistry performance.

## Use

Restart the previous server, then run from the repository root:

```powershell
python -m app.web_api --workbench
```

`python -m app.web_api` uses the same engine in Condition Desk. Neither launcher
needs the experimental flag. Omit a development-sample `--index` to select the
actual Full library. A previously configured `CONDITION_RECOMMENDER_INDEX` still
selects that custom library; remove it if Full is intended.

Full has 571,157 admitted observations and 561,283 eligible v2 projections;
Compact has 94,643 observations and 92,597 eligible projections. Both derived
artifacts are already built. No source conversion or full dataset rebuild is
needed for this default change. Chemistry definitions and serialized projections
are unchanged.

## Application and build contract

- `GenericConditionRecommender.from_path()` and `recommend_generic_conditions()`
  select shared-core retrieval by default. The desktop recommender and coworker
  shell inherit that loader. The convenience function and recommendation CLI
  default to the Full SQLite index.
- A source named `generic_index.sqlite` uses `generic_index.shared_core.sqlite`.
  Explicit companion paths remain supported. Loading verifies source identity,
  row count and chemistry definition binding. Missing or stale companions fail
  explicitly instead of switching engines.
- Review mode uses the separate review index's companion when that index exists.
  Verified safe reuse of trusted data uses the trusted companion. A source
  canonical JSON/manifest needs a companion built against that loaded source;
  its identity is not interchangeable with a persisted SQLite index's identity.
- Normal artifact builds, saved-batch combination and `generic_index_cli` build
  the companion automatically. Canonical-only conversion still omits indexes.
  For older custom SQLite libraries, build only the derived companion:

```powershell
python -m condition_recommender.shared_core_cli build `
  path/to/generic_index.sqlite path/to/generic_index.shared_core.sqlite `
  --workers 4
```

No API routes or chemistry schemas changed. Capabilities now identify the selected
engine as `shared_reaction_core.v2`. The historical result field
`recommendation_mode: experimental_shared_core` and pending-independent-review
warnings remain intact; a default setting must not rewrite uncertainty provenance.

For explicit comparison/rollback, the web runtime accepts
`CONDITION_SHARED_CORE_EXPERIMENTAL=0`; remove it and restart to restore v2.
Python accepts `use_shared_core=False`, the recommendation CLI accepts `--baseline`,
and the index CLI accepts `--baseline-only`. The historical build parameter
`experimental_shared_core=False` skips projections for a direct baseline build.
Low-level in-memory evaluation retains explicit engine selection and training-only
projections. These controls are temporary comparison facilities: remove the older
retrieval implementation after adjudication and replacement coverage establish
parity for the remaining supported contracts, preserving historical report files.

## Next priorities

1. **Adjudicate the graph-comparison disagreements.** Review the expanded packet
   and baseline-only examples with explicit accept/reject/uncertain decisions.
   Determine which excluded analogues are useful transfers and which involve
   different protected edits. Encode accepted relaxations as versioned graph
   rules with positive, negative and conflicting-evidence regressions. Avoid
   solving chemistry-policy disagreements with ranking weights.
2. **Measure realistic retrieval coverage.** Freeze a new publication/reaction
   grouped test set and a scaffold-disjoint set, train on the remaining larger
   corpus, and compare coverage, recipe recovery, independently judged precision,
   fallback levels, abstentions and latency. Keep the consumed 2,000-row panel
   as development evidence; do not reuse it as untouched validation.
3. **Close confirmed coverage gaps, then consolidate.** Prioritize missing
   observations and source-context abstraction only where the adjudicated cases
   justify them. General multievent edit graphs already exist; multievent source
   abstraction and some upstream state/stereo-only signatures remain incomplete.
   Rebuild derived projections when their chemistry contract changes, reconvert
   source records only when observation/admission changes, and retire the older
   engine once the replacement passes the agreed checks.

Current measured limitation: in the sparse-training frozen panel, grouped
coverage was 67.84% for v2 versus 79.65% for the earlier engine; scaffold-disjoint
coverage was 67.52% versus 77.01%. Scaffold top-five recipe recovery also exceeded
the original allowed regression. These are not Full-training performance estimates.
Default promotion preserves those results and leaves formal adjudication open.

## Verification

Current promotion checks and browser screenshots are under
`results/shared_core_default/`. The workbench smoke test starts without an
experimental activation flag, verifies the selected engine and actual library
counts, and exercises Full cyanation, Full aldehyde reduction and Compact cyanation.

- Full/Compact browser regression: **1 passed in 17.0 seconds**, with no page
  JavaScript errors. Full cyanation and reduction screenshots were inspected.
- Complete `pytest -q`: **1,458 passed in 383.33 seconds**. Coverage includes
  default activation, explicit baseline selection, stale/missing companion
  rejection, review-index pairing, automatic artifact builds and Condition Desk
  protocol generation with v2 evidence labels preserved.
- CLI smoke: the default index command built a paired projection artifact from
  the 68-observation development index; the default recommendation command
  returned three v2 cyanation recipes. This is an integration check, not a
  dataset performance estimate.
- Ruff Python name checks and `git diff --check` pass. No browser application
  source changed; the existing compiled client was used in the browser check.
- The temporary verification server on port 8026 was stopped after testing.
  Restart the existing user server to load the new default.
