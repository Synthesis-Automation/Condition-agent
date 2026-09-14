# Three-view condition retrieval

Condition recommendation now collects precedents through the complete reaction,
reactant-side graph, and product-side graph. These channels feed the same
original-query comparison, condition compatibility checks, source-context
recipe aggregation, and ranking. The new channel is visible as `reactant_side`
in candidate traces and recommendation evidence.

## Chemistry contract

The reactant view projects the before-side of each already validated shared
reaction-core graph. It retains reactive center anchors, atom states and
stereochemistry, bonds, local environments, and within-component center
distances. Qualified source/departure generalization comes from the existing
core levels. Product states and reaction edit events are omitted from this
candidate lookup view. The full before/after graph still decides eligibility.

This lets identical inputs or related reactant cores nominate observed
precedents. A competing transformation of those inputs is excluded when its
full reaction core contradicts the requested transformation. Such a precedent
does not silently supply conditions for a different product. Ar–I, Ar–Br,
mesylate and fragment-source analogues use the existing qualified graph rules;
there is no new list of reaction names or leaving-group substitutions.

The implementation uses persisted forward-facing and product-facing graph
views. It does not execute the forward or retro planner inside a condition
request. Existing retro precedent IDs remain an additional candidate provider.
No prediction is presented as an observed precedent, and no missing input or
product is invented. This change consumes a complete reaction query; partial
inputs still require a proposed transformation before condition recommendation.

## Retrieval and evidence

- Automatic mode retains the existing early stop when direct evidence has
  sufficient independent support. Otherwise it consults both side views.
- **Broader analogues** explicitly consults all three views. The same permitted
  core levels and original-query chemistry gates apply regardless of channel.
- Direct retrieval retains its 512-candidate capacity. The qualified union has
  a 768-candidate cap, reserving room for auxiliary providers. Product-side,
  reactant-side and supplied retro seeds share the remaining capacity in round
  robin order. Duplicate hits attach provenance without consuming a new turn.
- Each observed precedent is compared once and counted once. Multiple channels
  do not increase its evidence support or score. Independent-reference support
  and canonical recipe aggregation retain their existing semantics.
- Candidate truncation remains explicit. Expanded lookup cannot guarantee new
  eligible chemistry: a complete direct index may already contain every
  precedent allowed by the full-core comparison.

No ranking-profile control or chemistry weight was added. The Workbench shows
“Also found through reactant-side reaction matching” on relevant recommendations.

## Artifacts and versions

- Reaction signatures, shared-core comparison and existing core identities are
  unchanged.
- Reactant projection: `shared_core_reactant_projection.v1`, with a validated
  definition and RDKit-bound hash.
- Derived shared-core storage: `2.0`; its manifest now requires the reactant
  projection hash. Older companion artifacts are rejected explicitly.
- Retrieval: `shared_core_retrieval.v3@3.0`; traces also carry the reactant
  projection hash. The earlier v2 definition is retained as historical
  evaluation provenance, not another runtime selection path.

The canonical dataset-building path automatically includes the new lookup keys.
For existing datasets, only the derived companion indexes need rebuilding:

```powershell
python -m condition_recommender.shared_core_cli build datasets/literature/full/generic_index.sqlite datasets/literature/full/generic_index.shared_core.sqlite --workers 8 --resume
python -m condition_recommender.shared_core_cli build datasets/literature/compact/generic_index.sqlite datasets/literature/compact/generic_index.shared_core.sqlite --workers 3 --resume
```

Builds validate source binding and publish atomically. Raw reaction conversion,
condition normalization and retro/forward operator libraries do not need
rebuilding for this change. Restart running applications after updating code
and companion indexes.

## Verification

Focused regressions cover partner/map invariance, qualified source/departure
generalization, preserved stereochemistry and site geometry, excluded product
states, invalid graph colors, reactant-only candidate recovery, wrong-product
rejection, condition constraints, duplicate evidence, channel budgeting, and
stale manifests. Existing ambiguous/conflicting-evidence and artifact-resume
regressions continue to apply.

A before/after engineering comparison uses three previously consumed targets,
both local libraries, and automatic/broad scopes. Baseline retrieval code is
frozen from the pre-change Git commit. Raw responses, query protocol and
summary results are stored under `results/three_view_conditions_20260914/`.
This is a regression check, not an untouched chemical-accuracy evaluation.

Both operational companion indexes are rebuilt: **571,157 Full** and **94,643
Compact** observations. Ordered hashes of every projection's position, source
row hash and payload hash match the pre-change artifacts exactly. Source
identity, eligible counts and unavailable reasons also match. Only the new
lookup keys and their storage metadata were added.

Across all 12 before/after queries, the ranked recipe IDs, supporting precedent
IDs, qualified candidate counts and independent support counts are unchanged.
Automatic mode stops with sufficient direct support for these examples.
Broader mode adds side-view candidates, including competing transformations
that the common comparison excludes:

| Full query, broader scope | Candidates before | Candidates after | Qualified, both runs |
| --- | ---: | ---: | ---: |
| Aryl iodide cyanation | 512 | 520 | 20 |
| Indazole cyanation | 512 | 519 | 20 |
| Acetaldehyde reduction | 512 | 764 | 512 |

This demonstrates working additional retrieval and stable qualification for
these examples, not improved chemical recall. Timings are diagnostic only;
concurrent builds, cache state and execution order were not controlled.

- Final complete Python suite: **1,515 passed in 549.29 seconds**.
- Production frontend build passes; existing bundle-size warnings remain.
- Full/Compact Workbench test passes in 14.4 seconds: actual three-channel API
  traces, bounded and deduplicated candidates, reactant-side evidence text,
  recipe selection and zero page JavaScript errors. Both screenshots were
  inspected.
- Ruff error checks and `git diff --check` pass.

Evidence includes `regression_summary.json`, before/after projection
fingerprints, `full_suite_final.log`, `browser.log`, `frontend_build.log`, and
`workbench_full.png` / `workbench_compact.png` under the result directory above.

The browser regression can be run from `web/reaction_recommender` against a
running Workbench server after both companion indexes are rebuilt:

```powershell
$env:THREE_VIEW_CONDITIONS = '1'
$env:WEBUI_TEST_URL = 'http://127.0.0.1:8000'
$env:BROWSER_CHANNEL = 'msedge'
npx playwright test e2e/three-view-conditions.spec.ts --reporter=line
```
