# Scientific workspace implementation status

Date: 2026-09-26

The local scientific workspace is implemented in
[`chem_coworker/scientific_workspace`](../../chem_coworker/scientific_workspace).
It provides Python and CLI access to canonical operations, a content-identified
development baseline, immutable result artifacts, append-only notes and history,
fresh-session recovery, and deterministic replay. Its scientific operations remain
model-free. An optional Codex runtime adapter and local browser conversation
interface now answer natural-language questions using this same workspace.
Neither introduces new chemistry rules or a parallel recommendation engine.

See the [quickstart](Scientific_Workspace_Quickstart.md) for runnable commands.

## What is implemented

- Baseline recording: Git revision and worktree status, code/definition hashes,
  Python and chemistry dependency versions, registry validation, full checksums
  and available SQLite metadata for explicitly configured artifacts.
- Scientific access: reaction and molecule analysis, canonical condition
  recommendation, paginated indexed precedent records, procedure lookup, recipe
  normalization and assessment, route planning and issue-backed revision.
- Investigation history: objectives and constraints, full calls and results,
  structured execution errors, evidence references, hypotheses, decisions,
  open questions, attached scripts, and lifecycle transitions.
- Recovery: reopen without chat history; reject changed code or input identities;
  rehash data for replay; compare full results; rehydrate route objects by replay
  before a fresh session applies a revision.
- Tests: operation parity with domain outputs, provenance preservation,
  missing inputs, pagination, repeated procedure observations, corrupted
  artifacts, write conflicts, baseline drift, cancellation, and route revision.
- User conversations: `/scientific`, questions and follow-ups, saved thread
  identity, asynchronous progress, cancellation, runtime deadline, and linked
  evidence inspection. Runtime prompts and event logs are saved per turn.
- Concise chat interface: responsive history sidebar, bottom composer, inline
  state/elapsed time and expandable tool activity. Structures, scientific details,
  sources, and uncertainty are collapsed beneath the answer. Active work remains
  visible and cancellable from another chat; drafts survive navigation in the page.
  Progress reports actual recorded events; answers appear after final validation.
  Dark and light themes are available from the header, with a saved browser
  preference and dark mode as the default. Chemistry drawings keep their original
  element colors on light canvases in both themes.
- Optional runtime: native `codex exec` with structured final responses and exact
  thread resumption. Codex supplies the iterative tool loop and programmable
  environment; the scientific packages do not depend on its model service.
- Answer boundary: strict response schema, citation existence/type checks,
  explicit unreviewed status, and no-local-evidence status where applicable.
  These are traceability checks, not automatic verification of scientific prose.
- Structured answers: `scientific_answer.v2` links explicit molecules, reaction
  steps, routes, per-condition/yield attribution, claims, and captured sources.
  The UI renders reaction schemes and route dependencies and flags missing target
  products. Legacy saved answers remain readable without inferred route migration.
- Condition investigations: selected-precedent structural comparisons, canonical
  recipe assessment, source/procedure links, reference counts and missing fields;
  attributed recipe adaptations remain unreviewed proposals.
- Proposed-route investigations: recorded external step/route assessments, conversion
  of planner results into editable proposals, step/intermediate inspection, explicit
  branch replacement/extension, and same-target evidence comparisons. Revisions retain
  their source, declared unavailable-material constraints and assessment settings,
  and rerun all step/topology checks. Supplied recipes use canonical assessment.
- Recorded custom Python execution: snapshot code and cited inputs, deadline,
  JSON output and execution/log provenance. No automatic arbitrary-script replay.
- Answer correction: one same-thread retry after schema/evidence-reference failure,
  preserving rejected drafts, traces and usage. Other runtime failures are not retried.
- Local API: opt-in research profile, loopback launcher, same-origin/token checks,
  server-owned executable/model/data configuration, and one active local turn.

The existing agent remains the investigator. The original starter/comparison
scripts are deterministic utilities. The conversation feature uses a mature
external runtime rather than treating those utilities as an agent.

## User-testable conversation milestone

Launch with `python -m app.web_api --scientific-chat --port 8011` and open
`http://127.0.0.1:8011/scientific`. The [quickstart](Scientific_Workspace_Quickstart.md)
documents configuration, API routes, evidence inspection, and the live smoke test.

Development testing exercised an actual authenticated Codex runtime, not a mocked
model: a natural-language request for `CCBr.N>>CCN` ran recorded reaction analysis
and returned the C–Br cleavage/C–N formation interpretation with the inferred
correspondence caveat. A follow-up used the same thread and saved evidence to
distinguish structural interpretation from experimental feasibility. This validates
the question-to-answer integration; it is not a comparative intelligence result.

The initial trial exposed excessive full-result console output, so the prompt now
asks for selected fields and saved artifacts. Concurrent polling also exposed a
Windows file replacement sharing violation; atomic JSON writes now retry that
specific transient error while preserving atomicity.

The first trial is saved under
`results/ai_native/conversations/e7ee49317308447d9ec9dd777d9c6ba8/`.
Its pre-fix baseline is preserved; later source edits require a fresh investigation.
The reproducible live smoke report is written to
`results/ai_native/chat_smoke_report.json` when explicitly requested.
The post-fix two-turn smoke test passed in conversation
`d1424eb69bb54af09795c03cae03fd23` and produced that report.

A separate live condition-comparison conversation,
`707bec3963e6474c9c424aad99d99ce5`, successfully called canonical Compact retrieval
with `top_k=2` and `search_scope=broad`. It compared the two recipes, distinguished
20 observations/one independent reference from three observations/two independent
references, reported missing temperature/time, and retained review warnings.
The answer cited the complete call and an attached selected-field comparison.
This confirms dataset access through the agent runtime; condition-transfer
validity still requires the independent review described below. The underlying route
and revision workspace operations retain their earlier tests and development pilots.

Browser automation could not connect to the local browser in this session, so
visual acceptance remains pending. HTTP submission, polling, saved-answer retrieval,
and evidence endpoints were exercised against the running application and real
runtime. A browser page is included; this is still a local development preview.

## Development investigations

### Condition-investigation milestone

The next development phase adds `inspect_condition_precedents` and
`propose_condition_adaptation` to the existing recorded workspace. Domain comparison
lives in `condition_recommender`; procedure linkage and proposal history remain
application composition. The existing recommendation ranking and source data are
unchanged. Inspection counts apply to selected observations; exact observation IDs
link procedures, while reaction-only reports remain explicitly unassigned.

Direct recipe assessment now returns `reaction_recipe_assessment.v2`: unresolved
structure is `unknown`, invalid input is separate, and actual conflicts retain
their rule evidence. Observed unchanged functional groups now reach the existing
compatibility engine; previously direct assessment passed only the identity
signature and could omit spectator conflicts. See the
[contract/migration note](../new/condition_investigation_assessment_20260925.md).

An adaptation preserves the original and complete proposed recipe, every changed
field and reason, supporting artifacts, assumptions, risks and compatibility.
It always reports `transfer_status="not_established"`. Evidence-link validation
does not decide whether the cited source scientifically justifies the change.
Stage/declared-absence editing is explicitly unsupported instead of lossy.

Custom scripts can use `run_python` to snapshot code and JSON inputs, capture
output and exit status, and support a `computed` answer label. This is a trusted
local runner, not OS isolation or independent validation. Scripts are not
automatically replayed. Answer validation now allows one correction attempt
using the same runtime thread; rejected drafts remain visible in the history.
A second invalid answer fails. Cancellation, baseline drift and runtime failures
are not repaired; deadlines apply per runtime attempt.

Regression coverage includes observation/procedure separation, recipe attribution,
missing operating values, independent-reference limitations, recorded custom-code
success/failure/deadline/cancellation/input mutation, and bounded answer correction.
Concurrent testing also exposed a Windows transient sharing violation while
reading replaced state JSON; reads now have a bounded retry, without retrying
malformed JSON.

A real model-backed local-corpus trial completed on its first attempt in
conversation `f7c1a62e8a504a37bca7516bcf1ed0ef`. It retrieved two recipes, inspected
eight selected observations spanning three distinct reference IDs, and ran a
recorded Python count over the saved inspection. All eight lacked temperature,
time, concentration and atmosphere; the selected procedure catalog returned no
matching records. The answer rendered a sourced Markdown comparison table,
distinguished the selected counts from retrieval aggregates, and declined to
invent an adaptation. Seven source links were checked through the API. The local
report is `results/ai_native/condition_investigation_trial_report.json`.

This is an integration/development result, not an independently validated transfer.
The pilot baseline predates final custom-runner input-copy/non-finite-output
hardening; its evidence remains readable and a new investigation is required for
further work with current code. Those edge cases are covered by regression tests.

### Structured-answer milestone

New turns must return the versioned answer contract from `answer_contracts.py`.
Reported and computed labels require source IDs; computed statements require a
completed local call or replay. External sources require a saved excerpt attachment,
URL and locator. References are checked across all structured fields, including
sources not repeated in prose. Citation validity remains distinct from claim support.

Explicit molecule IDs connect reactants/products and route dependencies. The service
rejects duplicate/unknown IDs, cycles, missing dependencies and declared dependencies
without a shared intermediate ID. The interface retains invalid notation as an
inspectable drawing error and labels computed, reported, proposed, unknown and input
information separately. It does not present a drawable scheme as chemically verified.

A tafamidis regression fixture based on the user's quoted excerpt depicts the
acylation and ring-closure steps, with independently attributed conditions/yields.
Its constructed intermediate SMILES are labeled proposed in the fixture. The route
stops at the ester, and the UI flags the unproduced free-acid target. This fixture
tests presentation and provenance contracts, not patent accuracy or synthesis validity.

A live Codex trial in conversation `541e447aa5b4468caef3ee79c742750e` returned
five molecules, three steps, one route and five sources for the same target.
The API produced nine valid SVG documents, all five captured-evidence endpoints
were readable, and the saved answer artifact was unchanged by presentation.
The agent retained unresolved structural-analysis warnings for two transformations.
This exercised a source-based multistep answer, not autonomous route revision or
independent verification of the patent's chemistry.

The first turn was rejected because a derived-file attachment alone supported a
`computed` claim. An explicit correction in a follow-up produced the accepted
answer; the rejected turn and evidence remain saved. The prompt now explains this
boundary and asks the agent to validate its draft with the same contract/evidence
validators before submission. At that milestone there was no automatic server
retry; claim-level scientific review remains pending. The local API check report is
`results/ai_native/structured_answer_trial_report.json`.

After the prompt update, a fresh live reaction trial in conversation
`7ea4d1fb67df40708a7a6f58673c25c7` completed on its first turn. Its runtime trace
records a successful draft schema/evidence validation before the final answer.
The API returned three molecule drawings and one reaction scheme for
`CCBr.N>>CCN`, retaining missing conditions and yield. The compact check report is
`results/ai_native/structured_reaction_smoke_report.json`.

### Conversation presentation update

Saved answers now render Markdown tables, emphasis, headings, lists, code, and
external/evidence links. HTML and remote images supplied by a message are disabled.
Tables have a horizontally scrollable container for narrow screens. Questions and
answers display parseable SMILES as SVG cards generated by the existing visualization
package, retaining the exact notation and offering SVG downloads.

The user's saved tafamidis conversation was checked against the live API: one table,
one target drawing in the question, and two starting-material drawings in the answer
were returned. Its answer artifact reference remained unchanged. No model rerun or
scientific baseline change was needed. Seven new tests cover this example, fenced
SMILES, invalid notation, untrusted Markdown, scoped evidence links, immutable view
projection, and API integration. Browser automation remains unavailable, so visual
inspection of the live browser could not be completed.

The original scientific pilot's full local record is under
`results/ai_native/implementation_pilot_01/`. It is ignored by Git because it
contains local corpus evidence and generated artifacts. The manifest identifies
the exact scientific code and selected data; copied scripts and full results
are retained alongside the decision history.

### Condition transfer

Authored query:

```text
Clc1c(C)cccc1C.OB(O)c1ccccc1>>Cc1cccc(C)c1-c1ccccc1
```

The canonical recommender returned three candidate recipes from the Compact
corpus. The agent inspected ten indexed precedent records and ran a saved graph
and evidence-coverage comparison. The first recipe's reported support was 20,
but its independent evidence count was one. The other candidates had support
counts of three and two, with independent counts of two and one respectively.

The representative precursor graphs differ from the query. All three recipe
records lack temperature, time, concentration, and atmosphere; the procedure
catalog returned no matching record for the ten requested reaction IDs. These
are coverage findings about this query and selected snapshot, not proof that
the source publications lack the information.

The recorded decision retains the whole recipes as analogue evidence and
requests additional source evidence before a specific adaptation. No operating
details, mechanistic conclusion, or experimentally validated transfer are
invented. Existing shared-core uncertainty and independent-review warnings are
retained in the authoritative result.

### Route investigation

The initial acetamido-biaryl development target produced three one-step routes
and no typed refinement issues. This result was retained, and the agent recorded
why it did not exercise multistep revision.

A second authored target adds a substituted morpholinoethoxy aryl group:

```text
CC(=O)Nc1ccc(-c2ccc(OCCN3CCOCC3)c(C)c2)cc1
```

The bounded search returned a one-step route, two solved three-step routes, and
three partial routes. Here, "solved" means the planner's terminal predicate was
satisfied; it does not establish experimental feasibility or purchasability.
The three-step alternatives carry a selectivity warning at the final synthetic
step. Existing condition-selectivity repair proposals were unavailable because
the condition support was fallback-only or below the compatibility threshold.

The agent selected an alternate-disconnection investigation using the actual
route, step, and issue IDs. The workspace preserves the original result and
requires the domain refinement and verification checks on alternatives.

The refinement call completed with `improved_alternative_found` relative to the
selected three-step source route. It retained the one-step route with no detected
issue of the selected kind. That route was already present in the original
search; this demonstrates issue-guided selection and preserved lineage, not
discovery of a new superior route or improvement over the baseline's best route.
New three-step alternatives had more issues and were not accepted as improvements.
Experimental feasibility and the adequacy of the terminal-material evidence
remain review questions.

### Agent-proposed route revision milestone

The workspace now exposes the canonical external-step and complete-route assessors
through `assess_route_step` and `assess_route_proposal`. `prepare_route_proposal`
converts a selected planner tree without trusting its prior operator annotations.
`inspect_route_step` provides saved gates, structural neighbors, molecular audits
and supplied-recipe assessment. `revise_route_branch` explicitly replaces or adds
steps, preserves the original artifact, and reassesses the entire resulting route.
`compare_route_proposals` compares the same target under the same assessment options
and declared material constraints without selecting an automatic winner.

The new application records are `route_step_investigation.v1`,
`route_investigation.v1`, `route_step_inspection.v1`, and `route_comparison.v1`.
Core branch edits and declared-material results have their own versioned contracts.
Existing chemistry gates, signatures, operator admission and source datasets are
unchanged. No new chemistry rule or recommendation path is introduced. Invalid or
unsupported proposals remain inspectable; no reviewed/experimental status is granted.

An authored local-corpus pilot under `results/ai_native/route_revision_20260926/`
assessed an N-ethylbenzamide sketch with ethylamine declared unavailable, then added
a proposed ethylamine preparation from acetaldehyde and ammonia. The canonical
assessor returned `admitted_review_only` for both sketches. The declared-material
check changed from `violated` to `satisfied_for_declared_constraints`; actual stock
remained `not_assessed`. The unchanged amide step was assessed again, and the revised
record replayed identically from a fresh workspace instance. This is a predefined
development example, not autonomous discovery, independent chemistry validation,
or evidence of experimental feasibility. The executable example is
[`route_revision_pilot.py`](../../examples/ai_native/route_revision_pilot.py).

Twenty-two new regression cases cover canonical parity, invalid/ambiguous/conflicting
steps, explicit branch edits, broken topology, source preservation, inherited
constraints, fresh-session replay, planner-tree conversion, supplied-recipe checks
and optional forward assessment. The focused workspace/conversation run passed all
59 cases. The complete suite returned **1,644 passed, 2 existing registry failures**
in 333.64 seconds. Lint passed for the new modules, workspace, tests and example;
the package `__init__.py` retains two pre-existing unused evaluation-version imports.

The first live model-backed trial (`989d3225b9f94e8990712871b8945945`) stopped before
any recorded scientific call because the native Codex runtime reported a missing
session record. That failure is retained and is not counted as a successful agent
investigation. Direct operation, persistence and replay checks passed separately.

A fresh trial (`e0c594054adc4a23a33ec3d02c19f128`) selected both optional forward
and condition checks. Its first combined assessment ran for several minutes without
producing a completed call; it was cancelled through the API to keep the integration
trial bounded. Cancellation terminated its Python worker and freed the chat service.
Full-corpus latency for these combined optional checks remains a limitation; this
run is not counted as a successful route investigation.

The bounded live retry (`94e6eb9ab8464366ad5c6e79c0ab5564`) completed successfully
using the default assessment flags (`include_forward=False`, `include_conditions=False`).
The actual agent performed six recorded calls: original assessment, original step
inspection, branch revision, both revised-step inspections, and route comparison.
The original assessment took 4.985 seconds and the revision 5.203 seconds. Its answer
contained a before/after table and two explicit routes; all six cited artifacts were
readable through the API and all ten generated SVGs parsed successfully. The answer
distinguished analogue precedents from exact-substrate experimental evidence and
retained partial-mapping warnings, unrun checks, unspecified operating conditions and
unknown stock. No server answer-repair attempt was required. The report is
`results/ai_native/route_chat_trial_report.json`. This validates the bounded local
question-to-revision workflow, not independent scientific accuracy or comparative
agent superiority. Live browser visual acceptance remains pending.

Scientific source changes require a new chat/workspace baseline. Existing saved
answers remain readable; an old investigation is not silently migrated to new code.

## Phase and gate accounting

| Phase | Current state | Remaining requirement |
| --- | --- | --- |
| 0: baseline and preparation | Runnable development environment, manifest, registry audit, and local data configuration. | Resolve baseline registry defects before a release freeze; establish and control a new untouched scientific evaluation partition. |
| 1: existing-agent pilots | Recorded condition investigation and multistep route investigation, with original unsuccessful task choice retained. | Independent chemist assessment of usefulness and remaining scientific questions. |
| 2: scientific contracts/access | Local operations and CLI use the existing packages; direct-call parity covered; unresolved-signature assessment semantics corrected with regression coverage. | Independent review of the updated assessment contract; MCP is deferred until a selected client needs it. |
| 3: persistence | Local immutable artifacts, history, replay, saved user conversations, exact runtime thread resumption, cancellable background turns and progress. | Durable queue/crash adoption and remote storage are not implemented. |
| 4: condition investigations | Structured precedent inspection, procedure linkage, missingness, attributed adaptations and recorded custom calculations implemented. | Independent scientific review and source-supported transfer trials; recording a proposal does not validate it. |
| 5: iterative retrosynthesis | Issue-backed planner revision plus external proposal assessment, step inspection, branch edits, full reassessment, declared-material checks and recorded alternative comparison implemented. | Independent chemistry review, broader complex-route trials, actual stock/procedure evidence, and automatic condition-selectivity repair remain. No comparative superiority is claimed. |
| 6: comparative evaluation | Not performed. | Blind comparison, adjudication, matched-budget repeated runs, and untouched evaluation in roadmap order. |
| 7: release/deployment | A local question-to-answer browser preview is usable for development testing. | Production release is not performed: full scientific gates, multiuser isolation, consolidation, and selected remote/MCP deployment remain. |

No new untouched cases were opened or used to tune this implementation. The
authored pilots are development cases, not independent generalization evidence;
overlap with the training/source corpus has not been ruled out.

## Baseline defects and scientific limits

The registry audit reports two `UNDECLARED_AMBIGUOUS_IDENTIFIER` issues involving
`cas:76189-55-4` and `cas:98327-87-8`. The stored records have overlapping BINAP
aliases. Deciding whether to merge identities, correct stereochemical names, or
declare ambiguous aliases requires curation; this implementation does not
silently change those identities or mark the release baseline clean.

`assess_reaction_recipe()` now preserves unavailable signatures as `unknown`,
with no invented hard conflict. Its `compatible=False` value alone still must
not be interpreted as impossibility. Consumers should inspect status, analysis
warnings and unresolved requirements; rule coverage remains incomplete.

The workspace records explicit data files; callers must configure every external
input used by custom scripts. Ordinary execution detects file stat changes;
replay uses full checksums. Fingerprints identify an environment but do not
install old dependencies or restore historical source files automatically.

## Validation

Workspace, route-replay, conversation, structured-answer and condition-investigation
tests are included in the full suite. The 15 conversation tests cover real recorded chemistry behind a model test double,
follow-up recovery, invented citations, cancellation, baseline drift, API origin/token
boundaries, invalid IDs, unowned workers, Windows atomic-write retry, and actual
subprocess success/failure/missing-output/deadline handling, bounded answer correction,
and transient state-reading failures. Model-backed live
smoke tests are recorded separately above and are not counted as deterministic tests.

The selected Compact index passed the canonical integrity validator with
**94,643 rows and no reported integrity issues**. This is artifact integrity,
not independent scientific validation. A fresh workspace instance replayed the
reaction analysis after fully rehashing all five selected data artifacts; the
complete result matched. Ruff checks and documentation link checks passed.

Full suite after the condition-investigation update: **1,619 passed, 2 failed** in
285.20 seconds, including 21 additional regression cases and all seven existing
Markdown/SMILES presentation tests.
Ruff, JavaScript syntax, and local documentation link checks also passed. The failures are the
pre-existing registry checks
`test_validate_reports_current_registry_state` and
`test_registry_audit_reconciles_all_rows`, both caused by the two ambiguous
identifier issues described above. No new test failed.

After the chat interface update: **1,622 passed, 2 failed** in 319.49 seconds,
with the same two registry failures. The chat checks include 11 JavaScript
controller cases for submission, cross-chat cancellation, progress, drafts,
navigation races, and startup recovery, plus API activity/asset regressions.
JavaScript syntax and Ruff checks passed. The restarted local server served the
page, CSS, JavaScript, activity endpoint, and saved Markdown/structure answers
successfully. Browser control remained unavailable, so visual acceptance is still
pending; these checks do not constitute a live browser rendering test.

The dark/light theme follow-up also completed the full suite: **1,622 passed,
2 existing registry failures** in 367.86 seconds. Theme switching, preference
restoration, invalid saved values, and disabled browser storage were checked
with a JavaScript harness. Updated assets were verified on the running server.

The original pilot's local `review_packet.md` links the evidence needed for development
review, and `summary.json` provides its investigation context. That pilot's final
status is `insufficient_evidence`, with missing procedure evidence
and independent review explicitly outstanding. Automated tests validate software
behavior, not the scientific quality of the proposed conditions or routes.
