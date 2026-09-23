# Scientific workspace implementation status

Date: 2026-09-24

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
- Optional runtime: native `codex exec` with structured final responses and exact
  thread resumption. Codex supplies the iterative tool loop and programmable
  environment; the scientific packages do not depend on its model service.
- Answer boundary: strict response schema, citation existence/type checks,
  explicit unreviewed status, and no-local-evidence status where applicable.
  These are traceability checks, not automatic verification of scientific prose.
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
validity still requires the independent review described below. The new browser
runtime has not yet received a live multistep-route trial; the underlying route
and revision workspace operations retain their earlier tests and development pilots.

Browser automation could not connect to the local browser in this session, so
visual acceptance remains pending. HTTP submission, polling, saved-answer retrieval,
and evidence endpoints were exercised against the running application and real
runtime. A browser page is included; this is still a local development preview.

## Development investigations

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

## Phase and gate accounting

| Phase | Current state | Remaining requirement |
| --- | --- | --- |
| 0: baseline and preparation | Runnable development environment, manifest, registry audit, and local data configuration. | Resolve baseline registry defects before a release freeze; establish and control a new untouched scientific evaluation partition. |
| 1: existing-agent pilots | Recorded condition investigation and multistep route investigation, with original unsuccessful task choice retained. | Independent chemist assessment of usefulness and remaining scientific questions. |
| 2: scientific contracts/access | Local operations and CLI use the existing packages; direct-call parity covered. | Correct the unresolved-signature assessment semantics as separately reviewed domain work. MCP is deferred until a selected client needs it. |
| 3: persistence | Local immutable artifacts, history, replay, saved user conversations, exact runtime thread resumption, cancellable background turns and progress. | Durable queue/crash adoption and remote storage are not implemented. |
| 4: condition investigations | Evidence inspection, recipe resolution/assessment, and derived comparison usable. | A scientifically reviewed adaptation capability; the pilot lacked sufficient source detail to justify one. |
| 5: iterative retrosynthesis | Issue-backed alternate disconnection/realization, lineage, and repeated verification usable. | Broader independent chemistry review, richer agent-proposed step/route assessment, and condition-selectivity integration. |
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

`assess_reaction_recipe()` currently returns a hard-conflict field when a verified
signature is unavailable. The adapter returns this contract unchanged and
documents the limitation. This status alone is not evidence of chemical
impossibility. A domain-level semantic correction needs its own regressions and
migration note.

The workspace records explicit data files; callers must configure every external
input used by custom scripts. Ordinary execution detects file stat changes;
replay uses full checksums. Fingerprints identify an environment but do not
install old dependencies or restore historical source files automatically.

## Validation

Workspace, route-replay, and conversation tests: **33 passed** in the full suite.
The 13 conversation tests cover real recorded chemistry behind a model test double,
follow-up recovery, invented citations, cancellation, baseline drift, API origin/token
boundaries, invalid IDs, unowned workers, Windows atomic-write retry, and actual
subprocess success/failure/missing-output/deadline handling. Model-backed live
smoke tests are recorded separately above and are not counted as deterministic tests.

The selected Compact index passed the canonical integrity validator with
**94,643 rows and no reported integrity issues**. This is artifact integrity,
not independent scientific validation. A fresh workspace instance replayed the
reaction analysis after fully rehashing all five selected data artifacts; the
complete result matched. Ruff checks and documentation link checks passed.

Full suite after the presentation update: **1,584 passed, 2 failed** in
273.93 seconds, including all seven new presentation tests. Ruff, JavaScript syntax,
and local documentation link checks also
passed. The failures are the
pre-existing registry checks
`test_validate_reports_current_registry_state` and
`test_registry_audit_reconciles_all_rows`, both caused by the two ambiguous
identifier issues described above. No new test failed.

The original pilot's local `review_packet.md` links the evidence needed for development
review, and `summary.json` provides its investigation context. That pilot's final
status is `insufficient_evidence`, with missing procedure evidence
and independent review explicitly outstanding. Automated tests validate software
behavior, not the scientific quality of the proposed conditions or routes.
