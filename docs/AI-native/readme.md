# Scientific workspace and local agent conversations

Status: local development implementation; independent chemistry review remains pending.

Linear routes in structured answers now use a shared [three-column SVG scheme](../new/linear_route_svg.md),
with optional conditions and no compound captions. Branched or ambiguous routes
retain dependency and individual-step views. The same renderer is available to
local Python scripts and through `python -m visualization.route_cli`. Restart the
workspace server to load the presentation update; saved evidence is unchanged.

The workspace is organized into `core/`, `adapters/`, `runtime/`, `agent_context/`,
`answers/` and `views/`; the public Python exports and root CLI remain the entry
points. See [package organization](Scientific_Workspace_Core.md#package-organization)
for ownership and internal import migration. `agent_context/` contains Python
prompt/context code, `task_playbooks/` contains optional task strategies, and
`agent_instructions/` contains shared agent rules. Restart the server and create a new
investigation after this source reorganization; saved evidence remains readable.
No chemistry definitions, domain schemas or dataset indexes need rebuilding.

The [core boundary](Scientific_Workspace_Core.md) separates execution and evidence
recording from task advice and presentation. New investigations use
`scientific_baseline.v2`; each browser turn saves independently identified
application resources and runtime settings. Scientific code/data changes still
require a new investigation. Task-guide or presentation-instruction edits can enter
the next recorded turn without invalidating science; older v1 baselines keep their
original strict verification. The existing answer schema and browser layout remain
unchanged in this core phase.

New route assessments use lossless compact artifact storage when large diagnostic
sections are present. `scientific_artifact_storage.v1` wraps the call/replay payload
and identifies linked sections by literal JSON paths. Operator template-ID lists,
reaction signatures and per-molecule audits of at least 2 KiB move to immutable
`scientific_evidence_section.v1` artifacts; repeated identical sections share one
file within the investigation. Small records remain plain JSON. Structures, gates,
warnings, uncertainty and selected precedent reactions remain in the main payload.

`store.read_artifact(ref)` reconstructs the exact original scientific result and
verifies both section hashes and the complete expanded-payload hash. Existing
inspection, route revision, answer validation and replay use this expanded view.
`store.read_artifact(ref, expanded=False)` returns the stored compact representation
while still validating its linked evidence. Existing plain artifacts are unchanged
and readable. Copy the investigation's complete `artifacts/` directory when archiving;
a compact root file alone does not contain all its evidence. Scientific citations
continue to reference the parent call, not internal storage sections.

The web artifact endpoint defaults to the stored representation; append
`?expanded=true` for the original full JSON. A marker's `artifact_ref` can be read
through the same conversation-scoped endpoint to inspect just that section (its
data is under `value`). API clients needing the previous expanded response shape
must request `expanded=true` for new compact records. Scientific result and dataset
schemas are unchanged. No dataset rebuild or historical-file migration is needed.
Restart the server and start a new investigation after the code change; existing
saved conversations remain viewable. Compaction reduces repeated storage and response
size; it does not skip scientific computation or integrity verification.

[Fragment precedent search](Fragment_Precedent_Search_Design.md) is available as
an optional workspace operation over a prepared local index. It finds product
cores, distinguishes supported construction from retention, and preserves
unresolved evidence. Independent chemistry review remains pending.

For a target-derived fragment query, supply `target_smiles` to
`search_fragment_precedents`. The canonical matcher checks aromaticity, bonds,
stereochemistry and the requested ring topology against that target before index
access. A mismatch is recorded as an error, never a zero-hit result. Successful
searches retain versioned `target_validation` evidence in the result and summary.
Omitting the target remains supported for standalone fragment searches; it provides
no validation of membership in a particular target. Existing index and query
identities are unchanged; no index rebuild is needed.

Supporting reactions can now be inspected for each proposed route step with
`inspect_step_precedents(source_ref=..., realization_id=...)` for a saved
`disconnect_target` result. Use `step_id` for `assess_route_proposal` or
`revise_route_branch`; use only `source_ref` for `assess_route_step`. The operation
uses the canonical template lookup or the assessment's saved matches. It returns
actual source reaction structures, product comparisons, template edit context,
scoped counts, and available publication, observation and procedure records.
Pages contain up to five records (`limit=3` by default); follow `page.next_offset`.
Disconnection lookup retains at most the canonical top 20 template records;
assessment lookup is limited to its saved matches. Neither is a corpus-wide search.

Attach each inspection artifact to the corresponding answer step using
`precedent_refs: [event.artifact_ref]`. The service checks canonical reactant and
product identity, including specified stereo, before publication. The web UI draws
the same saved source records under **Precedent support**. The first source scheme,
reported conditions/yield and transfer cautions are visible. Copyable source SMILES,
IDs and detailed comparisons are under **Match details & cautions**;
procedures and search scope have separate disclosures. Missing inspection,
no retrieved precedent, and unavailable evidence are separate display states.
The optional field extends `scientific_answer.v2`; new runtime schemas require an
explicit list, while older saved answers load with an empty list. Historical
answers and reviews are not rewritten. For old answers missing links, the display
recovers an exact-structure saved inspection or previews source reactions already
present in a saved assessment. The preview is labelled as not having a detailed
inspection attached; conditions absent from that assessment remain missing. This
read path does not run scientific tools or claim that the agent reviewed the matches.

Finalization and service publication reject omitted supporting-reaction inspections
when saved calls contain evidence for the final reactant/product pair (including
specified stereo). The error gives a matching inspection ref to attach or exact
`inspect_step_precedents` arguments. Inspecting an earlier alternative or calling
`inspect_route_step` does not satisfy this check. A matching empty inspection cannot
hide known nonempty support. Steps without available local support may still be
presented with their literature citations and limitations.

Add `reference_catalog` to the artifact configuration to resolve publication
metadata; the local example pins the existing compact reference catalog. Missing
catalogs remain explicit. Conditions/yields retain observation IDs and source
uncertainty. Template records have no unique observation ID, so associated
experiments are labelled as reaction-ID joins, never merged or transferred to the
proposed step. These are template precedents, not automatic evidence of feasibility;
the agent must inspect substrate, functional-group and stereo differences before
claiming transfer support. Route admission and condition ranking are unchanged.

Use the existing agent's Python and shell access to investigate chemistry through
the canonical packages. The workspace adds saved evidence and replay without an
internal LLM controller. No model API key or MCP server is required for these
operations. An optional conversation adapter now connects the same workspace to
the installed Codex runtime and a local browser interface.

## Ask a question in the browser

From the repository root, using the project's Python environment:

```powershell
python -m app.web_api --scientific-chat --port 8011
```

Open **http://127.0.0.1:8011/scientific**. No frontend build is needed for this
page. The existing reaction workbench remains available at `/` when its frontend
has been built. This flag enables the research profile; the default focused
recommendation deployment does not expose the agent routes.

Try:

> Analyze CCBr.N>>CCN. What structural evidence supports the transformation,
> and what remains uncertain?

Then ask:

> Does that evidence prove experimental feasibility?

For conditions, use the **Explore conditions** example. For retrosynthesis,
use **Plan a synthesis** or supply a target SMILES with your constraints.
The agent chooses operations, examines results, and can run custom analysis
scripts. It can ask for missing structures rather than invent a reaction.

The page has a conversation sidebar and a compact message composer at the bottom.
The composer starts at one line and grows as you type a longer question.
The page opens in dark mode; use **Light mode / Dark mode** in the top-right corner
to switch. Your choice is saved in this browser. Scientific drawings retain a
light canvas so bonds and element colors stay readable in either theme.
Press **Enter** to send and **Shift+Enter** for a new line. While the agent works,
the conversation shows concise agent-written updates about the current check,
findings and next action, alongside its state and elapsed time. The prompt asks
for an opening update and further updates after meaningful findings or failures,
roughly once a minute during longer work when possible. These are public status
messages, not private reasoning or a validated final answer. Updates and recorded
actions appear together in chronological order within the conversation. Tool rows
show the command or script, search query, web-page target or scientific action;
expand a row inline for its command, inputs and recorded failure details. Each
action updates one row as it runs and finishes. The conversation has one main
scroll area, and scrolling back keeps your reading position. Older chats recover
details from saved runtime logs when available, without changing
their evidence. Restart the server and refresh the page after updating this code;
a browser refresh alone cannot update an already-running Python service.
The final answer appears when ready; the UI
does not stream intermediate drafts or estimate a completion percentage.
The completed timeline remains above the final answer after the turn ends.
If no new activity arrives for a minute, the page shows the time since the last
recorded action and reports that the runtime is still running. This status does
not invent an agent update or imply that a scientific check has succeeded.
Timestamped updates, tool lifecycle events, failure diagnostics,
answer-validation errors and turn exceptions remain available in the debug log.
New turns save this append-only log at `turns/<turn-id>/progress.jsonl`;
the full runtime event stream
and process output remain in `runtime.jsonl` and `runtime.stderr.txt` (and
`repair-1/` for a correction attempt). The debug log is an operational record,
separate from checksum-verified scientific evidence; its live download contains
only complete JSONL lines through the
`GET /api/v1/scientific/conversations/{id}/turns/{turn_id}/debug-log` endpoint.
The chat interface has no debug-log download button. Old chats can show recovered
updates from runtime logs but do not acquire a fabricated historical debug log.
Recorded scientific calls, source fetches/extraction, and custom-script outcomes
also enter activity and the debug log, including failures inside a command that
exits successfully. The server collects committed events during runtime polling
and when the turn ends, retaining their event sequence and original artifact
reference without duplicating them or changing scientific results.
Synthesis answers show the route first, one SVG per step, followed by a concise
explanation of the choice, strongest evidence and main uncertainty (normally one
short paragraph of two or three sentences). Detailed procedures and extended analysis are given when
requested, including in a follow-up; the full source evidence remains saved.
Linear synthesis routes use one three-column molecular SVG without compound
captions or numbers, followed by concise per-step rationale and evidence. The
agent supplies complete ordered `routes` and shared intermediate IDs. Source
reaction drawings are expandable; their reported observations, source identity,
and material cautions remain visible. Conditions already on the route arrows are
not repeated in the step body; their basis is visible and full attribution remains
in **Experimental details & sources**. Single-step schemes use compact
compound names, concise reagent/catalyst and solvent names above the arrow and a
percentage yield below when supplied. Amounts, temperatures, times and workup
instructions remain in the step details. Retrosynthesis plans
are shown in synthetic direction. Outside a continuous route, each step card presents its reaction SVG,
attributed **Why this choice** rationale, and **Precedent support**. The first source
experiment is visible with its own scheme, recorded conditions/yield and publication
link; additional precedents are expandable. Source observations stay separate.
Material step cautions, route gaps and open questions remain visible. Open
**Experimental details & sources** for the full proposed recipe and attribution.

For condition recommendations, independent steps with identical explicit reaction
SMILES share one target scheme; the first recipe is visible and **Alternative
conditions** is collapsed. Different transformations and dependent route steps are
never combined. Alternative routes are separate sections, with the first initially
open. Ordering comes from the agent; rendering adds no ranking.

Both target and precedent schemes use the shared `web_consistent` drawing preset
and a 75% display scale. Wide schemes scroll horizontally instead of shrinking
chemical labels to fit. Record IDs, source SMILES, full procedures, comparison data
and search scope live in details. **Notes & sources** retains additional findings,
molecule qualifications and captured excerpts. The molecular structures and route
dependencies remain in the saved answer data.

See [Presentation layer](Scientific_Workspace_Presentation.md) for the ownership map,
answer fields and evidence path. Restart the server after rendering code changes to
clear cached saved-answer presentations.
Answers render Markdown tables,
headings, emphasis, lists, code, and links. Wide tables scroll horizontally.
Raw HTML and remote Markdown images are disabled. Artifact hashes in prose become
readable citations: known external sources link to the original paper/patent URL,
and local-only results link to saved evidence within the conversation. The Sources
panel also retains access to captured excerpts. Technical fenced code remains literal.
Older answers without structured steps are not reconstructed from prose.
After updating the server code, restart the server and refresh the browser. Saved
answers acquire the new presentation without rerunning the agent. For v2 baselines,
changed prompt resources enter the next turn as recorded application context and
start a new agent thread over saved history. Scientific changes and older v1
baselines still require a new investigation.
Install `requirements-web.txt` when setting up a new
environment (the renderer uses `markdown-it-py`). Select
a saved conversation to continue it, including after a normal server restart.
The square **Stop investigation** button replaces Send while work is active.
It cancels the runtime process tree and preserves partial evidence. Switching
chats or selecting **New chat** keeps the active investigation visible above the
composer, with **Open chat** and **Stop** controls. You can write a draft while
waiting. Reloading the page discovers the active investigation automatically.

Recognized, parseable molecular SMILES in questions also appear as SVG structure
cards under **View structure(s)**, with their original notation.
Answers use the explicit reaction schemes without a duplicate molecule gallery.
Inline code, `smiles` code blocks, and recognizable plain-text notation are
supported. The existing RDKit-based visualization package produces the drawings.
Invalid strings remain readable in the message but receive no drawing; molecule
names alone are not converted to guessed structures. Up to 12 distinct structures
are shown per message. This is a presentation of the supplied structures, not
validation of the proposed chemistry.

The presentation is generated from saved messages and cached across polling.
Reloading an existing conversation renders its tables and structures without
rerunning the agent or changing its answer/evidence artifacts.

### Structured scientific answers

New agent turns use `scientific_answer.v2`. Alongside the explanation, the page
shows reaction schemes, condition details, yields, and inspectable source records.
Molecule identities and route dependencies are retained in the saved data without
extra display panels. Each molecule, step,
condition, yield, and standalone claim carries one of these labels:

| Basis | Meaning |
| --- | --- |
| Input | Supplied by the user |
| Reported | The agent attributes the statement to the cited source; not independently verified |
| Computed | A recorded local computation supplies evidence; execution does not establish experimental success |
| Proposed | An agent hypothesis or suggested adaptation |
| Unknown | Information remains unresolved or missing |

Source records point to checksum-verified local artifacts and specific record,
field, or example locators. External references require an HTTP(S) URL and a
saved source-excerpt attachment; the source panel offers both the capture and
the original URL. The service verifies reference existence/type, but does not
automatically prove that the cited text substantiates each claim.

The `computed` label requires a completed recorded workspace call, replay, or
`run_python` execution. Attaching a custom script/output preserves derived analysis, but does not
by itself establish recorded execution for that label. Such analysis can be
discussed in prose with its attachment and provenance limitation.

The agent can use `w.finalize_answer(draft_path, draft, findings=findings)` to
validate, record its explicit evidence self-review and save the current attempt's
fixed `answer-draft.json`. It fills only authoring boilerplate: empty lists, null
yield/source URL, the current schema version, and `needs_user_input=False`.
Scientific basis, structures, source locators, conditions and findings are never
inferred. Reported/computed objects still require valid source references.
The saved contract remains the full `scientific_answer.v2`; manually authored
complete drafts remain accepted under the same validation rules. For example:

```python
# draft_path is the exact path supplied for this runtime attempt.
import json

receipt = w.finalize_answer(draft_path, {
    "answer_markdown": "Please provide the target SMILES so I can draw and assess the route.",
    "needs_user_input": True,
})
print(json.dumps(receipt))
```

Scientific recommendations also supply concise findings for all five review areas.
A clarification such as the example has no invented review or chemistry. The helper
returns only this small acknowledgment for the agent's final message:

```json
{"schema_version":"scientific_answer_handoff.v1","answer_file":"answer-draft.json"}
```

The server reads the saved file and repeats the scientific answer schema,
evidence and baseline checks, including matching any recorded self-review to
the exact final draft. The handoff does not itself contain or approve scientific
claims. This avoids generating the full answer JSON again in the final message;
the public answer and saved scientific answer contract remain unchanged.
Each attempt has its own draft file; the correction attempt uses `repair-1/`.

A missing or malformed handoff/draft, or a schema/evidence failure, receives at
most one correction on the same agent thread. Rejected submissions and both
attempts remain saved. A second invalid submission fails visibly; runtime errors,
cancellation and baseline drift are not retried. Model and reasoning settings
are unchanged; this reduces duplicated output rather than promising a particular
response time.

Search advice remains optional: if a full-target query finds no construction
evidence, consider a deliberate query of the distinctive core with peripheral
substituents removed. Preserve topology, explain the relaxation and skip it when
it cannot change the plan. Similarly, stop repeated title/DOI or access searches
that produce no useful evidence; usually one alternative access path is enough
before switching sources or reporting a gap. Additional disconnections should
resolve a real planning question, rather than regenerate an already sourced step.

Route steps refer to explicit reactant/product molecule IDs and preceding step
IDs. Cycles, missing IDs, disconnected declared dependencies, and omitted route
dependencies are rejected. The page flags a route that has no declared terminal
product matching its target IDs. These are view-consistency checks, not chemical
identity matching, atom-balance validation, or an assessment of feasibility.
Malformed molecular notation stays inspectable with a drawing error; it is not
silently replaced by a guessed structure. Missing conditions and yields are
shown as missing rather than filled from model memory.

For example, ask: “Investigate a synthesis of `O=C(O)c1ccc2nc(-c3cc(Cl)cc(Cl)c3)oc2c1`.
Show the molecules, reaction steps, conditions and source evidence, and distinguish
reported facts from proposals.” A route ending at a methyl ester must disclose
that the free-acid target has not yet been reached.

Older saved answers remain readable with their existing Markdown/SMILES view;
the application does not reconstruct an attributed route from historical prose.
Restart the server and start a **new investigation** after this contract/code
update because existing fixed-baseline workspaces intentionally reject changed
code. Saved answers and original evidence remain readable. No historical manifest
is rewritten.

The shared contract is in
[`answer_contracts.py`](../../chem_coworker/scientific_workspace/answers/answer_contracts.py).
Its Pydantic schema validates the saved scientific answer. The runtime's strict
final-output schema is the small `scientific_answer_handoff.v1` acknowledgment;
it does not replace the answer schema or its validation. Conceptual answers use
empty object arrays. Unsupported steps may remain proposed or unknown, with
limitations; the schema never promotes an agent's structured answer into an
authoritative chemistry record.

Requirements and configuration:

- A current native Codex CLI supporting `exec --json`, `resume`, and
  `--output-schema`, authenticated through `codex login`. The adapter first
  checks PATH, then Windows IDE-bundled binaries when PATH is missing or obsolete.
  A legacy PATH installation does not silently substitute a reduced runtime.
- Use `--codex PATH_TO_EXECUTABLE` or `SCIENTIFIC_CODEX_PATH` for an explicit
  binary. Windows `.cmd`, `.bat`, and `.ps1` shims are not launched through a shell.
- The model defaults to local Codex configuration. Override it using
  `--agent-model MODEL_ID` only with a model available to your account/provider.
  The desktop and a separately installed CLI can have different model support.
  If the provider rejects the inherited model, restart with a current executable
  using `--codex`, or specify a supported model with `--agent-model`.
  Failed turns retain the provider's JSONL error in their public error message
  and `runtime-observations.json`; stderr is used when no JSONL reason exists.
  `runtime-request.json` records the selected executable and version. The
  adapter never retries a provider rejection with a different model.
- `--chat-artifacts FILE` selects the data configuration. The default is
  [artifacts.local.example.json](../../examples/ai_native/artifacts.local.example.json).
  Missing datasets are recorded, and dependent calls fail explicitly.
- `--chat-root DIRECTORY` defaults to `results/ai_native/conversations`.
  `--agent-timeout SECONDS` overrides the chosen profile, excluding baseline preparation.
  This is a per-attempt wall-clock bound, not a token or spending cap. A turn with
  one answer-correction attempt can use up to twice this runtime allowance.

Each question is processed by the actual Codex harness, with Python/shell access
and its available tools; the application does not implement a fixed LLM action
loop. Follow-ups resume the exact stored thread ID, never the globally latest
Codex session. If the adapter's requested runtime settings change, the next turn
starts a fresh thread and reads the same saved investigation evidence. Inherited
Codex configuration is not resolved or fingerprinted, so changes outside this
adapter require a new conversation for a clean comparison.
Questions and tool context can reach the configured model provider;
the browser page is local, but model inference is not necessarily local.
The default model identifier is inherited, not resolved into a pinned model by
this adapter. Runtime version, requested configuration, thread ID, usage, prompts,
event logs, and answers are saved. [Official non-interactive runtime documentation](https://learn.chatgpt.com/docs/non-interactive-mode)
describes the execution and resume protocol.

### Inspect and test a conversation

Under `<chat-root>/<conversation-id>/`, the standard scientific manifest,
`events/`, and `artifacts/` coexist with `conversation.json` and `turns/<turn-id>/`.
Each turn saves its prompt, schema, runtime JSONL, stderr, `answer-draft.json`,
final handoff, progress, and terminal state. The runtime's final message is the
small acknowledgment; the actual answer is read from the saved draft.
Successful answers record hashes of their trace files and
links to checksum-verified evidence. Citation validation checks existence and
artifact type; it does **not** establish that the prose correctly interprets
the evidence. Independent scientific review remains necessary.

To run the real model-backed smoke test against a running server:

```powershell
python -m examples.ai_native.chat_smoke_test --url http://127.0.0.1:8011 --follow-up --output results/ai_native/chat_smoke_report.json
```

This explicitly uses the configured model service. It submits a reaction
question, checks evidence retrieval, and checks that a follow-up uses the same
runtime thread. It is not included in ordinary pytest and is not a chemistry
quality benchmark.

The browser API uses the same `ConversationService` that Python applications can
instantiate. Runtime substitution happens at the `AgentRuntime` protocol;
Codex is currently the only production-provider adapter implemented here.
No MCP dependency is introduced.

### Research profiles and capability checks

Scientific chat now defaults to `--agent-profile research`:

| Profile | Requested reasoning | Requested native web search | Per-attempt deadline |
| --- | --- | --- | --- |
| `research` | high | live | 1800 seconds |
| `quick` | medium | cached | 300 seconds |
| `inherit` | inherited | inherited | 900 seconds |

No profile changes your model. `--agent-model`, `--agent-reasoning-effort`,
`--agent-web-search` and `--agent-timeout` override individual values. Model/provider
support varies; unsupported settings fail visibly rather than silently substituting
a different model. Example:

```powershell
python -m app.web_api --scientific-chat --agent-profile research --agent-timeout 1200 --port 8013
```

Open `http://127.0.0.1:8013/scientific` and start a **new chat**. Choose a free port
if another server is already running. Existing processes do not pick up this update.
`--agent-web-search disabled` disables that native tool only; it does not block
network use through Python or other tools.

Native browser/search access and command networking are separate. This adapter
keeps `workspace-write` and does not enable command network access. Where the
worker cannot download a URL directly, use an available native web reader and
save its actual returned passage with `capture_source`, retaining the failed
download artifact and its limitation. This does not turn an agent-supplied
passage into an independently fetched snapshot. See the runtime's
[network access configuration](https://learn.chatgpt.com/docs/agent-approvals-security).
Permission-denied fetches now carry `network_permission_denied` and explicit
browser-capture recovery guidance. Do not repeat the same blocked download;
retain its failed artifact and report a source-access gap if no reader can open it.

The runtime checks for a working `rg` on PATH or in installed runtime/editor
bundles and adds its directory only to the child process's PATH. It does not
install software or change the user's PATH. The detected path and version are
recorded in runtime metadata; if unavailable, use Python `pathlib`/`re` or
PowerShell `Get-ChildItem`/`Select-String`. Host discovery does not guarantee the
sandbox can execute the binary, so a denied invocation still requires a fallback.

Each attempt saves `runtime-request.json` and `runtime-observations.json`. The
former records requested settings and their origin; the latter counts observed
native search, command and MCP items, including completion/failure signals. An
observed search event does not by itself establish a useful search result. Effective
model/reasoning values remain unconfirmed when the runtime does not report them.
Each turn also saves `capabilities.json`: running interpreter, RDKit parsing probe,
optional PDF parser and configured file presence. File presence does not certify
index compatibility or corpus validation. Inspect the same local checks with:

```powershell
python -m chem_coworker.scientific_workspace capabilities results/ai_native/YOUR_INVESTIGATION
```

### Optional guides and lessons from previous runs

Two short main playbooks describe purpose/context, scientific questions, evidence
distinctions and stopping criteria:
[conditions](../../chem_coworker/scientific_workspace/task_playbooks/conditions.md) and
[retrosynthesis](../../chem_coworker/scientific_workspace/task_playbooks/retrosynthesis.md).
The agent chooses which advice to read and may reorder, repeat or replace suggestions.
Tool contracts own inputs and execution limits; domain validators and the saved-answer
contract own scientific and attribution checks; presentation owns answer formatting.

Specialized advice is available on demand:

| Guide name | Read when useful |
| --- | --- |
| [`conditions_screening`](../../chem_coworker/scientific_workspace/task_playbooks/conditions_screening.md) | Plan a diverse panel with explicitly weak-label evidence. |
| [`retrosynthesis_fragments`](../../chem_coworker/scientific_workspace/task_playbooks/retrosynthesis_fragments.md) | Find construction evidence for a chosen core. |
| [`retrosynthesis_revision`](../../chem_coworker/scientific_workspace/task_playbooks/retrosynthesis_revision.md) | Change a saved branch and compare complete reassessments. |
| [`retrosynthesis_forward`](../../chem_coworker/scientific_workspace/task_playbooks/retrosynthesis_forward.md) | Resolve a consequential product-competition question. |

The existing core discovers `task_playbooks/*.md` and advertises their names. It loads guide
text only on request or explicit inclusion; supporting advice adds no task router,
required sequence or new memory category. See [task ownership](Scientific_Workspace_Core.md#task-playbooks).

At investigation creation, the baseline freezes all discovered guides and up to three
applicable lessons per task (`general`, `conditions`, `retrosynthesis`). Selection
prefers exact-task advice, then general advice, and newer records within each
category. It verifies their recorded evidence and removes duplicate advice.
Procedural lessons stay pinned across follow-ups; newly published lessons affect
only new investigations. Browser turns explicitly snapshot current guides and
instruction resources in an immutable application-context event. Within a turn,
guide reads use that saved text. Headless workspaces retain the initial snapshot
until application context is explicitly recorded. Both snapshots carry checksums.

These are direct workspace methods, separate from chemistry operations:

| Method | Purpose |
| --- | --- |
| `w.task_guide("conditions")` | Read a main or supporting guide by its advertised name from the recorded context. |
| `w.recall_lessons("retrosynthesis", limit=3)` | Read up to three pinned lessons for the selected task. |
| `w.record_lesson(task, advice, applies_when, evidence_refs, scope="code")` | Save concrete procedural advice supported by actual artifact references from this investigation. |
| `w.retire_lesson(lesson_id, reason, evidence_refs)` | Retire recalled advice for future investigations when recorded evidence contradicts it. |
| `w.publish_lessons()` | Publish pending records after a Python/CLI investigation; repeated publication does not duplicate them. |

An agent may record **zero to three lessons per turn**. Useful lessons concern
tool usage, recovery and investigation practices; a successful answer alone is
not supporting evidence. Records retain applicability, source investigation,
supporting artifacts and version scope. The default `code` scope requires the
same scientific code manifest and recorded environment versions. Use
`environment` only for practices independent of code changes; it still requires
matching recorded OS/Python/dependency versions. It does not freeze tool availability;
check current runtime diagnostics before applying advice about an unavailable tool.

The project store is `results/ai_native/lessons.jsonl`. The conversation server
publishes already-recorded lessons and retirements after the turn, including a
failed turn with useful evidence. It does not manufacture lessons from logs.
Python/CLI users call `w.publish_lessons()` explicitly when finished. Keep the
source investigation and artifacts available: unverifiable records are skipped.

Lessons are **agent-authored, unreviewed advice**, not scientific evidence or
automatic chemistry learning. They cannot change graph rules, scoring, admission,
recipes or evidentiary standards; those changes still require the normal code
and chemistry validation process. The feature is development-only. Do not feed
untouched evaluation cases into this learning path; lesson retrieval and changes
are disabled for baselines outside the declared development partition.

### Literature and final evidence review

The [chemistry investigation guide](Chemistry_Investigation_Guide.md) maps the owned
instructions for exact target identity, primary-source experiments, local graph and
condition checks, and an explicit challenge of the weakest claim. The agent chooses
its own tool order. Full-text access remains dependent on the source; no subscription
service, OCR engine or automatic molecular-image extraction is bundled.

The Python workspace provides recorded source tools:

```python
source = w.fetch_source("https://patents.google.com/patent/EP0076530A2/en",
                        title="EP0076530A2")
passage = w.inspect_source(source.artifact_ref, query="LXVIII", limit=4000)
print(passage)  # inspect context and continue reading when needed
if passage["query_found"] and passage["text"].strip():
    excerpt = w.record_source_excerpt(
        source.artifact_ref, start=passage["location"]["start"],
        end=passage["location"]["end"], locator="Captured passage around LXVIII",
    )
```

`fetch_source` preserves original bytes, URL/redirects, retrieval date, extracted
text, checksum, extraction status and explicit failures. Public HTTP(S) downloads
are bounded to 8 MiB; passages are bounded to 16,000 characters. HTML/text and
text-based PDFs are supported (`pypdf` in `requirements-web.txt`). PDF page ranges
are retained; scanned pages, inaccessible documents and missing parsers stay explicit.
Figures and formula drawings require original-source inspection. An exact excerpt
check proves membership in the captured text, not correctness of a chemistry claim.

`w.capture_source(text, url=..., title=..., locator=...)` records text obtained
through another tool as `agent_supplied_excerpt`; it does not claim the workspace
fetched or authenticated that URL. New external answer sources must match their
captured original/final URL and contain actual text. Failed fetch artifacts can
support an access limitation but cannot masquerade as a read source passage.
Legacy source attachments remain supported with their original limitations.

Route assessments accept recorded literature sources and exact excerpts through
`evidence_refs`, retaining them as `literature_provenance` with acquisition and
extraction limitations. Failed downloads and empty source text cannot serve as
route literature provenance. Captured text is not proof of a chemical claim:
source attribution does not override structural, topology or condition gates.

For recommendations the prompt requests `w.record_evidence_review(draft, findings)`
covering source identity, structure/stereochemistry, conditions/yields, route
completeness and counterevidence. Each finding has `area`, `claim`, `assessment`,
`evidence_refs`, and `reason`; use `not_checked`/`not_applicable` explicitly when
needed. `supported`, `partial` and `conflicting` require evidence references.
The record is bound to the exact draft and current turn. Answers record whether
that self-review exists; absence is not disguised as success. It remains an
**agent self-review**, not an independent chemist assessment or semantic verifier.

`w.call_summary(event)` now prints a brief decision view: execution status,
key warnings, candidate or step status, and non-passing gates. It gives an
artifact pointer for inspecting omitted evidence.
The `scientific_call_brief.v2` view reads the saved result directly: counts refer
to the full saved collections, and non-passing gates are selected before any
preview limit. Large route previews prioritize non-actionable or warned steps.
Repeated empty cautions and field catalogues are omitted. Structure strings and
evidence IDs remain usable; exceptionally long strings carry a truncation flag.
Condition summaries include ingredient names, reported operating fields and
source IDs. Procedure summaries show which observations have text, with its
length and saved record index; they do not print entire procedures. Unresolved
recipe identities stay visible even when a component preview is limited.
Omitted quantities/operating fields mean they are not displayed, not zero or
experimentally unnecessary. Inspect the saved recipe/procedure before using it.
`w.call_summary(event, detailed=True)` retains the larger projected summary with truncation and
collection counts when needed; neither view replaces the immutable full result.
Start with the brief view after a recorded call. When a decision needs
more detail, inspect the relevant artifact fields or collection slice instead
of printing the entire result again. For an existing condition-recommendation
call:

```python
print(w.call_summary(event))
print(w.inspect_artifact(
    event.artifact_ref, path=("result", "recommendations"), offset=0, limit=3,
))
# Inspect a selected recommendation only when its full recipe is needed.
print(w.inspect_artifact(
    event.artifact_ref, path=("result", "recommendations", 0, "resolved_recipe"),
))
```

Paths use literal dictionary keys and list indices. Inspection returns a bounded
preview with pagination, explicit truncation and relevant surrounding statuses,
errors and warnings; `limit` is 1–20. Follow the returned path and page information
when details are missing. Ancestor context uses short cautions and explicit pointers
for nested details, with its own budget, so a scalar lookup does not print enclosing
result trees or lose the selected value to surrounding warnings.
The full checksum-verified artifact remains available
through `w.store.read_artifact(ref)`. The prompt supplies task guidance, an operation
overview and targeted signature lookup, so reading the whole README or source files is useful
only when a specific question requires it. Smaller observations do not justify
skipping evidence checks or hiding uncertainty.

Reuse saved audits for unchanged structures. Put executable analysis in
`if __name__ == "__main__":` blocks so importing helper definitions does not rerun
earlier calls, and prefer saved scripts to deeply nested shell quoting. For an
unfamiliar operation, inspect only its catalog entry. `w.store.note` accepts
`hypothesis`, `decision`, `question`, `limitation`, and `review`; use `decision` for
branch selection. After a recorded environment-wide network denial, use available
browser research and capture the inspected passages instead of repeating direct
downloads against different URLs.

Calls record baseline/evidence checking, operation and serialization timings plus
serialized result bytes, so later performance changes can target measured costs.
Failed phases are not reported as successfully timed; artifact persistence is
outside these timings.

Run the separate bounded real-agent research integration check explicitly:

```powershell
python -m examples.ai_native.research_smoke --run-live
```

It uses an isolated conversation and a 300-second outer budget, requests one native
search, source capture/excerpt, recorded target analysis and final self-review.
A failed direct download can be followed by one native source open and honest
capture of returned text; it never replaces missing source text with memory.
Reports under `results/ai_native/research_smoke/reports/` preserve missing or failed
capabilities. This uses your configured model account. It does not measure whether
Codex outperforms ChatGPT. [Development cases](../../examples/ai_native/research_development_cases.json)
provide review criteria for stereo ambiguity, proposed conditions, unsupported
chemistry, literature access and conflicting sources. These are known development
cases, not the untouched chemistry evaluation.

### Local execution limits

This is one trusted user's development environment, with one active turn per
service instance. The launcher accepts only loopback hosts; browser mutations
require the page's session token and pass same-origin checks. Runtime settings
and data paths are configured on the server, not supplied by browser requests.
The runtime explicitly uses `workspace-write`, with no sandbox-bypass flags.
This is not a multiuser isolation boundary or a production hosting setup.

Scientific code/definition changes invalidate an investigation's baseline. Start a
**new investigation** after scientific implementation changes; do not rewrite a saved manifest to
make an old conversation pass validation. The agent is instructed to keep new
analysis files inside its investigation, and source/data baseline checks run
before and after a successful turn. These checks detect drift; they do not
restore files or independently enforce immutability of every local resource.

Graceful shutdown requests cancellation. After an abrupt process crash, the UI
marks an unowned in-progress turn as interrupted and preserves its lock. Inspect
the worker PID in `turn.json` and `.conversation.lock` before removing a stale
lock; do not remove it while the worker is alive. Starting a new investigation
is always preferable when ownership is uncertain. There is no durable queue,
automatic crash adoption, remote storage, or cross-machine Codex thread migration.

API additions, enabled only with the scientific service in the research profile:

| Method | Path | Purpose |
| --- | --- | --- |
| GET | `/api/v1/scientific/config` | Runtime/data availability and page session token |
| GET | `/api/v1/scientific/conversations` | Saved conversations |
| GET | `/api/v1/scientific/activity` | Worker-owned active turn and recent tool events, across chats |
| POST | `/api/v1/scientific/turns` | Submit question, optional conversation ID; returns 202 |
| GET | `/api/v1/scientific/conversations/{id}` | History and current progress |
| POST | `/api/v1/scientific/conversations/{id}/cancel` | Request cancellation |
| GET | `/api/v1/scientific/conversations/{id}/artifacts/{ref}` | Read verified artifact JSON |
| GET | `/api/v1/scientific/conversations/{id}/turns/{turn_id}/debug-log` | Download a snapshot of the turn's timestamped JSONL progress/error log |

## Start an investigation

Run from the repository root in the project's Python environment. Python 3.10+
and RDKit are required. The current environment also has the local Compact
condition index, its shared-core companion, the operator library, and stock
index used by the example configuration. These generated artifacts are not
included in Git; availability on another checkout must be checked.

Copy and adjust [the artifact configuration](../../examples/ai_native/artifacts.local.example.json)
when using other data. Paths in that file are relative to `--repository` unless
absolute. The configuration explicitly selects Compact rather than silently
switching datasets; use Full only by selecting its matching index and companion.

```powershell
python -m chem_coworker.scientific_workspace init results/ai_native/my_investigation --objective "Investigate condition transfer" --artifacts examples/ai_native/artifacts.local.example.json
python -m chem_coworker.scientific_workspace catalog results/ai_native/my_investigation
```

Initialization hashes selected data files, code, and definitions, records versions
and a registry audit, and refuses to overwrite an existing investigation. Large
data files take time to hash. Missing configured inputs are recorded explicitly;
operations requiring them produce recorded errors. Successful initialization
does not establish that all tools or chemistry are validated.

For a structure-only investigation, omit `--artifacts`. Analysis and registry
operations remain usable. Code/definition edits require a new investigation;
continuing against changed inputs must not silently reuse the old baseline.

## Invoke an operation

Save an input file such as `results/ai_native/reaction_request.json`:

```json
{"reaction_smiles": "CCBr.N>>CCN"}
```

```powershell
python -m chem_coworker.scientific_workspace run results/ai_native/my_investigation analyze_reaction --input results/ai_native/reaction_request.json
```

The response includes a `sha256:...` artifact reference. Retrieve the full result
or resume the investigation from another session:

```powershell
python -m chem_coworker.scientific_workspace show results/ai_native/my_investigation sha256:REPLACE_WITH_RETURNED_HASH
python -m chem_coworker.scientific_workspace summary results/ai_native/my_investigation
python -m chem_coworker.scientific_workspace replay results/ai_native/my_investigation sha256:REPLACE_WITH_RETURNED_HASH
```

Replay rehashes selected data and compares the complete scientific result. It
does not replay a model's private reasoning. A changed result is reported as a
mismatch; the earlier evidence is retained. Ordinary calls check code hashes,
runtime versions, and data size/modification identity; full data hashing occurs
at initialization and replay. This is a local reproducibility check, not a
tamper-proof or multi-user security boundary.

## Prepare and search fragment precedents

An optional first step, `suggest_search_fragments`, needs only a target structure,
with no index or retrosynthesis call. It proposes up to five overlapping search
regions (complete ring systems, contextual variants, scaffolds and functional
regions). The agent chooses whether any is useful; it can also supply its own core.

```python
event = w.run("suggest_search_fragments", {
    "target_smiles": "CC(=O)c1ccc2c(c1)COc1ccccc1-2", "limit": 5,
})
print(w.call_summary(event))
# Inspect the selected candidate before making a separate precedent-search call.
print(w.inspect_artifact(event.artifact_ref, ("result", "candidates", 0)))
```

Each candidate includes its query, target atom IDs, omitted atoms, boundary bonds,
structural descriptors and cautions. IDs use the returned canonical target;
`target_atoms[].input_atom_index` links to the original parsed input. To extract
an agent-chosen region, repeat the call with `selected_atom_ids=[...]`. Selections
must be connected and preserve full ring systems and required valence/stereo
context; invalid selections are rejected rather than silently expanded.

Ordering is a transparent structural heuristic, not rarity or synthetic difficulty.
Boundary bonds describe query extraction, not recommended disconnections. Simple
cores carry a broad-query caution; only a subsequent search establishes breadth
in the indexed corpus. Suggestions do not automatically search any candidate.

Build once offline from canonical observations (not the condition-admitted index):

```powershell
python -m condition_recommender.fragment_search build --source datasets/literature/full/combined_records.jsonl.gz --procedure-catalog datasets/literature/full/experimental_detail_catalog.jsonl.gz --output results/ai_native/indexes/fragment_precedents.sqlite
```

The example artifact configuration already names this path as `fragment_index`.
Restart the scientific chat server and begin a new investigation after updating
code or replacing an index; existing investigations keep their original baseline.
An absent index produces an explicit capability error and is never built by an
agent call. `--max-records N` is available for development pilots; results expose
their restricted `prefix_pilot` coverage.

With an open workspace `w`:

```python
event = w.run("search_fragment_precedents", {
    "query": "c1ccc2c(c1)COc1ccccc1-2", "limit": 5,
})
print(w.call_summary(event))
print(w.inspect_artifact(event.artifact_ref, ("result", "hits", 0, "matches")))
```

SMILES queries allow peripheral substitution and preserve the represented ring
systems. For deliberate flexibility, supply `query_format="smarts"` and explicit
constraints; partial ring queries require `topology="subgraph"`. Query a complete
distinctive core rather than an unrestricted common fragment such as biphenyl.

Inspect `search_status`, `source_scope`, count precision, and the returned hit's
literal `inspect_paths`. `complete` covers the indexed scope; `too_broad` requests
refinement; `partial` records a budget limit. Neither missing data nor a timeout
establishes that a core cannot be synthesized. Construction evidence currently
uses validated supplied maps; other correspondence remains unresolved. Conditions
remain reported observations, not recommendations for a new target.

Procedures join by exact observation ID or explicitly unassigned reaction scope.
Long text is a pageable list of chunks (300 characters each), with hashes and
offsets; per-procedure truncation at 60,000 characters is explicit. Saved worker
diagnostics include stage timings, errors, and stderr. The default workspace
deadline is 10 seconds, with a maximum of 30 seconds including worker startup.

## Evaluate agent use of fragment tools

The main retrosynthesis playbook identifies decision-changing gaps; optional
[fragment advice](../../chem_coworker/scientific_workspace/task_playbooks/retrosynthesis_fragments.md)
helps select a core, inspect construction evidence and explain transfer. The agent
can skip this advice or choose its own core; suggested initial searches remain flexible.

For an explicitly live development comparison, use a fresh output directory:

```powershell
python -m examples.ai_native.fragment_agent_comparison --run-live --output results/ai_native/fragment_agent_comparison_new --timeout 240
```

This uses the configured agent model with its existing sandbox, two structure-only
development cases, and the existing artifact configuration in
`examples/ai_native/artifacts.local.example.json`. It selects `fragment_index`,
`retro_library`, `condition_index` and `procedure_catalog` when configured so both
arms can inspect ordinary precedent records and procedures. Supply `--artifacts PATH` to use your own
mapping, or `--case cyclic_ether` for one pair. Model calls can incur usage.

Each case gets fresh agent-only and fragment-assisted threads with the same
scientific baseline and budgets, alternating arm order across cases. The former
is instructed to avoid both fragment tools and their underlying APIs; this is a
prompt-based ablation, not an access-control boundary. Both arms may use other
available research tools. Shared lesson recall/publication is disabled, and agents
are instructed not to inspect other runs or evaluation answers. One runtime
attempt per arm uses the production investigator prompt and answer validators;
the chat service's repair loop is not included.

`comparison.md` and `comparison.json` retain every outcome, latency, recorded calls,
returned construction observations and detected arm violations. Individual runs
keep `answer.md`, `trial_report.json`, workspace evidence and runtime logs.
Construction hits are not counted as useful inspected precedents automatically.
Review the cited source, transfer argument and unsupported steps before filling
the pending manual-review fields. Timeouts are not zero-quality chemistry scores.
These cases are neither independently reviewed nor guaranteed absent from model
training; this comparison does not satisfy untouched-evaluation release gates.

The [initial live validation](../../results/ai_native/fragment_agent_comparison/validation.md)
exercised core selection, broad-query refinement and source/transfer interpretation,
but all four 240-second attempts timed out before final submission. A separate
saved-evidence follow-up completed successfully. Use an explicit larger budget
(for example `--timeout 600`) for a follow-up comparison; faster or better completed
answers have not yet been established.

## Deterministic fragment-guided retrosynthesis POC

The [deterministic development POC](Deterministic_Fragment_Retrosynthesis_POC.md)
connects target-derived fragment queries and observed construction witnesses to
the canonical single-step retrosynthesis engine. It records source-operator
compilation, rejection evidence, witness-directed searches and repeat checks,
without an LLM. Its default example is `Fc(cn1)cc2c1c(c3ccccc3OC)n[nH]2`.

```powershell
python -m examples.ai_native.fragment_guided_retrosynthesis_poc --output results/ai_native/my_fragment_retro_poc
```

This produces local JSON evidence and an HTML molecular review report. It is a
bounded single-step investigation, not a complete route or independent evaluation.

The regular research workbench also exposes **Fragment-guided retro**. Start with
`python -m app.web_api --workbench --build --port 8000`, select the mode and use
the default **Assisted fragment research** workflow: keep a full target, choose
or edit a separate fragment query, search, revise deliberately, inspect procedures,
select source observations and assess transfer. Query revisions, notes, results
and errors remain in bounded browser history and can be exported. The baseline
is optional. Discovery requires the prepared fragment index; transfer also needs
the selected compact/full operator library. **Automatic POC comparison** retains
the earlier automatic selection and baseline experiment. Both workflows share
the domain implementation and chemistry gates. See the
[POC documentation](Deterministic_Fragment_Retrosynthesis_POC.md) for limits and
HTTP contracts. Interactive calls run once; use the CLI for pinned, recorded
investigations and repeated-run checks.

## Focused molecular inspection

Two optional, local RDKit tools help answer structural questions without loading
an index or running forward prediction:

- `compare_molecules`: compare a target and a precedent, returning a strict common
  core, unmatched atoms, attachment boundaries, coverage and possible alignments.
  Supply `core_smiles` to inspect your chosen core without an MCS search. Otherwise
  MCS defaults to a two-second search timeout (maximum five); identical constitution
  skips the search. Timeout results remain partial, and symmetric alignments are
  explicitly ambiguous. Automatic search returns one core, not all equivalent cores.
- `inspect_reactive_sites`: inspect selected atoms and their neighborhood using the
  existing motifs and steric/electronic descriptors. Other detected sites stay
  visible, without claiming experimental competition or selectivity. Omit the
  selection to inspect all sites. Potential unassigned atom and double-bond stereo
  are included alongside specified stereo.

```python
comparison = w.run("compare_molecules", {
    "left_smiles": "Cc1ccccc1", "right_smiles": "Clc1ccccc1",
    "core_smiles": "c1ccccc1",
})
print(w.call_summary(comparison))
# Complete atom pairs and attachment differences remain in the saved result.
print(w.inspect_artifact(comparison.artifact_ref,
                         path=("result", "alignments"), limit=1))

site = w.run("inspect_reactive_sites", {
    "smiles": "CCBr", "selected_atom_ids": [2], "radius": 1,
})
print(w.call_summary(site))
```

Atom IDs are zero-based positions in each **returned canonical SMILES**, consistent
with `suggest_search_fragments`; `input_atom_index` links to the parsed original
SMILES. Inspect `result.molecule.atoms` (or `result.left.atoms` / `result.right.atoms`)
before selecting unfamiliar structures. Salt mixtures, radicals and wildcards are
rejected explicitly; nothing is silently stripped, neutralized or tautomerized.
Charge, isotopes, aromaticity and ring topology remain constrained. Original SMILES
and warnings are saved, including ignored atom-map labels.

These comparisons are structural evidence, not reaction atom mapping or proof that
conditions transfer. `core_smarts` is a search pattern; strict checks apply to the
returned atom pairs. Stereo equivalence is assessed only for otherwise identical
graphs; differing graphs retain separate stereo inventories and an explicit
unassessed comparison. A changed CIP label is not interpreted as mechanistic
inversion. Full descriptors, evidence and version metadata remain in call artifacts;
the agent sees concise summaries first. These tools do not satisfy chemistry-review
or untouched-evaluation release gates.

## Assess starting materials

`assess_starting_material` helps the agent decide whether to expand a route leaf.
It tries exact identity in `condition_registry`, exact product-component identity
in the optional `fragment_index`, then a molecular-weight fallback. It is a
planning policy, not a confirmation of stock or synthetic accessibility.

```python
event = w.run("assess_starting_material", {
    "smiles": "CCO",
    "mw_threshold": 200.0,
    "allow_registry_stop": True,
    "allow_literature_stop": True,
    "allow_mw_stop": True,
    "unavailable_starting_materials": [],
})
print(w.call_summary(event))
```

The threshold defaults to 200 g/mol from the validated
`condition_recommender/definitions/starting_material_policy.v1.json`; equality
does not pass. MW is RDKit `Descriptors.MolWt` for the complete supplied form.
Disabled stages are skipped, and later stages are skipped after a stopping
decision. Explicitly unavailable structures override every stopping rule.
Invalid structures, query atoms and radicals are rejected. Atom-map labels are
ignored; specified/unspecified stereo, isotopes, charges, salts and tautomers
retain distinct identities. No component is silently removed or neutralized.

The exact product lookup uses the existing SQLite unique SMILES index and bounded
observation joins. It does not load the substructure library and needs no index
rebuild. Product-component occurrence alone does not establish an isolated
preparation; returned observation/reaction/reference IDs support further inspection.
A disconnected salt query is not matched to only its organic product component.
Missing, incompatible or failed indexes are distinct from a completed no-match
lookup and remain visible even when the MW fallback permits stopping.

Results use `starting_material_assessment.v1` with a versioned effective policy,
descriptor provenance, evidence and explicit `stop_expansion` / `stop_reason`.
All accepted leaves have `status="assumed_terminal"` and `availability="unknown"`.
Registry stops assume obtainability from curated membership; literature and MW
stops leave a route partially resolved pending preparation/supply evidence. The
agent must retain those assumptions in its answer. Route admission is unchanged.
Restart the workspace server and start a new investigation after installing the
code; historical investigations and index artifacts are unchanged. This is a
development capability, not an independent chemistry-review release gate.

## Available operations and ownership

`catalog` prints exact callable signatures. No arbitrary import or code execution
is accepted through the operation dispatcher.

| Operation | Owner and input | Evidence / limitation |
| --- | --- | --- |
| `analyze_reaction` | `reactive_taxonomy`; `reaction_smiles` | Complete analysis, versions, edits, interpretations, ambiguity, warnings. |
| `analyze_molecule` | `reactive_taxonomy`; `smiles` | Graph-derived target audit; reactive sites are hypotheses. |
| `compare_molecules` | `reactive_taxonomy`; `left_smiles`, `right_smiles`, optional `core_smiles`, `timeout_seconds` | Strict structural cores, differences, coverage, ambiguous alignments and stereo scope. No dataset dependency or reaction-map claim. |
| `inspect_reactive_sites` | `reactive_taxonomy`; `smiles`, optional canonical `selected_atom_ids`, `radius` (0–3) | Focused motifs and existing descriptors, other sites and assigned/unassigned stereo. No experimental selectivity prediction. |
| `assess_starting_material` | `condition_recommender` composes registry identity and taxonomy descriptors; `smiles`, optional `mw_threshold`, three `allow_*_stop` flags, `unavailable_starting_materials` | Ordered registry / exact product / MW stopping policy with evidence and warnings. Optional `fragment_index`; no substructure scan or verified availability claim. |
| `recommend_conditions` | `condition_recommender`; `reaction_smiles`, optional `top_k`, `search_scope` | Canonical shared-core results with compatibility, ranking, and provenance unchanged. Requires `condition_index` and `shared_core_index`. |
| `generate_weak_label_screening_array` | `condition_recommender`; `reaction_smiles`, optional `array_size` (default 24, currently 1–250), `source_reaction_type_hint` | Diverse intact recipes from separate weak-label observations after graph-query and compatibility checks. Requires `weak_label_records`; its sibling recipe catalog is automatically baseline-pinned as `weak_label_recipe_catalog`. No structural condition index required. Preserves unverified-source warnings, recipe IDs and source row numbers; may return fewer recipes. |
| `search_fragment_precedents` | `reactive_taxonomy` graph/evidence rules and `condition_recommender` discovery index; `query`, optional `query_format`, `topology`, `limit`, `timeout_seconds`, `target_smiles` | Bounded product-fragment discovery, per-embedding changes, exact source/procedure joins, and saved inspection paths. Supply `target_smiles` for target-derived queries to validate their semantics before scanning. Requires prebuilt `fragment_index`; never expands a route or rebuilds data. |
| `suggest_search_fragments` | `reactive_taxonomy.search_fragments`; `target_smiles`, optional `limit` (1–5), `selected_atom_ids` | Optional, overlapping, target-derived queries with atom provenance and boundaries. No index, corpus call, retro, mapping, forward check or mandatory workflow. |
| `get_precedents` | Canonical index; `reaction_ids`, optional `offset`, `limit` | All indexed fields, distinct observation IDs, admission/condition status, missing IDs, pagination. Indexed records are reduced representations of source data. |
| `get_procedures` | Configured `procedure_catalog`; `reaction_ids` | All matching procedure observations, including missing fields. No invented procedure text. |
| `inspect_condition_precedents` | Canonical condition index, optional procedure catalog; `reaction_smiles`, `reaction_ids`, optional `offset`, `limit` | Selected-observation structural differences, full recipes, compatibility, publication counts, missing fields and exact procedure links. No new ranking or transfer claim. |
| `propose_condition_adaptation` | Registry and canonical recipe assessment; `source_ref`, `observation_id`, `components`, `operating_conditions`, `change_reasons`, `evidence_refs`, `assumptions`, `risks` | Preserve the inspected recipe and an explicitly proposed replacement, attributed changes and compatibility. See the example below. |
| `resolve_recipe` | `condition_registry`; typed `components` and optional operating values | Canonical identities, contextual roles, raw identifiers, uncertainty, provenance. |
| `assess_recipe` | `condition_recommender`; `reaction_smiles`, resolved `recipe` | Existing compatibility result; no yield or experimental-success prediction. |
| `disconnect_target` | Existing single-step coworker and `core_retrosynthesis`; `target_smiles`, optional search limits and `include_conditions` | One target only: validated strategies, concrete precursor realizations, precedent IDs and warnings. Requires `retro_library`; conditions default off and require condition artifacts when enabled. No stock index, recursive expansion or internal LLM review. |
| `assess_route_step` | Canonical external-proposal assessment; `proposal`, optional `include_conditions`, `evidence_refs` | Structural, operator, precedent and compatibility gates; optional forward challenge is separate. Requires `retro_library`; no stock index required. Supplied resolved recipes are assessed separately from retrieved conditions. |
| `assess_route_proposal` | Same assessor plus route topology; `proposal`, optional `unavailable_starting_materials` and assessment options | Retains invalid/unsupported proposals for inspection; declared material constraints are checked against graph-matched leaves. |
| `inspect_route_step` | Saved proposal `source_ref`, `step_id` | Step gates, graph-matched upstream/downstream steps, molecular audits and supplied-recipe assessment. |
| `inspect_step_precedents` | Saved `source_ref`, optional `realization_id` or `step_id`, `offset`, `limit` | Actual supporting reactions for a disconnection realization or assessed step. Use `realization_id` for disconnections, `step_id` for route assessments/revisions, and neither for a single-step assessment. Follow `page.next_offset`; the shared answer contract validates inspection links for final step structures. |
| `assess_route_step_forward` | Saved route `source_ref`, eligible `step_id`, decision-changing `question`, optional `timeout_seconds` (1–30; default 30) | Optional single-step product-competition challenge using a prebuilt, baseline-pinned `forward_library`. A killable worker records stages and timeout/error; no route-wide prediction or library rebuilding. Does not upgrade the saved route's admission. |
| `revise_route_branch` | Core explicit route edit plus complete reassessment; `source_ref`, `remove_step_ids`, `replacement_steps`, `reason`, `risks`, optional `assumptions`, `evidence_refs` | Preserves the source, inherits material constraints and condition settings, reassesses all steps and topology. Optional forward challenges are separate and not inherited. Empty removals can extend a leaf branch. No automatic improvement or admission claim. |
| `compare_route_proposals` | Recorded results; 2–5 `source_refs` | Same-target comparison with identical settings/constraints; separate gates and missing evidence, no synthetic route score. |

Recipe input example:

```json
{
  "components": [
    {"raw_identifier": "ethanol", "identifier_type": "name", "source_field": "user"}
  ]
}
```

Do not fill unreported operating values with defaults. `source_field` preserves
provenance; the registry decides identity and roles. Proposed recipes remain
hypotheses even after normalization or compatibility assessment.

`assess_recipe` distinguishes `conflict`, `unknown` and `invalid_input`. Inspect
`status`, `hard_conflicts`, `warnings` and `unresolved_requirements`;
`compatible=False` alone does not establish chemical incompatibility. Precedent
inspection counts describe the selected indexed scope, not independent publications.
Procedure links distinguish exact observation IDs from unassigned reaction-level records.

Retrosynthesis in the scientific workspace uses **single-step calls only**. The
agent owns multi-step planning: it selects a concrete realization, chooses the next
intermediate, records alternatives and branch links, avoids cycles, and decides
when to stop or broaden the investigation. It must not invoke the built-in multistep
planner through custom scripts. The standalone planner remains available outside
this workspace; `plan_routes`, `revise_routes`, and `prepare_route_proposal` are no
longer workspace operations. Historical artifacts remain readable, but retired
planner calls cannot be replayed through the current operation catalog.

```python
first = w.run("disconnect_target", {
    "target_smiles": "CC(=O)Nc1ccc(-c2ccccc2)cc1",
    "top_k": 3, "max_templates_to_apply": 40, "max_candidates_to_validate": 10,
})
print(w.call_summary(first))
strategies = w.store.read_artifact(first.artifact_ref)["result"]["strategies"]
# Inspect representatives and alternate_realizations, then explicitly choose an
# intermediate for another disconnect_target call. No branch is expanded for you.
```

Record selected strategy/realization IDs, evidence and stopping reasons in notes.
Assemble chosen steps into `assess_route_proposal`, citing the single-step call
artifacts in `evidence_refs`. Use `revise_route_branch` for explicit branch edits
and `compare_route_proposals` for alternatives. Terminal-material availability
requires separate evidence; no candidates means only that this bounded local
search found no supported disconnection. It does not prove synthetic impossibility.

### Focus a single-step search on one bond

`disconnect_target` contract v2 adds an optional required disconnection bond.
It uses the existing operator library and chemistry validation. First inspect
the target using `inspect_reactive_sites` or `compare_molecules`; use the returned
canonical SMILES and its zero-based atom IDs. Atom-map labels and input SMILES
positions are not these IDs. `focus_target_smiles` binds the selection to the
inspected target and is required when supplying a bond. A stale target, invalid
atom IDs, or a nonexistent bond produces an explicit invalid request result.

```python
# For canonical CCNC, IDs 2 and 3 are the N-methyl bond.
event = w.run("disconnect_target", {
    "target_smiles": "CCNC",
    "required_disconnection_bond": [2, 3],
    "focus_target_smiles": "CCNC",
    "top_k": 3,
    "max_templates_to_apply": 40,
    "max_candidates_to_validate": 10,
})
print(w.call_summary(event))
```

The selected bond must be absent in precursors and present as a validated formed
bond in the forward reaction. Both selected atoms must have precursor contributors.
Ring closure can qualify without splitting the molecule; bond-order changes alone
do not qualify. The policy is `disconnection_bond_focus.v1@1.0`, with target and
check contracts `disconnection_bond_focus.v1` / `disconnection_bond_check.v1`.
Final ambiguous/conflicting mapping evidence cannot satisfy the constraint.

Focused generation preserves RDChiral mapped outcomes until the constraint has
been checked, including symmetric outcomes that otherwise share precursor SMILES.
The check runs before the validation-budget cutoff and is confirmed against the
final reaction observation after existing forward/signature checks. Every returned
focused candidate carries `bond_focus_check.status="verified"` and a typed formed
edit witness. No family name is required. This is structural evidence, not proof
of experimental feasibility. Existing selectivity and compatibility cautions remain.

`bond_focus` echoes the target and selected bond. `search_diagnostics` records
per-level template/validation counts, focus rejections, unresolved checks, budget
exclusions and bounded rejection examples. A zero-result search retains its focus
and diagnostics; it is not silently relaxed. Template retrieval is unchanged, so
the desired operator can still fall outside the template budget. Limits remain
per specificity level; focused calls can try more levels than unrestricted calls.
Ordinary calls retain their generation and ranking behavior. No preserved-core
constraints, automatic precursor expansion or automatic focus selection are added.

For a recorded development comparison against an existing library:

```powershell
python -m examples.ai_native.focused_retrosynthesis_pilot --output results/ai_native/focused_retro_dev --library results/operator_retrosynthesis_poc/full_scale_v3/compact/operator_library_v3.json.gz
```

The ten authored probes compare a fixed single template pool under identical total
template and validation limits, independently check selected-bond compliance, and
report target-specific recovery, candidate counts and observational timings.
The runner uses recorded workspace Python execution; reports and corpus-derived
structures stay under `results/ai_native/`. It does not measure autonomous-agent
route quality or satisfy independent chemistry-review / untouched-evaluation gates.
Restart the server and start a new investigation; existing libraries need no rebuild.

The authored targets and single-step `arguments` are in
[development_cases.json](../../examples/ai_native/development_cases.json).
`start_pilots.py` records the first disconnection; subsequent planning is left to
the agent. These development examples do not satisfy independent review gates.

## Investigate an agent-proposed route

Start a **new chat** after upgrading this scientific code. Existing answers and
artifacts remain readable, but their fixed scientific baseline cannot be resumed
against changed code. Example question:

> My target is `O=C(NCC)c1ccccc1`. My initial sketch uses benzoic acid and
> ethylamine for the final amide step, but ethylamine cannot be purchased.
> Investigate making that intermediate, assess the revised route, and compare it
> with the original. Use local evidence and keep selectivity and condition gaps explicit.
> Start with the default structural checks; leave condition retrieval and forward
> checking for a follow-up and report those checks as not run.

The agent now has the proposal, inspection, branch revision and comparison operations
above. Route drawings and the before/after explanation use the existing structured
answer format. A gate may be `unresolved`, `not_run`, or `out_of_scope`; these do not
mean either experimental success or impossibility. The tool's completed execution
does not convert a proposed step into a reported experiment.

For direct Python access in an initialized workspace:

```python
original = workspace.run("assess_route_proposal", {
    "proposal": {"target_smiles": "O=C(NCC)c1ccccc1", "steps": [{
        "external_step_id": "amide", "target_smiles": "O=C(NCC)c1ccccc1",
        "precursor_smiles": "O=C(O)c1ccccc1.NCC",
    }]},
    "unavailable_starting_materials": ["NCC"],
})
revised = workspace.run("revise_route_branch", {
    "source_ref": original.artifact_ref,
    "remove_step_ids": [],
    "replacement_steps": [{
        "external_step_id": "amine", "target_smiles": "NCC",
        "precursor_smiles": "CC=O.N",
    }],
    "reason": "Investigate making the unavailable intermediate",
    "risks": ["Hypothesis: selectivity, actual stock and operating conditions remain unverified"],
})
comparison = workspace.run("compare_route_proposals", {
    "source_refs": [original.artifact_ref, revised.artifact_ref],
})
print(workspace.call_summary(comparison))
```

Each proposed step requires `external_step_id`, `target_smiles` and dot-separated
`precursor_smiles`. Optional supplied mapping is independently validated. Optional
`proposed_conditions` is a resolved recipe from `resolve_recipe`; assessment of that
recipe remains separate from analogue-condition retrieval (`include_conditions`).
Condition retrieval defaults to false and is inherited during revision. To replace an existing step ID, explicitly
include it in `remove_step_ids`. Disconnected remnants are reported as invalid;
they are never silently deleted. Unchanged downstream steps are reassessed as well.

Start with the standard structural assessment, then investigate the gap most likely
to change the answer. In particular, missing halogen or oxygen contributors require
source inspection and accurate reactant structures, not broad forward prediction.
Never invent a donor or mapping merely to pass a gate.

Forward prediction is optional and applies to one eligible step of a saved route
assessment/revision. Use it when competing products could change a route decision:

```python
challenge = workspace.run("assess_route_step_forward", {
    "source_ref": original.artifact_ref,
    "step_id": "amide",
    "question": "Could a competing product make this step unsuitable for the proposed route?",
    "timeout_seconds": 30,
})
print(workspace.call_summary(challenge))
```

The optional `forward_library` artifact is configured in the same artifact JSON as
`retro_library` and pinned at investigation creation. The example configuration
selects `results/operator_retrosynthesis_poc/full_scale_v3/compact/forward_operator_library_v1.json.gz`.
Build the library separately before starting an investigation; scientific chats do
not build libraries or run forward prediction across the entire route. Missing
forward data does not prevent standard route assessment or single-step retrosynthesis.
The worker's 30-second maximum covers loading and prediction, records stage timings,
and stops unfinished work automatically. Manual process polling is unnecessary.
Inspect a timeout or error as an unresolved check; do not report it as passed.
The challenge is separate evidence and does not rewrite the original route assessment.
Structural validation inside `disconnect_target` and the normal route assessor is
unchanged. Historical route-wide forward calls remain readable, but the workspace
no longer permits enabling that expensive path through `include_forward=True`.

The declared material constraint checks route leaves only. Making an intermediate
can remove the need to purchase it, but does not prove that its new precursors are
available. Actual stock lookup and experimental validation remain separate.

Reproducible authored development example, including fresh-session replay:

```powershell
python -m examples.ai_native.route_revision_pilot --output results/ai_native/route_revision_01 --library results/operator_retrosynthesis_poc/full_scale_v3/compact/operator_library_v3.json.gz
```

This script supplies a predefined hypothesis; it is a contract demonstration, not
an agent-quality comparison or untouched chemistry evaluation.

## Investigate and propose condition changes

```python
inspection = workspace.run("inspect_condition_precedents", {
    "reaction_smiles": reaction_smiles,
    "reaction_ids": selected_precedent_ids, "offset": 0, "limit": 10,
})
inspection_ref = inspection.artifact_ref
result = workspace.store.read_artifact(inspection_ref)["result"]
print(result["distinct_reference_count"], result["page"]["next_offset"])
```

Follow `next_offset` when needed. Distinct-reference counts describe the selected
page, not whole-corpus support or independently replicated experiments. Full
indexed recipes and raw procedure records retain their provenance. Procedures
without an observation ID remain reaction-level evidence; do not transfer their
conditions to every observation. Missing data remains missing.

Only after inspecting actual supporting evidence, prepare the complete replacement
component list and operating values. In the following example, the proposal and
rationale variables must come from that investigation; no default temperature,
quantity, time or yield is supplied:

```python
proposal = workspace.run("propose_condition_adaptation", {
    "source_ref": inspection_ref, "observation_id": selected_observation_id,
    "components": proposed_components,  # ConditionComponentInput dictionaries
    "operating_conditions": proposed_operating_conditions,
    "change_reasons": reasons_by_changed_field,
    "evidence_refs": supporting_artifact_refs,
    "assumptions": explicit_assumptions, "risks": unresolved_risks,
})
```

This operation preserves the original recipe, normalizes the replacement through
the registry, records before/after values for every changed bucket or operating
field, and assesses compatibility. Reasons must cover exactly those changed
fields. The proposal is always agent-authored and unreviewed; reference checks do
not establish that an adaptation is justified. Staged-protocol/declared-absence
editing is not supported yet. Abstain when evidence cannot support a change.

## Custom analysis and persistent notes

For **recorded execution**, save a Python script inside the investigation. It
receives input/output JSON paths in `sys.argv[1:3]`. For example, `count.py`:

```python
import json
from pathlib import Path
import sys

request = json.loads(Path(sys.argv[1]).read_text("utf-8"))
inspection = request["evidence"][request["parameters"]["inspection_ref"]]["result"]
count = sum(not row["missing_operating_fields"] for row in inspection["precedents"])
Path(sys.argv[2]).write_text(json.dumps({"observations_with_all_four_operating_fields": count}), "utf-8")
```

Run it through the workspace:

```python
event = workspace.run_python(
    "count.py", {"inspection_ref": inspection_ref},
    evidence_refs=(inspection_ref,), timeout_seconds=60,
)
print(workspace.call_summary(event))
```

The CLI equivalent is `run-python WORKSPACE count.py --input parameters.json
--evidence sha256:REFERENCE --timeout 60`. Parameters are a JSON object. The
runner snapshots inputs and script text, records output/log hashes, rejects
nonzero exits and invalid output, and limits execution to 1–120 seconds. Custom
scripts are trusted local code with your OS permissions, not a security sandbox.
Keep external dependencies in the baseline and do not mutate source/data.
Custom execution is not automatically replayed or chemically validated.

Use Python directly when a question requires a new comparison. Call the same
domain packages and save the script, inputs, outputs, and assumptions. The
workspace does not run attached scripts automatically.

```python
from chem_coworker.scientific_workspace import ScientificWorkspace

workspace = ScientificWorkspace("results/ai_native/my_investigation")
event = workspace.run("analyze_molecule", {"smiles": "CCO"})
workspace.store.note(
    "hypothesis", "The proposed transfer needs additional evidence.",
    evidence_refs=(event.artifact_ref,),
)
workspace.store.attach_file(
    "results/ai_native/my_analysis.py",
    description="Exploratory comparison; assumptions documented in the script",
    evidence_refs=(event.artifact_ref,),
)
workspace.store.set_status("insufficient_evidence", "Need an independent procedure")
```

The CLI also provides `note`, `attach`, and `status`. Notes are explicitly
agent-authored and unreviewed. An agent-written `review` note is not independent
chemist sign-off. Lifecycle states include `active`, `completed`,
`insufficient_evidence`, `failed`, `cancelled`, and `budget_exhausted`. Resume a
stopped investigation with an explicit `active` transition and reason.

`execution_status=completed` means an operation returned normally. Inspect the
domain's `valid`, `status`, warnings, and limitations separately. A normal return
with insufficient scientific evidence is not an execution failure or chemical
success. Exceptions are saved as errors. Ctrl+C saves a cancelled call when
handled by the running process; a forcibly killed process cannot guarantee that.

Storage is a local single-writer directory with atomic JSON writes and immutable
content-addressed artifacts. Concurrent writers fail explicitly. After a crash,
inspect a leftover `.writer.lock` before manually removing it. No durable
background-job queue or remote multi-user execution service is implemented yet.

## Development pilots

```powershell
python -m examples.ai_native.start_pilots --output results/ai_native/new_pilot
python -m examples.ai_native.inspect_condition_transfer results/ai_native/new_pilot sha256:REPLACE_WITH_RECOMMENDATION_HASH
```

The starter records initial calculations. The external agent inspects the
evidence and decides subsequent calls; the script is not a replacement agent
runtime. These are authored development tasks, not untouched evaluation cases.
No new held-out partition or independent review is claimed.

Current baseline limitations include the two ambiguous BINAP identifiers in the
registry audit and shared-core results awaiting independent review. The direct
recipe assessor now reports unresolved signatures as `unknown` and invalid
input separately; `compatible=False` alone does not mean a chemical conflict.
The workspace preserves these contracts and records the limitations. It does
not treat incomplete validator coverage as proof that a proposed reaction is
impossible. Scientific fixes must be separately versioned and reviewed.

See [implementation status](Scientific_Workspace_Implementation_Status.md) for
the pilot findings, validation, and remaining release gates.
