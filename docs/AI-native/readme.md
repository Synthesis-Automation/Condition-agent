# Scientific workspace and local agent conversations

Status: local development implementation; independent chemistry review remains pending.

[Fragment precedent search](Fragment_Precedent_Search_Design.md) is available as
an optional workspace operation over a prepared local index. It finds product
cores, distinguishes supported construction from retention, and preserves
unresolved evidence. Independent chemistry review remains pending.

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
Reaction schemes appear beneath the concise answer, one SVG per step, with compact
compound names, concise reagent/catalyst and solvent names above the arrow and a
percentage yield below when supplied. Amounts, temperatures, times and workup
instructions remain in the step details. Retrosynthesis plans
are shown in synthetic direction. Open **Step details & evidence** for full names,
condition text and sources. The step heading retains its declared
reported/proposed/computed status; drawings do not validate feasibility.
Scheme previews prefer 60% of the native SVG width and shrink further to fit a
narrow card, keeping the complete reaction visible. **Download SVG** retains the
native vector drawing. Molecule cards and reaction schemes use the shared
`web_consistent` drawing preset (about 30 pixels per bond before preview scaling).
Larger molecules expand the SVG canvas; standalone molecule cards keep their
intrinsic size and provide scrolling. Restart the server after rendering code
changes to clear cached saved-answer presentations.
Alternative routes are separate expandable sections; the first is initially open,
without implying that it is scientifically preferred. Step evidence, molecule
galleries, SMILES, route connections, and uncertainty remain expandable. Missing
conditions or yields are marked as missing in the details, and step cautions remain visible.
Answers render Markdown tables,
headings, emphasis, lists, code, and links. Wide tables scroll horizontally.
Raw HTML and remote Markdown images are disabled. Artifact hashes in prose become
readable citations: known external sources link to the original paper/patent URL,
and local-only results link to saved evidence within the conversation. The Sources
panel also retains access to captured excerpts. Technical fenced code remains literal.
Older answers without structured steps are not reconstructed from prose.
After updating the server code, restart the server and refresh the browser. Saved
answers acquire the new presentation without rerunning the agent; start a new chat
for investigations after the runtime prompt changes its recorded baseline.
Install `requirements-web.txt` when setting up a new
environment (the renderer uses `markdown-it-py`). Select
a saved conversation to continue it, including after a normal server restart.
The square **Stop investigation** button replaces Send while work is active.
It cancels the runtime process tree and preserves partial evidence. Switching
chats or selecting **New chat** keeps the active investigation visible above the
composer, with **Open chat** and **Stop** controls. You can write a draft while
waiting. Reloading the page discovers the active investigation automatically.

Recognized, parseable molecular SMILES in questions and answers also appear as
SVG structure cards under **View structure(s)** or **Structures & scientific
details**, with their original notation and a **Download SVG** link.
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
shows explicitly identified molecules, reaction schemes, route dependencies,
condition details, yields, and inspectable source records. Each molecule, step,
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

The agent writes the complete `scientific_answer.v2` to the current attempt's
fixed `answer-draft.json`, validates that saved draft and records its evidence
self-review. It then returns only this small acknowledgment:

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
[`answer_contracts.py`](../../chem_coworker/scientific_workspace/answer_contracts.py).
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

Two short versioned guides offer questions, examples, pitfalls and stopping
considerations: [conditions](../../chem_coworker/scientific_workspace/guides/conditions.md)
and [retrosynthesis](../../chem_coworker/scientific_workspace/guides/retrosynthesis.md).
The agent may skip a guide or reorder, repeat or replace its suggestions. The
guides do not impose a required workflow. Tool contracts still enforce scientific
inputs, structural checks and evidence requirements.

At investigation creation, the baseline freezes both guides and up to three
applicable lessons per task (`general`, `conditions`, `retrosynthesis`). Selection
prefers exact-task advice, then general advice, and newer records within each
category. It verifies their recorded evidence and removes duplicate advice.
Follow-up turns read this same snapshot; newly published lessons affect only
new investigations. The snapshot includes guide contents, lesson records and a
checksum so an investigation never silently acquires changed guidance.

These are direct workspace methods, separate from chemistry operations:

| Method | Purpose |
| --- | --- |
| `w.task_guide("conditions")` | Read the optional frozen conditions or retrosynthesis guide. |
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

The agent receives a [chemistry investigation guide](Chemistry_Investigation_Guide.md)
that connects exact target identity, primary-source experiments, local graph and
condition checks, and an explicit challenge of the weakest claim. It chooses its
own tool order. Full-text access remains dependent on the source; no subscription
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

Operation summaries now project relevant fields, preserve scientific statuses,
and disclose truncation/collection counts. Complete results remain immutable.
Route summaries show concise per-step status, evidence tier, admission and warnings;
inspect the selected step for its structures and complete gate evidence. Ordinary
short warning lists do not repeat a separate collection-metadata record.
Start with `w.call_summary(event)` after a recorded call. When a decision needs
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

Code/definition changes invalidate an investigation's baseline. Start a **new
investigation** after implementation changes; do not rewrite a saved manifest to
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

## Available operations and ownership

`catalog` prints exact callable signatures. No arbitrary import or code execution
is accepted through the operation dispatcher.

| Operation | Owner and input | Evidence / limitation |
| --- | --- | --- |
| `analyze_reaction` | `reactive_taxonomy`; `reaction_smiles` | Complete analysis, versions, edits, interpretations, ambiguity, warnings. |
| `analyze_molecule` | `reactive_taxonomy`; `smiles` | Graph-derived target audit; reactive sites are hypotheses. |
| `recommend_conditions` | `condition_recommender`; `reaction_smiles`, optional `top_k`, `search_scope` | Canonical shared-core results with compatibility, ranking, and provenance unchanged. Requires `condition_index` and `shared_core_index`. |
| `search_fragment_precedents` | `reactive_taxonomy` graph/evidence rules and `condition_recommender` discovery index; `query`, optional `query_format`, `topology`, `limit`, `timeout_seconds` | Bounded product-fragment discovery, per-embedding changes, exact source/procedure joins, and saved inspection paths. Requires prebuilt `fragment_index`; never expands a route or rebuilds data. |
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
