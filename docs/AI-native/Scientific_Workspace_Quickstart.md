# Scientific workspace and local agent conversations

Status: local development implementation; independent chemistry review remains pending.

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

For conditions, use the **Compare conditions** example. For retrosynthesis,
use **Investigate a route** or supply a target SMILES with your constraints.
The agent chooses operations, examines results, and can run custom analysis
scripts. It can ask for missing structures rather than invent a reaction.

The page shows preparation/running state, recent tool activity, the final answer,
uncertainties, and links to full evidence JSON. Answers render Markdown tables,
headings, emphasis, lists, code, and links. Wide tables scroll horizontally.
Raw HTML and remote Markdown images are disabled; evidence links stay scoped
to their conversation. Install `requirements-web.txt` when setting up a new
environment (the renderer uses `markdown-it-py`). Select
a saved conversation to continue it, including after a normal server restart.
**Stop** cancels the runtime process tree and preserves partial evidence.

Recognized, parseable molecular SMILES in questions and answers also appear as
SVG structure cards, with their original notation and a **Download SVG** link.
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

The `computed` label currently requires a completed recorded workspace call or
replay. Attaching a custom script/output preserves derived analysis, but does not
by itself establish recorded execution for that label. Such analysis can be
discussed in prose with its attachment and provenance limitation. The agent is
instructed to validate its draft before submission; the service repeats the
checks. A rejected answer remains a failed turn with saved evidence, and a
follow-up can correct it. There is no automatic server-side repair loop yet.

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
Start a **new investigation** after this contract/code update because existing
fixed-baseline workspaces intentionally reject changed code. Original evidence
remains available for inspection. No historical manifest is rewritten.

The shared contract is in
[`answer_contracts.py`](../../chem_coworker/scientific_workspace/answer_contracts.py).
Its Pydantic schema also supplies the runtime's strict output schema, avoiding a
second definition. Conceptual answers use empty object arrays. Unsupported steps
may remain proposed or unknown, with limitations; the schema never promotes an
agent's structured answer into an authoritative chemistry record.

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
  `--agent-timeout SECONDS` defaults to 900, excluding baseline preparation.
  This is a wall-clock bound, not a token or spending cap.

Each question is processed by the actual Codex harness, with Python/shell access
and its available tools; the application does not implement a fixed LLM action
loop. Follow-ups resume the exact stored thread ID, never the globally latest
Codex session. Questions and tool context can reach the configured model provider;
the browser page is local, but model inference is not necessarily local.
The default model identifier is inherited, not resolved into a pinned model by
this adapter. Runtime version, requested configuration, thread ID, usage, prompts,
event logs, and answers are saved. [Official non-interactive runtime documentation](https://learn.chatgpt.com/docs/non-interactive-mode)
describes the execution and resume protocol.

### Inspect and test a conversation

Under `<chat-root>/<conversation-id>/`, the standard scientific manifest,
`events/`, and `artifacts/` coexist with `conversation.json` and `turns/<turn-id>/`.
Each turn saves its prompt, schema, runtime JSONL, stderr, final response, progress,
and terminal state. Successful answers record hashes of their trace files and
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
| POST | `/api/v1/scientific/turns` | Submit question, optional conversation ID; returns 202 |
| GET | `/api/v1/scientific/conversations/{id}` | History and current progress |
| POST | `/api/v1/scientific/conversations/{id}/cancel` | Request cancellation |
| GET | `/api/v1/scientific/conversations/{id}/artifacts/{ref}` | Read verified artifact JSON |

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

## Available operations and ownership

`catalog` prints exact callable signatures. No arbitrary import or code execution
is accepted through the operation dispatcher.

| Operation | Owner and input | Evidence / limitation |
| --- | --- | --- |
| `analyze_reaction` | `reactive_taxonomy`; `reaction_smiles` | Complete analysis, versions, edits, interpretations, ambiguity, warnings. |
| `analyze_molecule` | `reactive_taxonomy`; `smiles` | Graph-derived target audit; reactive sites are hypotheses. |
| `recommend_conditions` | `condition_recommender`; `reaction_smiles`, optional `top_k`, `search_scope` | Canonical shared-core results with compatibility, ranking, and provenance unchanged. Requires `condition_index` and `shared_core_index`. |
| `get_precedents` | Canonical index; `reaction_ids`, optional `offset`, `limit` | All indexed fields, distinct observation IDs, admission/condition status, missing IDs, pagination. Indexed records are reduced representations of source data. |
| `get_procedures` | Configured `procedure_catalog`; `reaction_ids` | All matching procedure observations, including missing fields. No invented procedure text. |
| `resolve_recipe` | `condition_registry`; typed `components` and optional operating values | Canonical identities, contextual roles, raw identifiers, uncertainty, provenance. |
| `assess_recipe` | `condition_recommender`; `reaction_smiles`, resolved `recipe` | Existing compatibility result; no yield or experimental-success prediction. |
| `plan_routes` | Existing multistep coworker and `core_retrosynthesis`; `settings` | Bounded route alternatives, whole-route checks, issues and repair proposals. Requires `retro_library`, `stock_index`, and condition artifacts when conditions are enabled. |
| `revise_routes` | `core_retrosynthesis`; saved `source_ref`, typed `intent` | Issue-backed alternate disconnection/realization, original-result preservation, rerun and verification. Fresh sessions replay the source to recover domain objects. |

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

Route settings follow `MultistepRetrosynthesisRequest`; see the bounded example
in [development_cases.json](../../examples/ai_native/development_cases.json).
Its `route_revision_followup` entry records the more substituted target used
after the initial target returned only one-step routes. Submit its `settings`
object through `plan_routes` to reproduce that search against a selected baseline.
An external agent selects a supported repair proposal and supplies its actual
`source_route_id`, `source_step_id`, `objective`, `method`, and `issue_ids` in
`intent`. Fabricated IDs and unsupported method combinations are rejected by
the owning domain contract. This revision adapter currently supports alternate
disconnections and realizations; condition-selectivity repair remains a direct
domain capability, not a supported workspace revision method.

## Custom analysis and persistent notes

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
registry audit, shared-core results awaiting independent review, and the supplied
recipe assessor classifying an unresolved reaction signature as a conflict.
The workspace preserves these contracts and records the limitations. It does
not treat incomplete validator coverage as proof that a proposed reaction is
impossible. Scientific fixes must be separately versioned and reviewed.

See [implementation status](Scientific_Workspace_Implementation_Status.md) for
the pilot findings, validation, and remaining release gates.
