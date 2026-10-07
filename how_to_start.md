# Start the local web interfaces

Run these commands in **PowerShell**, using the Python environment where the
project's chemistry dependencies are installed. Keep the server terminal open;
press **Ctrl+C** to stop it.

```powershell
Set-Location "C:\Git-softwares\Condition-agent"
```

## Quick start

After updating scientific tools, restart the agent server and create a new
investigation. The [literature preparation and compact inspection
workflow](docs/new/literature_preparation_and_console_efficiency.md) preserves
source evidence while reducing large console dumps and one-off summary scripts.
It also covers concrete-candidate selection, upstream source searches, proposed
recipe checks, captured scheme images and the new answer-attachment helpers.
The [agent log review fixes](docs/new/agent_log_review_fixes_20261004.md) add
source quantity-conflict checks, visible missing-atom diagnostics, nested input
help and explicit aliases for methyl iodide and acetic acid.

Choose one server command. The separate ports below let both interfaces run at
the same time in different PowerShell terminals.

| Interface | PowerShell command | Browser URL |
| --- | --- | --- |
| Regular research web UI, including fragment search | `python -m app.web_api --workbench --port 8000` | http://127.0.0.1:8000/ |
| Agent web UI / scientific workspace | `python -m app.web_api --scientific-chat --port 8011` | http://127.0.0.1:8011/scientific |
| Focused condition recommendation UI | `python -m app.web_api --port 8000` | http://127.0.0.1:8000/ |

The default command opens **Condition Desk**. Use `--workbench` to expose the
research tools, including **Fragment search**. These are different deployment
profiles; changing the browser URL alone does not switch profiles.

## First-time preparation

Use the existing project environment; `requirements-web.txt` supplies the web
dependencies, not a complete fresh chemistry environment. Python 3.10+ and RDKit
are required. Node.js 24.14.1+ is needed to build the regular React web UI.

```powershell
python --version
python -c "import rdkit; print(rdkit.__version__)"
python -m pip install -r requirements-web.txt
```

If this checkout uses a `.venv`, select its Python before running commands. For
example, put it first on PATH for the current PowerShell session:

```powershell
# Only when .venv already exists and contains the project dependencies.
$env:Path = "$(Join-Path (Get-Location) '.venv\Scripts');$env:Path"
```

Install the regular UI's JavaScript dependencies once, or after its lockfile
changes:

```powershell
node --version
npm.cmd --version
Push-Location .\web\reaction_recommender
npm.cmd ci
Pop-Location
```

`npm.cmd` avoids PowerShell script-shim execution-policy issues. All subsequent
Python commands below run from the repository root.

## Regular web UI and fragment search

Build the browser assets and start the workbench:

```powershell
python -m app.web_api --workbench --build --port 8000
```

Open **http://127.0.0.1:8000/**. Select **Fragment search**, click
**Cyclic ether example**, and then **Search fragments**. You can also enter your
own connected core as SMILES or SMARTS. Results include reaction drawings,
construction/retention evidence, source references and available procedures.
In SMILES mode, use **Draw** or **Edit drawing** to draw the
core in the molecule editor. Click **Use drawing**, then **Search fragments**.
SMARTS queries are edited as text.

To start from a **whole target**, enter or draw it in SMILES mode, click
**Suggest fragments**, inspect the highlighted regions, then click **Use candidate**
and **Search fragments**. Suggestions work without an index and never launch a
search automatically. You can still enter your own core directly.

To combine fragment discovery with **single-step retrosynthesis**, select
**Fragment-guided retro**. **Assisted fragment research** is the default workflow:

1. Enter or draw the full target, or load **Fragment research example**.
2. Click **Suggest strategic regions**, choose a region, or draw/type your own
   query in the separate core editor. The full target remains unchanged.
3. Click **Search chosen fragment** and inspect source reactions and procedures.
   If useful hits are missing, edit the query or choose **Suggest simpler queries**,
   record your reason, and search again. Searches never broaden automatically.
4. Select up to six source observations and click **Assess selected precedents on
   target**. Results distinguish local construction evidence, whole-source
   admission, structural differences and verified single-step proposals.
5. Export **research JSON** to retain query revisions, notes, source choices,
   results and errors. History holds 20 attempts in browser memory; export before
   closing or refreshing the page.

Discovery needs the prepared fragment index; transfer also needs the selected
operator library. The unrestricted baseline is optional in this workflow.
Choose **Automatic POC comparison** in the Workflow control to run the existing
automatic fragment selection and baseline comparison.
Operator libraries default to `results/operator_retrosynthesis_poc/full_scale_v3`;
set `CORE_RETROSYNTHESIS_LIBRARY_ROOT` to use another directory with `compact/` and
`full/` subdirectories containing `operator_library_v3.json.gz`.

Results show the baseline, direct source transfers, witness-directed proposals,
source admission failures, incomplete search exclusions and actual search work.
JSON export retains the complete result for comparing examples. Guided searches
use more total work; additional precursor sets are a descriptive coverage result,
not a validated accuracy gain. See the
[POC documentation](docs/AI-native/Deterministic_Fragment_Retrosynthesis_POC.md)
for evidence gates and the CLI for recorded investigations.

On later starts, omit `--build` unless the frontend source has changed:

```powershell
python -m app.web_api --workbench --port 8000
```

The default fragment index is
`results/ai_native/indexes/fragment_precedents.sqlite`. To use a different
**already prepared** index:

```powershell
python -m app.web_api --workbench --port 8000 --fragment-index "C:\data\fragment_precedents.sqlite"
```

Alternatively, set `$env:FRAGMENT_PRECEDENT_INDEX` before starting the server.
The web UI does not build the index during a search. See
[fragment index preparation](docs/AI-native/readme.md#prepare-and-search-fragment-precedents)
if it is missing. Interactive API documentation is at
**http://127.0.0.1:8000/api/docs**.

## Agent web UI / scientific workspace

The agent interface uses an installed, authenticated native Codex runtime. The
adapter discovers a suitable executable from PATH or supported Windows IDE
installations. See the [runtime requirements](docs/AI-native/readme.md) if
discovery or authentication fails.

```powershell
python -m app.web_api --scientific-chat --port 8011
```

Open **http://127.0.0.1:8011/scientific** and start a new conversation.
**This page does not need Node.js or a frontend build.** It is served directly
by the Python app.

Use **Workspace mode** in the top bar before sending the first message:

| Mode | Project tools/data | Task guidance | Required answer format |
| --- | --- | --- | --- |
| 1. Pure agent | Disabled | Disabled | None |
| 2. Tools and data | Enabled | Disabled | None |
| 3. Normal (default) | Enabled | Enabled | Existing structured answer |
| 4. Tools + formatting | Enabled | Disabled | Existing structured answer |
| 5. Tools + guidance | Enabled | Enabled | None |

The mode is saved with the conversation and every turn. Start a **New chat**
to change it; follow-ups keep the same mode and runtime thread when possible.
Old conversations remain in normal mode. The turns API accepts an optional
`mode` value (`pure_agent`, `tools_only`, `normal`, `tools_formatting`, or
`tools_guidance`); omitting it uses normal
for new conversations and the saved mode for follow-ups. A mode change within
an existing conversation returns HTTP 422.

Pure mode sends only the user question, runs outside the checkout, omits the
project Python path and disables local shell execution, local image reading,
MCP servers, plugins, personal memory and computer/browser integrations. Native
web search still follows the server's research profile. This is a model plus
native web-search baseline; it does **not** retain a general-purpose local shell.
The native tool dispatcher stays enabled so web search can execute. Local image
access uses the CLI's `features.view_image` flag. Mode policy v2 corrects an earlier
v1 dispatcher setting that prevented native tool calls despite requesting live
search. Restart the server and use a fresh chat when comparing against a v1 run.
A current native CLI with the required isolation controls is required; unsupported
CLIs fail instead of silently running with local access. These controls use the
[official Codex configuration reference](https://developers.openai.com/codex/config-reference).

Tools-and-data mode keeps deterministic scientific operations, dataset baselines
and evidence recording. Its prompt contains an API reference, without task
playbooks, learned advice, answer-authoring or presentation instructions. Saved
guidance is empty, and lesson recall/publication is disabled. Because this mode
deliberately has source access, it is an application-guidance ablation, not a
filesystem security boundary against an agent deliberately reading guide files.
Tools + formatting adds the same answer-authoring instructions, presentation
profile, structured handoff, evidence validation and single repair attempt as
normal mode, while keeping task playbooks and procedural learning disabled.
The formatting treatment therefore includes the existing evidence and self-review
requirements; it is more than a cosmetic change to the display.

Tools + guidance enables the same task guides and pinned procedural learning as
normal mode, without answer-authoring instructions, a presentation profile or a
required output schema. Pure agent, Tools and data, and Tools + guidance save
final text verbatim without scientific-answer validation or repair. Their answers
are marked as not validated. Markdown still renders in the browser in every mode.
Disabled instruction resources are omitted from each turn's saved context.
Modes without task guidance suppress inherited project instructions, personal
memory, skills and plugins; Tools + guidance follows Normal's native configuration.
Policy v3 records the independent guidance and formatting treatments. Restarting
on this version starts a new runtime thread when the saved configuration differs.

For comparisons, use the same questions, explicit model, reasoning/search settings,
deadline and dataset snapshot in fresh chats, and repeat trials. Mode 1 versus 2
measures adding project capabilities (including local execution). Modes 2, 4, 5
and 3 cover all four combinations of task guidance and structured output:
compare 2 versus 4 or 5 versus 3 for formatting; compare 2 versus 5 or 4 versus 3
for guidance. Guidance comparisons also include the inherited native configuration
described above. Modes 3 and 4 allow one answer-repair attempt, so compare total
usage and elapsed time as well as answer quality. These controls enable experiments;
they do not establish efficacy. Modes 3 and 5 retain procedural learning, so record
and hold their guidance and lesson snapshots fixed when comparing repeated trials.

The agent defaults to the `research` profile, which requests high reasoning,
live web search and a 1,800-second per-attempt deadline. It inherits the model
from the local runtime configuration. Example with a shorter deadline:

```powershell
python -m app.web_api --scientific-chat --port 8011 --agent-profile research --agent-timeout 1200
```

To select a native executable explicitly, replace the example path:

```powershell
python -m app.web_api --scientific-chat --port 8011 --codex "C:\path\to\codex.exe"
```

Use the actual native executable, not a `.cmd`, `.bat` or `.ps1` wrapper.

The default agent dataset configuration is
`examples/ai_native/artifacts.local.example.json`. It includes a `fragment_index`
entry. To use your own configuration:

```powershell
python -m app.web_api --scientific-chat --port 8011 --chat-artifacts "examples\ai_native\artifacts.local.example.json" --chat-root "results\ai_native\conversations"
```

`--fragment-index` configures the regular workbench. The **agent's** fragment
index is selected by `fragment_index` inside the `--chat-artifacts` JSON file.
After changing datasets or the configuration, restart the server and start a
new conversation so its investigation records the new baseline.

Condition screening is available directly in the agent conversation. For example:

> Generate 96 diverse weak-label screening recipes for
> `Brc1ccccc1.CN>>CNc1ccccc1`. Preserve full recipes, source IDs and missing
> values, and export JSON plus a CSV review table.

The built-in `generate_weak_label_screening_array` operation supports 1–250
recipes (default 24), returning fewer when the compatible recipe pool is smaller.
The default artifact configuration includes
`"weak_label_records": "datasets/weak_label/v2.1_cleaned.csv"`; custom configurations
must include this entry to enable screening. Its sibling
`v2.1_cleaned.condition_recipes.jsonl.gz` is automatically recorded and checked
in the investigation baseline. Screening works without the structural condition
indexes. Queries require a graph-supported transformation; source reaction
structures remain unverified, and historical yields are not yield predictions.
Keep these suggestions separate from the agent's structural recommendations.

After updating the code or artifact configuration, restart the agent server and
start a new conversation. Existing investigations retain their frozen baseline.

Conversations, investigation artifacts and per-turn `progress.jsonl` debugging
logs are saved under `results/ai_native/conversations` by default. The agent
server supports loopback hosts only; keep the default `127.0.0.1`.

## Both interfaces from one server

`--scientific-chat` also enables the regular research workbench. After installing
the JavaScript dependencies, this command prepares and serves both:

```powershell
python -m app.web_api --scientific-chat --build --port 8011
```

- Regular workbench: **http://127.0.0.1:8011/**
- Agent workspace: **http://127.0.0.1:8011/scientific**
- API documentation: **http://127.0.0.1:8011/api/docs**

Omit `--build` on subsequent starts when the frontend has not changed.

## Restart and troubleshooting

- **Stop/restart:** press Ctrl+C in the server terminal, then rerun the selected
  command. These launch commands do not enable automatic Python reload.
- **Port already in use:** stop the server you started on that port, or choose
  another port, for example `--port 8012`, and update the browser URL.
- **Fragment search is missing:** restart with `--workbench` (or
  `--scientific-chat`). Rebuild with `--build` after frontend changes and refresh
  the browser.
- **Frontend build is missing:** run the JavaScript preparation above, then add
  `--build`. The agent page at `/scientific` can run without that build.
- **Fragment index unavailable:** check the prepared SQLite path; configure
  `--fragment-index` for the workbench or the artifact JSON for the agent.
- **Python import errors:** confirm that `python` refers to the project's
  environment, then install the missing dependencies there. `requirements-cli.txt`
  alone is not the web dependency set.
- **Agent runtime unavailable:** check the native executable and its
  authentication, or use `--codex` to select it explicitly.
- **`helper_unknown_error: setup refresh had errors`:** inspect the newest
  `%USERPROFILE%\.codex\.sandbox\sandbox.*.log` for the underlying Windows
  error. The 2026-10-07 conversations failed because a running process locked
  `node_repl.exe` during runtime read/execute permission validation (OS error 32).
  Finish active Codex tasks, close unused Codex sessions/apps holding that runtime,
  and restart the agent server before retrying. A new chat alone does not release
  the lock; do not disable sandboxing to hide this failure.
- **`Scientific runtime changed` in a fresh chat:** compare the recorded/current
  fields in the error and ensure the server and worker use the same Python and
  installed dependencies. Windows baseline identity now uses Python's build
  architecture, avoiding WMI and missing `PROCESSOR_*` variables in restricted
  workers. Restart the server and start a new investigation after this fix;
  existing baselines remain frozen. Older errors omit the differing fields and
  cannot establish which dependency or platform field differed.
- **Agent exits with a model rejection:** open the failed turn's details. The
  runtime now reports the provider error from `runtime.jsonl`, even when stderr
  is empty. An older CLI can reject a model selected in your desktop Codex
  configuration. Restart with a current native executable using `--codex`, or
  explicitly choose a model supported by that CLI's login, for example:

  ```powershell
  python -m app.web_api --scientific-chat --port 8011 --agent-model gpt-6-sol
  ```

  Model availability depends on the account and provider. The adapter does not
  automatically substitute another model. After a code update, start a new
  investigation to record the current baseline.

List all startup options:

```powershell
python -m app.web_api --help
```

More details: [regular web UI](web/reaction_recommender/README.md) ·
[scientific workspace](docs/AI-native/readme.md).
