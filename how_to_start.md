# Start the local web interfaces

Run these commands in **PowerShell**, using the Python environment where the
project's chemistry dependencies are installed. Keep the server terminal open;
press **Ctrl+C** to stop it.

```powershell
Set-Location "C:\Git-softwares\Condition-agent"
```

## Quick start

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
