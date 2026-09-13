# Reaction Condition Recommender Web UI

Use **one launcher**: `python -m app.web_api`. It serves **Condition Desk** and
the matching recommendation API at `http://127.0.0.1:8000/`.
The separate `app.conditions_web` launcher has been removed.

For research tools, stop the running server and use
`python -m app.web_api --workbench`. The same root URL then serves the research
workbench with its API enabled. `/index.html` and `/workbench.html` redirect to
the selected root page, so a restricted server cannot open an incompatible UI.
Run only one server on port 8000. Use `--build` to rebuild browser assets before
startup, or `--host` and `--port` to change the listening address.

The new UI downloads full result JSON, individual automation protocol handoffs,
and selected screening sets. Automation exports are explicitly planning drafts:
they preserve missing quantities and operations and require a robot-specific
adapter before execution. See
[`docs/new/recommendation_only_web_ui.md`](../../docs/new/recommendation_only_web_ui.md)
for the API, environment configuration, evidence policy, and export contract.

Local React/Ketcher client for the repository's canonical structure-first
condition recommender. The browser only owns interaction and presentation; all
reaction analysis, chemistry filtering, retrieval, scoring, recipe
aggregation, and SVG rendering stay in the standalone Python packages behind a
versioned FastAPI boundary.

Prerequisites are Python 3.10+ with the repository chemistry dependencies and
Node.js 24.14.1+.

## Condition search and priorities

In the Workbench, **Search scope** controls retrieval independently of ranking:

- **Same reactive handle** excludes the new related-handle tier and stops the
  signed search before broad core/edit-graph neighbours.
- **Automatic broadening** (default) checks permitted related handles when the
  close same-handle pool has fewer than two independent evidence units.
- **Broader analogues** checks permitted related handles even when close support
  is sufficient. It still applies chemistry and condition compatibility gates.

The initial rule covers a single observed aromatic C-I/C-Br substitution in
either direction, with a formed C-C, C-N, C-O or C-S bond. The partner chemistry,
remaining edits, center states and departing fragments must match. It does not
enable a universal halogen reactivity ladder or substitute a different reagent
source. Ambiguous/unsigned queries retain their qualified fallback paths; they
cannot generate these verified related-handle hypotheses.

**Prioritize** offers Balanced, Closest chemistry and Strongest supporting
evidence. **Advanced weights** retains detailed relative weights. The advanced
**Independent evidence target** changes the support threshold, not the number
of duplicate source rows required. Top results is a maximum: a supported
analogue tier may return fewer recipes instead of forcing a broader search.

Results show ordinal match levels: **1 Same reaction**, **2 Close precedent**,
**3 Related handle**, **4 Broader analogue**. These are evidence distances, not
success probabilities. Level 1 requires the same normalized reactants and
products; matching reaction facets alone is Level 2. Broader matches remain
below closer matches regardless of ranking weights. The display limit can hide
lower-tier recipes; the result JSON retains the attempted retrieval tiers.

The request API accepts `search_scope: "same_handle" | "automatic" | "broad"`.
Recommendation result schema 4.1 adds `search_scope`, `match_level`,
`match_label`, `match_details`, and broadening provenance in retrieval traces.
Existing SQLite index schema 6.5 is reused; **no dataset or index rebuild is
needed**. Restart the server after updating the code and browser assets:

```powershell
python -m app.web_api --workbench --build
```

Implementation and validation details:
[related-handle retrieval and simpler ranking](../../docs/new/related_handle_retrieval_20260913.md).

## Development

From the repository root, start the API:

```powershell
python -m pip install -r requirements-web.txt
python -m app.web_api
```

In another terminal, start the client:

```powershell
cd web/reaction_recommender
npm install
npm run dev
```

Open `http://127.0.0.1:5173/`. Vite proxies `/api` to the local service on port
8000. The default recommendation index is
`datasets/literature/full/generic_index.sqlite`; override it with
`CONDITION_RECOMMENDER_INDEX` or `python -m app.web_api --index <path>`.
The weak-label mode uses `datasets/weak_label/v2.1_cleaned.csv` and its
paired recipe catalog by default. Override the CSV with
`CONDITION_RECOMMENDER_WEAK_LABEL_RECORDS`.

The experimental **Coupled two-step strategies** mode is enabled when both the
validated-departures route-step operator library and frozen v1 strategy panel
are available. Override their default result paths with
`CORE_RETROSYNTHESIS_COUPLED_LIBRARY` and
`CORE_RETROSYNTHESIS_COUPLED_PANEL`.

## Single-port local build

```powershell
cd web/reaction_recommender
npm run build
cd ../..
python -m app.web_api
```

Open `http://127.0.0.1:8000/`. FastAPI serves the compiled client and the
versioned `/api/v1` routes from the same local process. Interactive API
documentation is available at `http://127.0.0.1:8000/api/docs`.

## Drawing editor browser checks

After building the frontend and starting `python -m app.web_api`, run from this
directory:

```powershell
npx playwright install chromium
npm run test:e2e
```

Set `WEBUI_TEST_URL` for a different local server port. For a server started with
`--workbench`, also set `WEBUI_TEST_PROFILE=research_workbench`. To use an installed Edge
or Chrome instead of downloading Chromium, set `BROWSER_CHANNEL` to `msedge` or
`chrome`. The checks cover drawing import/export on both pages and interrupted
editor downloads or render failures. An editor failure now stays inside the
dialog, with recovery guidance, while existing SMILES input remains available.

Refresh an already-open browser tab after rebuilding the frontend so it uses
the current asset filenames.

## Research workbench workflow (`--workbench`)

The precedent selector displays the loaded index's actual record count. An
explicit custom `--index` is labeled “Custom index”, so a development sample
cannot be mistaken for Full. For the opt-in Full/Compact shared-core artifacts,
see the [larger-corpus validation report](../../docs/new/shared_reaction_core_rollout_20260913.md).

- Draw, clear, load, paste, and export reaction SMILES with Ketcher.
- Validate product-fragment source requirements before recommendation.
- Retrieve and rank chemically compatible canonical condition recipes.
- Run a separate weak-label condition mode as either ranked fallback recipes
  or a diversity-selected screening array. The UI shows the graph-derived
  reaction-type hint, matched participants, source-label evidence, and
  persistent warnings that the precedent reactions are not structure verified.
- Apply declarative ranking profiles or transparent custom weights.
- Predict possible products from dot-separated starting materials with the
  target-blind Forward synthesis mode. A canonical condition-recipe JSON object
  can be supplied to apply hard compatibility checks before ranking.
- Include or exclude intermolecular self-reactions. When enabled, a
  bifunctional input can occupy multiple operator roles as separate assumed
  equivalents; these pathways are labelled, stoichiometrically traced, and
  modestly penalized rather than treated as ordinary single-equivalent paths.
- Select a structured condition profile instead of writing JSON: reaction
  strategy, transition-metal family, redox environment, and medium are loaded
  from a versioned backend catalog. Every structural ranking adjustment and
  uncertainty notice appears in the candidate evidence. Canonical recipe JSON
  remains available under Advanced options for expert, substance-resolved
  compatibility checks.
- Inspect forward graph-validation evidence, alternative operator/template
  pathways, competition groups, atom correspondence, source support, and route
  audit disposition, or export the complete result as JSON.
- Analyze a molecule or reaction with the same deterministic featurization used
  by the Qt tool, including motifs, reactive sites, reaction-core evidence,
  partner roles, mapping provenance, and the canonical nested analysis.
- Inspect score traces, structural matches and mismatches, cautions, conditions,
  yields, fallback levels, and precedent provenance.
- Test promoted v1 two-step operator pairs against an arbitrary target. Each
  logical strategy exposes its intermediate, terminal precursors, two physical
  reactions, validation statuses, training support, and ordinary one-step
  fallbacks. This mode is explicitly experimental and requires chemist review.
- Export the complete versioned recommendation or feature result as
  JSON.
- Preview and download the selected recommendation's versioned synthesis
  protocol JSON, including registry substance IDs, CAS numbers, quantities,
  operating conditions, observed operations, and execution-readiness gaps.

The forward endpoint uses a prebuilt `forward_operator_library_v1.json.gz` next
to the selected retrosynthesis library when available. Otherwise the API derives
and process-caches a source-round-tripped forward library on first use; that
first request can take longer than later requests.

The UI intentionally exposes no arbitrary file paths or upload endpoints. Local
dataset identity and access remain server configuration concerns.
