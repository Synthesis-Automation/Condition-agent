# ZBS chemistry recommender: focused recommendations and automation handoff

ZBS chemistry recommender is the default browser page. It combines the canonical generic
and weak-label recommenders into one chemist-facing shortlist. Retrosynthesis,
forward synthesis, and feature-analysis controls are absent from this page.

## Run

Build the frontend from the repository root:

```powershell
cd web/reaction_recommender
npm ci
npm run build
cd ../..
python -m app.web_api
```

Open `http://127.0.0.1:8000/`. This is the single launcher; `app.conditions_web`
has been removed. The default CLI and ASGI application both use the explicit
recommendation route allowlist.

For development and parity review, stop the server and restart with
`python -m app.web_api --workbench`. This selects the research frontend at `/`
and enables its matching API. Old `/index.html` and `/workbench.html` URLs
redirect to the selected page. Only compiled `/assets` are mounted as static
files; the alternate HTML page cannot be opened against an incompatible API.
HTML responses use `Cache-Control: no-store` to avoid stale entry documents.
Run `python -m app.web_api --build` when browser assets need rebuilding.

Linux process-manager command (put an HTTPS reverse proxy in front):

```bash
python -m uvicorn app.web_api.main:app --host 127.0.0.1 --port 8000 --workers 1
```

Inside a private container network, bind to `0.0.0.0` instead. Configure the
existing runtime data paths:

```bash
export CONDITION_RECOMMENDER_LIBRARY=/srv/conditions/data/literature
export CONDITION_RECOMMENDER_LIBRARY_MODE=full
export CONDITION_RECOMMENDER_WEAK_LABEL_RECORDS=/srv/conditions/data/weak_label/v2.1_cleaned.csv
```

The literature root contains `full/generic_index.sqlite` and optionally
`compact/generic_index.sqlite`. Shared reaction core v2 is the default in both
the focused page and the workbench; each SQLite index needs its matching
`generic_index.shared_core.sqlite` companion in the same directory. Normal builds
generate both files, and the current Full/Compact libraries already contain them.
See [default selection and validation status](shared_reaction_core_default_20260913.md).
Keep reference and experimental-detail catalogs
beside their index. Keep the weak-label catalog
`v2.1_cleaned.condition_recipes.jsonl.gz` beside its CSV.

RXNMapper is used when available. Set `CONDITION_RECOMMENDER_USE_RXNMAPPER=false`
to disable it for this combined endpoint; this can reduce query coverage.
No model installation or chemistry-free mapping substitute is performed.

## Chemist workflow

1. Draw or paste a product-specified reaction, or select **Try an example**.
   The example picker randomly chooses a complete Suzuki, C–N coupling, or amide
   formation involving heteroaryl, protected, or substituted substrates. It avoids
   repeating the current example and displays the selected reaction type.
2. Confirm missing fragment sources when the existing completion workflow asks.
3. Choose the maximum recipes per tab and the reaction match scope, then
   select **Find conditions**. Automatic broadening uses the same canonical
   search as Workbench. Results appear in **Literature-based** and
   **Screening suggestions** tabs with source-specific recipe counts.
4. Select a row in the table on the left to view its materials, reported amounts,
   temperature, time, and atmosphere on the right. All returned recipes are
   listed in source rank order. Each tab remembers the last viewed row; the first
   recipe opens initially. Multi-stage records retain their original order.
5. Review **Details & references**, open by default, for structures, references, procedures,
   compatibility evidence, historical yield summaries, cautions, and missing
   preparation details.
6. Copy or export the viewed condition set, or use table checkboxes to select
   intact recipes for a screening-selection JSON bundle. Download full results JSON
   for the complete source audits, including failures and abstentions.

The page uses a light gray background, white cards, subtle borders, and green
actions, following the restrained style of Google Classroom. The shared
Workbench drawing editor opens directly below the header; there is no separate
hero, decorative reaction, or introductory step list. The results table summarizes
catalysts/reagents, solvents, and reported historical yields; the detail pane
shows the complete selected recipe. On narrow screens the table scrolls within
its panel and details appear below it. Tabs support arrow, Home, and End keys;
row-number buttons also work with the keyboard.
The chemistry summary retains the structure-derived reaction equation and
reactive-group labels; nearby groups can be expanded when available. It does
not infer a named reaction from the displayed text.

Screening suggestions use the canonical weak-label engine's screening mode,
including its existing recipe-diversity selection. They remain explicitly
unverified source-structure evidence. The screening tab remains available even
when structural results fill their result limit. Each tab preserves its own
engine rank; a recipe supported by both sources appears in both tabs but
shares one selection and is exported once. Cross-source selection remains the
chemist's choice and never generates new reagent combinations.
If only screening options are available, that tab opens automatically. Empty
sources retain their status and show no stale recipe details. A new search clears
the viewed rows and export selection.

## Source combination contract

`POST /api/v1/conditions/recommend` accepts:

```json
{
  "reaction_smiles": "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
  "completion_choices": [],
  "top_k": 10,
  "search_scope": "automatic"
}
```

The application calls the canonical engines; chemistry filters, admissibility,
and scoring remain owned by their existing packages. Generic results precede
weak-label results, retaining engine rank. Scores are never averaged or
compared across engines. Policy bounds and evidence descriptions live in
`condition_recommender/definitions/combined_recommendation.v1.json`.

`top_k` is optional (1–50 recipes per source, defaulting to the existing policy
when omitted). `search_scope` accepts `same_handle`, `automatic`, or `broad`
and defaults to `automatic`; it controls structural retrieval. The weak-label
request always uses `weak_label_screening`. No research routes are exposed.

Only identical registry-identified full recipe payloads for the same effective
query are coalesced. Equality includes quantities, stages, and uncertainty;
core identity alone does not establish procedural equivalence. Distinct source
provenance can conservatively keep apparently equivalent recipes separate.
Every supporting recommendation and its warnings remain independently
attributed. No support counts are summed across sources.

Generic review or fallback results are described as broader structural matches,
not verified signatures. Weak-label suggestions retain their evidence
qualification and their engine's requirement for a structurally supported
query. With explicit completion choices, weak-label retrieval is skipped
because that engine does not support those choices. The skip is reported.

An unavailable source does not erase valid results from another source. A
source that abstains cannot contribute any stale recommendations. Library
failures appear as source status; private filesystem paths are only logged on
the server. The source outputs in full JSON retain typed uncertainty and
definition-version evidence supplied by the canonical engines.

## Automation JSON: planning input, not robot commands

Each **Export recipe** file conforms to
`condition_recommender/definitions/automation_handoff.v1.schema.json` and embeds
the existing registry `synthesis_protocol.v1.schema.json` contract. It contains:

- deterministic handoff and condition-option identifiers;
- original query and effective reaction used for the protocol;
- parsed reaction inputs and intended outputs;
- resolved condition identities, CAS numbers, reported amounts and units;
- the intact canonical recipe, including definitions, stages, and provenance;
- reported operating setpoints and the existing `maintain_conditions` stages;
- source-specific recommendation evidence and available source procedures;
- missing required fields and an explicit execution status.

All exports have `execution_ready: false`,
`execution_status: "requires_robot_adapter_and_review"`, and `robot_target: null`.
The nested protocol retains `execution_readiness: "review_required"`.
The handoff does not invent reaction scale, absolute dosing quantities,
addition order, vessel settings, mixing, quench, or workup. Reported quantities
are precedent data, not automatically scaled instructions for a new substrate.
No raw procedure text is treated as executable code.

The screening-selection JSON conforms to
`condition_recommender/definitions/screening_selection.v1.schema.json`. It is
a versioned envelope containing the selected
individual handoffs under `handoffs`; each handoff keeps its own checksum and
provenance. Selection never joins multiple recipes into one experiment.

To actually run a robot, implement and validate an adapter for its specific
command schema, inventory identifiers, equipment, and limits. This requires
the robot/software contract and target scale. That adapter must reject
unresolved preparation fields and must not interpret these downloads as
already authorized, executable jobs. There is no robot dispatch endpoint.

## Validation and rollout

Run `pytest -q` and the frontend production build. New regression tests cover
source ordering, exact deduplication, differing operating variants, review and
conflicting evidence, partial failures, completion choices, protocol fidelity,
deterministic exports, JSON schema validation, and research-route isolation.

Dataset contents and chemistry definitions are unchanged. The focused endpoint
now selects the existing weak-label screening mode, and forwards the chemist's
structural scope and result limit to the canonical engines. This UI does not
establish production chemistry validity: the
roadmap's independent review and untouched evaluation gates still apply.

The 2026-09-16 refresh was checked with the real Full and weak-label libraries:
the supplied indole chloride/heteroaryl boronic-acid query returns ten structural
recommendations and ten weak-label screening suggestions. Exact recipe
coalescing yields nine distinct structural options and ten screening options.
Desktop and 390-pixel mobile screenshots are saved under
`results/zbs_tabs/`. Browser checks cover all three example choices,
no immediate example repeats, the page title, real results, tab switching,
per-source row order and remembered details, selection/export deduplication,
keyboard controls, search resets, partial source failures, and mobile overflow.
The focused API/combination/weak-label suite has 26 passing tests.
All three displayed examples pass structure and completion validation; the
C–N and amide examples also return both structural and screening results from
the local libraries. The complete `pytest -q` suite passes (1,546 tests in
396.09 seconds), as do all eight browser checks, the frontend production build,
and Python name checks. Browser checks used installed Edge in headless mode;
no additional browser installation was required.

If the installed SQLite library uses an older schema, the result reports
`INDEX_REBUILD_REQUIRED`; rebuild recommendation artifacts using the current
canonical conversion/index tooling before using that library. The application
does not rewrite metadata to bypass chemistry validation.
