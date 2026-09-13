# Condition Desk: focused recommendations and automation handoff

Condition Desk is the default browser page. It combines the canonical generic
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
Condition Desk and the workbench; each SQLite index needs its matching
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

1. Draw or paste a product-specified reaction; use the example to explore.
2. Confirm missing fragment sources when the existing completion workflow asks.
3. Select **Find conditions**. Three options are shown initially; more can be
   revealed without rerunning chemistry.
4. Review materials, reported amounts, temperature, time, atmosphere, and
   cautions. Multi-stage records remain ordered and are never flattened.
5. Expand **Precedents & details** for structures, references, procedures,
   compatibility evidence, historical yield summaries, and qualifications.
6. Copy a condition set, download its automation handoff, or select intact
   recipes for a screening-selection JSON bundle. Download full results JSON
   for the complete source audits, including failures and abstentions.

The initial screening feature is a chemist-selected set of retrieved recipes.
It does not claim to optimize experimental diversity or generate new reagent
combinations. A cross-source diversity selector is a later, separately
validated recommendation policy.

## Source combination contract

`POST /api/v1/conditions/recommend` accepts:

```json
{
  "reaction_smiles": "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
  "completion_choices": []
}
```

The application calls the canonical engines; chemistry filters, admissibility,
and scoring remain owned by their existing packages. Generic results precede
weak-label results, retaining engine rank. Scores are never averaged or
compared across engines. Policy bounds and evidence descriptions live in
`condition_recommender/definitions/combined_recommendation.v1.json`.

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

Each **Download automation JSON** file conforms to
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

Dataset contents, chemistry definitions, and existing retrieval policies are
unchanged. This UI does not establish production chemistry validity: the
roadmap's independent review and untouched evaluation gates still apply.

If the installed SQLite library uses an older schema, the result reports
`INDEX_REBUILD_REQUIRED`; rebuild recommendation artifacts using the current
canonical conversion/index tooling before using that library. The application
does not rewrite metadata to bypass chemistry validation.
