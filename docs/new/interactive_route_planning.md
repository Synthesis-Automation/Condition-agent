# Interactive route planning

Implemented 2026-10-07 as a Workbench mode. This is a human-selected route
editor over the existing single-step engine. Chemistry admission, operator
definitions, ranking, and independent evaluation release gates are unchanged.

## Use

Launch `python -m app.web_api --workbench --build`, then select **Interactive
route planning**. Draw or paste one connected target and click **Start planning**.
Select a molecule in the route tree, click **Find disconnections**, inspect the
strategies and their alternate precursor choices, and click **Use this step**.
Every precursor is added as a separate molecule occurrence. Select any precursor
to continue. Equal molecules on separate branches remain independently editable.

The right panel retains the chosen step and candidate alternatives. Changing
the choice replaces only that molecule's downstream branches. **Undo** and
**Redo** recover route edits, including starting-material designations. A new
edit clears redo history. **Remove step and its branches** leaves that molecule
unresolved. **Use as starting material** is an explicit user decision; **Reopen
for planning** reverses it. Empty search results never terminate a branch.

**Find step conditions** calls the existing condition service for the selected
reaction only. **Check supplier stock** is available when a local supplier
portfolio is configured. It retains exact matches, source snapshot provenance,
and check time separately from the user's stopping choice. Saved stock checks
are not live inventory. Stock matches do not automatically change the route.

The tree collapses by branch and scrolls normally over molecule drawings. On
narrow screens the candidate panel appears below the route tree.

## Persistence and evidence

The current plan autosaves in this browser at the same origin. A different port
or browser has different storage. Storage failure is shown with an export
instruction; edits remain available in memory. **Export session** retains the
route, alternative results, search settings, reported condition and stock
evidence, and undo/redo history. **Import session** validates molecular identity
and topology before replacing the current plan. An invalid import preserves the
current plan. **New target** starts a new session, so export first to keep both.

**Export route** produces the selected canonical `ReactionRouteTree` schema 2.0.
It marks steps as predicted and terminal leaves as user-designated. The full
session remains the artifact for detailed search and condition evidence.

Session validation reconciles candidate products, complete precursor sets,
reaction sides, occurrence identities, and ancestor cycles using parsed molecular
graphs. Repeated sibling molecules are allowed. Search-reported signature checks
and all other evidence remain reported evidence, including after import; restoring
a session does not rerun operator reconstruction, forward prediction, conditions,
or stock checks. The UI discloses this and offers explicit refresh actions.

Search results are retained with normalized settings and reused for an identical
molecule and settings within a session. **Search again** explicitly reruns the
engine. Changing search settings does not erase chosen steps or their original
evidence. Mode changes, molecule selection, and cancellation discard late browser
responses. Cancelling a request stops waiting; an already-running synchronous
server search may continue to completion without modifying any server session.

## Ownership and API

- `core_retrosynthesis/interactive_planning.py` owns deterministic session edits,
  graph identity checks, occurrence/cycle rules, and canonical route export.
- `app/web_api/planner.py` composes those edits with the existing runtime search,
  condition, and exact-stock services.
- `InteractivePlanner.tsx` owns presentation, selection, autosave, import/export,
  and cancellation. It contains no molecular identity or reaction-admission rules.

`POST /api/v1/retrosynthesis/planner` accepts an `action` plus a session and
action-specific arguments. Actions are `start`, `restore`, `search`, `select`,
`remove`, `stop`, `reopen`, `undo`, `redo`, `conditions`, and `stock`. `select`
references a retained `search_id`, `strategy_index`, and `realization_index`;
it does not accept an arbitrary replacement reaction. The response uses the
existing API envelope and contains `session`, `route_tree`, and `summary`.
Invalid actions or molecular/topological conflicts return 422; unavailable
runtime resources return 503. The endpoint is absent from Condition Desk's
restricted deployment.

New schemas are `interactive_planning.v1` for domain sessions and
`interactive_planning_browser.v1` for the browser export wrapper. Existing
reaction signatures, definitions, route schema, and search APIs are unchanged.
The first release allows 200 molecule occurrences, 20 reaction levels, 30 undo
and redo states, and 100 retained searches. Session responses are limited to
18 MB and imports to 20 MB; a rejected action leaves the previous session intact.
Browser storage may impose a smaller limit, reported by the autosave notice.

## Validation

Backend regressions cover branch replacement, sibling preservation, duplicate
molecule occurrences, cycle prevention, contradictory candidate structures,
invalid imports, deterministic route identity, and user stopping semantics.
Real Suzuki, C–N, C–O, and C–S operator searches round-trip through selection and
restore. HTTP tests exercise deployment isolation and exact supplier evidence.

After `npm run build`, run `npm run test:planner` in
`web/reaction_recommender`. This starts a dedicated local fixture server using
the real session API and predictable search results. Set `BROWSER_CHANNEL` to
`msedge` or another installed Playwright-compatible browser as needed. Browser
tests cover three-step branching, alternative choices, undo/redo, save/reload,
import/export, cycles, empty results, cancellation, errors, and responsive layout.
Fixture candidates are UI test data and make no chemical validity claim.
