# Interactive route planning

Implemented 2026-10-07 as a Workbench mode. This is a human-selected route
editor over the existing single-step engine. Chemistry admission, operator
definitions, ranking, and independent evaluation release gates are unchanged.

## Use

Launch `python -m app.web_api --workbench --build`, then select **Interactive
route planning**, the last analysis mode. Draw or paste one connected target and
click **Start planning**.
Select a molecule in the tree and click **Find disconnections**. Each viable
reaction proposal appears as a numbered reaction node beneath that molecule;
its required precursor molecules appear beneath it. Strategies and their
alternate precursor choices are separate reaction alternatives. Search adds
alternatives without selecting a route. Select any precursor, including one on
an unselected branch, to continue exploring it. Equal molecules on separate
branches remain independently editable.

Click a numbered reaction node or **Use this step** to select that reaction
and its connecting path back to the target. Green edges and reaction nodes
highlight the selected route; a blue molecule outline identifies the inspector's
current molecule. Selecting another alternative preserves all explored branches
and their local choices. The right panel retains candidate details, reported
validation, conditions, and stock evidence for the inspected molecule.

**Clear route choice** deselects the molecule's reaction while keeping every
alternative. **Remove selected alternative** deletes only that reaction and its
precursor subtrees; other alternatives remain. **Undo** and **Redo** recover
search expansions and route edits, including starting-material designations.
A new edit clears redo history. **Use as starting material** is an explicit
decision after clearing a route choice; it retains explored alternatives but
ends the selected route at that molecule. **Reopen for planning** reverses it.
Empty search results never terminate a branch.

**Find step conditions** calls the existing condition service for the selected
reaction only. **Check supplier stock** is available when a local supplier
portfolio is configured. It retains exact matches, source snapshot provenance,
and check time separately from the user's stopping choice. Saved stock checks
are not live inventory. Stock matches do not automatically change the route.

The tree runs from top to bottom, alternating molecule and reaction nodes.
Reaction alternatives are choices; precursors under one reaction are all
required. Molecule cards and structure previews keep a fixed size as branches
grow or the window narrows. Scroll horizontally and vertically to explore the
tree, collapse alternatives to reduce its extent, or choose **Focus molecule**
to return to the current molecule. The view keeps that molecule visible after
layout changes. On narrow screens the candidate panel appears below the tree.

## Persistence and evidence

The current plan autosaves in this browser at the same origin. A different port
or browser has different storage. Storage failure is shown with an export
instruction; edits remain available in memory. **Export session** retains the
selected route, all explored alternative subtrees, search settings, condition and stock
evidence, and undo/redo history. **Import session** validates molecular identity
and topology before replacing the current plan. An invalid import preserves the
current plan. **New target** starts a new session, so export first to keep both.

**Export route** produces the selected canonical `ReactionRouteTree` schema 2.0.
It marks steps as predicted and terminal leaves as user-designated. The full
session remains the artifact for detailed search and condition evidence.

Session validation checks every alternative, even outside the selected route.
It reconciles candidate products, complete precursor sets,
reaction sides, occurrence identities, and ancestor cycles using parsed molecular
graphs. Repeated sibling molecules are allowed. Search-reported signature checks
and all other evidence remain reported evidence, including after import; restoring
a session does not rerun operator reconstruction, forward prediction, conditions,
or stock checks. The UI discloses this and offers explicit refresh actions.

Search results are retained with normalized settings and can supply choices for
an identical molecule elsewhere in the session. **Search again** explicitly
reruns the engine and adds new alternatives, retaining existing branches.
An identical result does not duplicate alternatives or erase their edits.
Changing search settings does not erase chosen steps or their original
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
`clear`, `remove`, `stop`, `reopen`, `undo`, `redo`, `conditions`, and `stock`. `select`
references a retained `search_id`, `strategy_index`, and `realization_index`;
it does not accept an arbitrary replacement reaction. The response uses the
existing API envelope and contains `session`, `route_tree`, and `summary`.
Invalid actions or molecular/topological conflicts return 422; unavailable
runtime resources return 503. The endpoint is absent from Condition Desk's
restricted deployment.

Sessions now use `interactive_planning.v2` and browser exports use
`interactive_planning_browser.v2`. Each molecule stores `alternatives` (choice
reference plus children), one optional selected `choice`, and
`expansion_warnings`. There is no duplicate selected-children field in JSON.
Canonical route export projects only the selected alternatives. The `clear`
action is new; `select` preserves alternatives and activates the connecting
ancestor choices; `remove` deletes only the selected alternative.

The importer migrates v1 sessions and all undo/redo roots into v2. Existing
browser storage is read at its original key and rewritten as v2 on restore.
Only v2 is emitted and edited; there is no separate v1 runtime. Keep the one-way
reader while supported user exports may contain v1; remove it only after a
documented import-format retirement, retaining the migration regression fixture.
Reaction signatures, chemistry definitions, canonical route schema, and search
APIs are unchanged.

Limits remain 200 molecule occurrences across **all** alternatives, 20 reaction
levels, 30 undo and redo states, and 100 retained searches. Search retains all
candidate evidence but does not attach an alternative that would create an
ancestor cycle or exceed tree limits. The inspector lists these omissions;
no incomplete precursor set is attached. Session responses are limited to
18 MB and imports to 20 MB; a rejected action leaves the previous session intact.
Browser storage may impose a smaller limit, reported by the autosave notice.

## Validation

Backend regressions cover alternative preservation, ancestor-path selection,
selected-only export, clearing and removal, v1 migration, duplicate
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
