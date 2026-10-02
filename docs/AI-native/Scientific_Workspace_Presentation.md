# Scientific workspace presentation

The default answer is a chemist-facing reaction card: **scheme → rationale →
precedent support**. Conditions use one target scheme, a preferred recipe, and
expandable alternatives. Routes use ordered cards with explicit intermediates.
The agent chooses the content and order; the UI does not rank chemistry.

## Where to change it

| Responsibility | Definition |
| --- | --- |
| Default answer content and level of detail | `chem_coworker/scientific_workspace/presentation/default.md` |
| How the agent authors and cites structured answers | `chem_coworker/scientific_workspace/agent_instructions/answer_authoring.md` |
| Typed answer fields and attribution checks | `chem_coworker/scientific_workspace/answers/answer_contracts.py` |
| Saved evidence binding and projection | `answers/step_precedents.py`, `answers/condition_precedents.py` for validation; `views/precedents.py`, `views/condition_precedents.py` for display |
| SVG generation and safe source links | `app/web_api/scientific_presentation.py` |
| Card order, expandable details and browser behavior | `app/web_api/scientific_chat.js` |
| Spacing, type, responsive layout and themes | `app/web_api/scientific_chat.css` |
| SVG layout and molecular drawing presets | `visualization/definitions/annotated_scheme.v1.json`, `visualization/definitions/render_styles.v1.json` |

Task guides own scientific strategy. Deterministic domain packages own chemistry,
compatibility and retrieval. Presentation consumes their recorded results and does
not introduce a second chemistry or scoring path. Explicit user requests may change
answer detail; attribution and evidence validation still apply.

## Answer and evidence contract

`scientific_answer.v2` adds two optional fields on each step:

- `rationale`: an attributed claim with `text`, `basis`, `source_ids`, `limitations`.
  Reported/computed rationale follows the existing citation requirements.
- `condition_precedent_refs`: saved, completed `inspect_condition_precedents` call
  references. Both sides must match the displayed reaction after canonicalization,
  including specified stereochemistry and chemical form.

Existing `precedent_refs` still identify route-support inspections. Condition links
do not replace required route inspections. A step may have both. Older v2 answers
default to no rationale and no condition links; the renderer does not invent them.
The runtime JSON schema declares these fields explicitly, using `null` and `[]`
when absent. `finalize_answer` supplies those defaults through the typed contract.

The condition inspection operation now has adapter contract version **2**. It saves
optional reference-catalogue metadata alongside the unchanged domain comparison.
No endpoint changes. The saved `condition_evidence_comparison.v1` result keeps its
existing chemistry schema with additive citation metadata. Catalogue access is shared
with route inspection through `source_catalogs.py`.

The data path is: saved answer → validated inspection references → saved source
observations → server-generated SVGs/links → browser cards. Rendering does not rerun
retrieval, consult the live catalogue, or edit evidence. Condition procedures join
by exact observation ID; unassigned reaction-level procedures are labeled separately.
Route template support may have several experiments, each displayed separately.
Source yields never become target yields through presentation.

## Visible and expandable information

Each card shows its scheme, basis, rationale, material cautions and first inspected
precedent. The source card includes reported conditions/yield, publication link,
recorded match/difference context and transfer limitations. Additional precedents,
full procedures, raw notation, record identifiers, comparison diagnostics and search
scope are expandable. Missing evidence, empty searches and failed evidence loading
remain distinct. Captured literature sources remain usable without local precedents.

For condition alternatives, only identical explicit reaction strings on independent
steps are grouped. The first step is visible; alternative recipes are collapsed.
Different transformations are not inferred equivalent. Routes retain each dependent
step and their original ordering. The first route opens initially; later refreshes
preserve the reader's expansion choices.

Target and source schemes share a readable display scale, with horizontal scrolling
on narrow screens. Important caveats are outside collapsed details. All source text
uses text nodes; the existing server sanitization owns rendered Markdown and URLs.

Restart the web server and refresh after upgrading. Saved answers receive the new
layout without another agent run. Existing scientific baselines detect the changed
adapter code; start a new investigation to generate new scientific calls. Later
presentation-profile edits can enter v2 investigations as recorded application context.
