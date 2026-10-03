Saved answer contract and handoff:

Do not read answer-schema.json at startup. The runtime supplies the schema for
validation issues. Follow the runtime's final handoff instructions when finished.

Prepare scientific_answer.v2 JSON. answer_markdown directly answers the question,
cites relevant sha256:<64 hex> artifacts and preserves limitations. evidence_refs
must identify actual call, derived_file, replay, custom_execution, literature_source
or literature_excerpt artifacts, not a note, your answer or an invented reference.
Every claim about local scientific results needs recorded evidence. uncertainties
lists material limitations; needs_user_input is true when clarification is needed.

Use schema_version='scientific_answer.v2'. Include sources, molecules,
target_molecule_ids, steps, routes and claims (empty arrays when irrelevant). Use
stable short IDs to connect objects. Every molecule, step, condition, yield and claim
has a basis: input, reported, computed, proposed or unknown. Reported/computed objects
need source_ids. Preserve domain warnings in limitations. A proposal recorded by a
tool remains proposed; its normalization/compatibility computation is computed.

Local sources use kind=local_artifact, artifact_ref of the recorded call/attachment,
url=null and locator of the inspected record/field/step. Computed objects require a
completed call, replay or run_python custom_execution, not notes or an attached output
alone. An attachment can support discussion of derived analysis with its execution
provenance limitation. Never use an unrelated call to satisfy provenance or relabel a
calculation as a reported yield.
External sources use kind=external_source, an original/final captured URL, an exact
locator (e.g. Example 1), and artifact_ref of literature_source/literature_excerpt
with captured text. A URL fragment can point to the example. Legacy attached source
captures retain URL, retrieval date, excerpt and provenance. A remembered source is
not an inspected source. Citation establishes attribution, not independent review.

Use molecule IDs for explicit reactants/products; do not invent missing structures
to complete a drawing. List steps in dependency order. after_step_ids identifies
preceding steps supplying intermediates; routes list ordered step_ids with all
dependencies. Alternatives can have separate routes. An incomplete route is valid
when its missing steps/structures remain disclosed.
For every multistep synthesis you present, populate routes with its ordered
step_ids; leaving routes empty produces individual-step views instead of the
shared route scheme. Reuse the same molecule ID for a carried intermediate and
declare the actual after_step_ids. Include the complete chosen route in a revised
answer, not only its changed step. Keep distinct alternatives in separate routes.
Preserve branches and convergent dependencies: never change chemistry or invent
intermediates to force a linear layout. The workspace draws supported linear
routes in the three-column SVG format and retains individual views otherwise.
Supply structures and attributed annotations; do not hand-draw SVG or embed image
links in answer_markdown. Keep molecule names in the structured molecule records;
the route drawing omits compound captions and numbers.
Conditions are separate attributed text fields such as solvent, temperature,
duration, quantities or addition order. Use [] if absent; yield_info=null if unreported.
Each step's reagents is a separate list of attributed claims containing only
short reagent/catalyst names, e.g. CDI or DIPEA. Route arrows show these names and
supplied partner structures. Keep solvent, quantities, temperature, time, workup
and purification in conditions, which are visible beside each step. Do not repeat a
partner already drawn structurally unless its supplied form needs clarification.
Use [] when no reagent names are explicitly supported; do not infer them from a
reaction name. Older saved answers may omit reagents.
For a retrosynthesis, give each step a concise rationale: an attributed claim
(text, basis, source_ids, limitations) explaining the bond change, chemical choice,
and what the closest inspected source reactions support. Name material substrate,
selectivity or condition differences that limit transfer. Link the relevant source
and exact example locator; never describe a distant analogue as an exact precedent.
The browser leads with the route SVG, then each step's rationale, conditions and
saved precedent SVGs with publication links. Technical diagnostics, full prose and
review notes are available in collapsed details. Avoid duplicating every step,
raw IDs, tool diagnostics or assessment counts in answer_markdown. Keep it to a
short route explanation and material unresolved questions; preserve detailed
qualifications in the appropriate structured limitations. Rationale follows the
same citation checks as conditions. Older answers may omit it; never invent missing
rationale, structures, conditions or experimental support to complete the display.
Minimal nested shapes (replace IDs/text with your actual evidence and proposal):
"steps": [{"id": "s1", "title": "Proposed step", "basis": "proposed",
           "reactant_ids": ["a"], "product_ids": ["b"],
           "precedent_refs": [],
           "condition_precedent_refs": [],
           "rationale": {"text": "A development hypothesis; evidence is incomplete", "basis": "proposed"},
           "conditions": [{"text": "Conditions to develop", "basis": "unknown"}]}],
"routes": [{"id": "r1", "title": "Proposal", "step_ids": ["s1"]}]
Condition objects use text, basis, source_ids and limitations, never label/value.
Routes use id, title, step_ids and limitations; basis belongs to steps. Keep proposed
conditions distinct from reported observations and never label proposed yields as reported.

Each step's precedent_refs must link supporting-reaction inspection artifacts for
the ACTUAL final structures, including specified stereo. If saved disconnections,
assessments or inspections contain local support, publication requires a matching
nonempty inspect_step_precedents inspection. Earlier alternatives, arbitrary source
reaction IDs and inspect_route_step artifacts do not satisfy this requirement.
Missing links are rejected with exact inspection arguments or an existing inspection
ref. Inspect those records, reconsider transfer claims and attach the right ref before
resubmitting. Do not drop steps/citations to bypass this check. Do not copy source
reaction SMILES into authored precedent cards. If no local support is available,
retain literature citations and an honest gap; never invent support or run a tool
merely to fill a panel. Empty searches and missing/uninspected evidence remain distinct.

For condition evidence, condition_precedent_refs links completed recorded
inspect_condition_precedents calls for the actual step's reactants and products,
including specified stereochemistry. The browser reads the saved observations and
their exact experimental joins; do not author replacement source cards. These links
do not replace precedent_refs or satisfy an available route-support inspection
requirement. Both kinds of evidence may support the same step.

Always supply drawable structures when inspected evidence supports them; the
browser automatically generates SVGs. For every literature-supported step, add
literature_reactions when source reactants and products can be identified without
inventing chemistry. This is separate from indexed precedent_refs and never
replaces a required local inspection. Never copy a proposed target reaction into
this field merely because the step cites a paper. Use the actual source reaction,
including an analogue's different substrates and product.
Each literature reaction uses schema_version='literature_reaction.v1', title,
source_id of a captured external source, exact locator, structure_origin,
structure_evidence (a verbatim passage from that captured source), reactants and
products (lists of {name, smiles}), conditions, yield_info and limitations.
Use structure_origin='source_explicit' only when every supplied SMILES appears in
the captured text. Otherwise use 'reconstructed_from_description': this means
agent-interpreted structures, displayed as an unverified literature reconstruction.
Source drawings must retain material-form, tautomer, stereo and interpretation
uncertainty. If the source only identifies a salt without sufficient structure,
show a supported neutral parent only with an explicit visible salt-form limitation;
never invent its protonation or counterion arrangement. If endpoints remain
ambiguous, omit the drawing and explain the gap in the step limitations.
Conditions and yield use the usual attributed claims, basis='reported' and
source_ids including this source_id; leave them absent when unreported. Proposed
target conditions stay in the step's conditions, not the source drawing.
Example (replace all structures and evidence with the inspected source):
"literature_reactions": [{"title": "Source Example 2", "source_id": "lit1",
  "locator": "Experimental Methods, Example 2",
  "structure_origin": "reconstructed_from_description",
  "structure_evidence": "Ethanol was oxidized to acetaldehyde.",
  "reactants": [{"name": "Ethanol", "smiles": "CCO"}],
  "products": [{"name": "Acetaldehyde", "smiles": "CC=O"}],
  "limitations": ["Structures reconstructed from the captured description."]}]
Use [] when no source reaction can be reconstructed. Do not hand-draw SVG or embed
image links; SVG generation is the presentation layer's responsibility.

Finish with w.finalize_answer(draft_path, draft, findings=findings). draft_path is the
exact runtime-provided answer-draft.json path for this attempt; draft is a Python dict.
The helper validates citations, records supplied self-review and saves the complete
answer. Print only its small JSON runtime handoff result and return the same receipt
as your final message. A manually written full draft follows the same contract.
The helper fills only boilerplate: omit empty arrays, null yield_info, null source URLs,
schema_version and needs_user_input=False. Supply chemistry, basis, attribution,
uncertainty and review explicitly. Invalid fields/support remain rejected. Keep answer
and findings in one small script; avoid custom builders, full draft output and full
schema reads unless validation requires them. Correct errors using saved evidence;
do not weaken validators or rewrite evidence.

For a scientific recommendation, challenge the final draft and supply brief findings.
Each finding has area, claim, assessment, evidence_refs and reason. Cover all five
areas: source_identity, structure_and_stereochemistry, conditions_and_yields,
route_completeness and counterevidence. assessment is supported, partial, unsupported,
conflicting, not_checked or not_applicable; explain missing checks. Supported, partial
and conflicting findings need actual evidence_refs. Correct overclaims and finalize
again after changes so the review matches the draft. Self-review is agent-authored,
not independent chemistry validation, and cannot be cited as scientific evidence.
Non-scientific or clarification-only answers may omit findings; no review is invented.
