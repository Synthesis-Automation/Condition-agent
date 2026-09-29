# Retrosynthesis investigation adviser v4

This is an optional menu, not a required workflow. Skip, reorder, repeat, or replace
suggestions according to the question and available evidence. Scientific contracts,
recorded provenance, and honest uncertainty still apply.

Questions worth resolving:

- Is the target's connectivity, chemical form and specified stereochemistry correct?
- Is there an exact published route, a useful analogue, or an unsupported hypothesis?
- Which intermediate deserves another expansion, and what makes a terminal material
  acceptable under the user's constraints? Vendor listings do not confirm stock.

Useful recorded calls (with an open workspace `w` and the chosen target):

```python
audit = w.run("analyze_molecule", {"smiles": target_smiles})
cut = w.run("disconnect_target", {"target_smiles": target_smiles, "top_k": 3})
```

Reuse an existing audit for the same structure. Keep executable script bodies under
`if __name__ == "__main__":` so importing a helper does not repeat scientific calls.
Print `w.call_summary(event)` first; inspect selected fields only when needed.

When an unfamiliar core is the synthesis bottleneck, the following optional
pattern can help. An exact supported route or a clear disconnection may already
answer the question; skip fragment discovery when it adds no useful evidence.

- State the uncertainty first: which core needs a construction precedent, and why?
  Explain the chemical question in a short progress message. A complex-looking
  fragment or high structural score is not automatically the difficult part.
- Select your own connected core, or use `suggest_search_fragments(target_smiles=...)`
  for graph-validated regions. Candidates overlap and are not precursors. Inspect
  atom boundaries and cautions; keep distinctive ring topology, junctions and
  heteroatom positions. Avoid an unrestricted common biphenyl query when the
  distinctive chemistry lies elsewhere. No suggestion call needs an index or retro.
- Usually start with one or two informative queries using
  `search_fragment_precedents(query=..., target_smiles=target, limit=5)` when
  `fragment_index` is configured. Supply the target for every target-derived core,
  including broadened queries: the exact search semantics must still match it.
  A mismatch fails before scanning and is not evidence of zero precedents.
  This is an advisory starting budget, not a mandatory sequence or hard search cap.
  Search another region only to resolve a remaining question. Reuse unchanged
  queries and saved results. `too_broad` suggests adding distinctive context;
  `partial` has incomplete coverage; zero hits apply only to this index.
  Refinement with explicit SMARTS must state what was relaxed. There is no general
  graph-edit similarity search, and altered ring topology is not an exact-core hit.
- If a whole-target query returns only carried-through hits or no construction
  evidence, consider one core query with peripheral substituents removed. Preserve
  distinctive ring topology and heteroatom positions; explain the relaxation.
  This is useful when the bottleneck is building that core, not a reason to search
  generic rings or every suggested fragment. Skip it when it cannot change the plan.
- Inspect the selected hit's local bond-change witnesses, observation/reference ID,
  exact procedure link and source scope. A product containing the core is not
  evidence that this reaction constructed it. Distinguish construction, modification,
  boundary changes, retention and unresolved correspondence; groups can overlap.
  A construction witness can establish only part of a core, not necessarily its
  entire synthesis. Follow the original source when needed; preserve missing data.
- Explain transfer to this target: which observed bonds/operation could be reused,
  which substituents, reactive handles or stereochemistry differ, and what remains
  a hypothesis. Keep reported chemistry separate from proposed adaptations.
  Use single-step retrosynthesis for gaps when helpful, then assess explicit steps
  and route topology through the existing operations. A retrieved hit is not an
  automatically validated route or a reason to run forward prediction.

Record consequential choices with `w.store.note('decision', text, evidence_refs=(ref,))`.
Before describing a chosen route step as precedent-supported, call
`inspect_step_precedents(source_ref=..., realization_id=...)` on its saved
disconnection, or select `step_id` from a saved route assessment/revision. Inspect
source reactions, comparisons and caveats; follow `page.next_offset` if needed.
Include the returned inspection refs in that answer step's `precedent_refs` so the
reader sees the same evidence. Template support and high similarity do not establish
experimental feasibility. Distinguish uninspected steps from bounded searches that
retrieved no supporting reactions. Literature-only evidence keeps its captured
citations; do not invent template support or run an inspection merely to fill a panel.
When fragment evidence affects the answer, show a concise comparison of selected
core, inspected construction evidence, proposed transfer and unresolved gap. Cite
specific records/fields, not just the search call as blanket support. If no useful
precedent emerges, report the bounded search and continue with another hypothesis.

The agent owns multi-step planning. Select concrete single-step realizations and
choose the next precursor yourself. Avoid cycles and repeated expansions; record
branch choices, alternatives, constraints and stopping reasons. An exact literature
route may make further searches unnecessary. Do not invoke the built-in multistep
planner, including through custom scripts.
Expand an intermediate only to resolve a decision-changing gap. Do not repeat a
disconnection to regenerate an already explicit, source-supported step.
Use `w.store.note("decision", text, evidence_refs=(ref,))` for branch choices. Other
valid note kinds are `hypothesis`, `question`, `limitation`, and `review`.

Assemble steps with `external_step_id`, `target_smiles` and dot-separated
`precursor_smiles`. Use `assess_route_step` or `assess_route_proposal`; extend or
replace explicit branches with `revise_route_branch`, then compare alternatives
when useful. Start with the standard structural assessment. Inspect its weak steps
before adding optional computation. Unsupported reconstruction can reflect library coverage, not chemical
impossibility. Never fabricate atom maps, atom donors, or structural correspondence.
For bromination or oxidation, inspect the reported bromine or oxygen source and represent known contributing
reactants accurately; if the donor is missing or ambiguous, retain that limitation
instead of inventing a reagent to satisfy an atom-source check.

Forward prediction is useful only when competing products could change a route
decision. For that question, optionally call `assess_route_step_forward` with a saved
route `source_ref`, one eligible `step_id`, the decision-changing `question`, and
`timeout_seconds=30` (the maximum). It loads the prebuilt, baseline-pinned
`forward_library`; it never rebuilds libraries or checks the whole route. The worker
records its stages and stops at the deadline. No process polling is needed. Missing
libraries, timeout or failure leave the question unresolved. Skip this check when the
main gap is missing atom contributors, source evidence, or starting-material supply.
This optional challenge is separate from the structural validation in single-step
retrosynthesis and standard route assessment; those checks remain in place.

If direct retrieval fails, inspect the primary source in the browser and preserve
actual visible text with `w.capture_source(text, url=url, locator=locator)`. Record
an exact passage with `w.record_source_excerpt`. Pass its artifact reference in
`evidence_refs` as provenance; citations do not override graph or topology gates.
Failed downloads remain debugging records, not supporting literature.
After an environment-wide network denial, use browser capture for subsequent sources
as well; changing URLs does not justify repeating the same denied download path.
Usually try one alternative access path for an unavailable paper or supplement;
then switch to accessible primary evidence or state the gap. Continue only with
a specific reason another search could change the conclusion. Repeated title/DOI
variants without new evidence are a signal to stop, not a fixed search quota.

Stop with a supported route or a clearly bounded partial proposal. Identify unresolved
steps/leaves, availability assumptions, condition gaps and the next useful check.
Structural admission is not experimental feasibility; a claimed published yield or
stereoisomer needs an inspected source that actually supports it.

Show the route as explicit reaction steps with concise reagent/catalyst labels,
then explain the choice, strongest source and main uncertainty in one short
paragraph of two or three sentences (usually 60-90 words). Do not repeat the long
target name, list the steps again, or append a catalogue of caveats. Detailed procedures and extended analysis wait
for a user request. Preserve the full evidence and review for that follow-up;
avoid duplicating the route in prose, tables or claims.
Use `w.finalize_answer(draft_path, draft, findings=findings)` to fill empty-field
boilerplate, validate citations, save the explicit self-review and write the full
answer. The agent still supplies chemistry, attribution and review findings.
Before finalization, attach `inspect_step_precedents` artifacts for the final
step structures whenever saved calls contain support. Earlier alternatives and
`inspect_route_step` calls do not satisfy this requirement. A rejection names the
step and supplies inspection arguments or an existing inspection ref; inspect the
records, reconsider transfer claims, and attach the ref before submitting again.
