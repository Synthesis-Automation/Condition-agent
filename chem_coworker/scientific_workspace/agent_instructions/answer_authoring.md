Saved answer contract and handoff:

For synthesis routes, spend the investigation on route generation and precedent
inspection. The answer needs only route SVGs, reaction SMILES and precedent support.
Use w.route_answer once; the application assembles IDs, shared intermediates,
citations and scientific_answer.v2, and draws the SVGs. Do not write SVG, a custom
builder or diagnostic-panel text. A short summary may explain the preferred route,
its shared bottleneck and why an alternative changes that risk; do not duplicate
the full route in prose. Step-by-step rationales remain optional.

```python
draft = w.route_answer({
    'target_smiles': target_smiles,
    'summary': 'Preferred proposal and its main unresolved bottleneck.',
    'routes': [{'title': 'Top pick', 'steps': [
        {'reaction_smiles': reactants + '>>' + product,
         'support': [{'ref': saved_ref, 'locator': 'Example 1'}]}
    ]}]
})
print(w.finalize_answer(draft_path, draft))
```

Steps default to proposed. Keep missing support missing. Do not invent structures,
conditions, source identities or yields. Use reported basis only for the actual
inspected experiment. Preserve material gaps or conflicts in limitations, briefly;
use route/global limitations for a shared unresolved preparation or supply gap.
The application displays saved admission failures without agent-authored diagnostics.
If a later candidate addresses an unresolved step, reuse its assessment_proposal
(including saved_candidate) in an explicit step reassessment or route revision.
Inspecting a precedent alone does not update an earlier assessment. Retain remaining
experimental and stereochemical gaps; do not invent a mapping to clear a warning.

Support refs must be saved captured literature/excerpts or completed scientific
calls. Capture exact inspected passages with w.capture_source or
w.capture_source_file; a bare URL is insufficient. Titles and URLs come from the
saved source. A support note is optional. Previously recorded nonempty inspections
attach to proposed steps with exactly matching structures including stereo. Explicit
support refs remain authoritative; recipe checks attach only when supplied. If local support
exists, publication still requires a matching
nonempty inspect_step_precedents inspection. Reuse saved matching evidence.

Put steps in synthetic order. Unique exact intermediate links are automatic;
after_steps optionally names earlier one-based positions for ambiguous producers.
Missing steps remain gaps. No feasibility or atom correspondence is inferred.

Conditions are optional: omit them unless requested or already useful to the route.
When recommending a concrete recipe, assess_proposed_recipe checks its actual
components and operating conditions; attach that saved ref in support and retain
unresolved coverage. Structural analogues alone do not validate a recipe. Source
conditions and yields belong to their source experiments, not to the proposed target.
No separate rationale, five-area self-review, literature reconstruction or full
procedure is required. Do not run scientific tools merely to populate answer fields.

Use w.help('route_answer') only if needed. Do not read answer-schema.json at startup.
Other tasks can use w.answer_template(answer_markdown). Finalize once to
the exact runtime-provided answer-draft.json path and return only the handoff receipt.
On validation failure, fix the indicated fields using saved evidence; do not weaken
validators or rewrite evidence. w.answer_preflight is optional for debugging.
