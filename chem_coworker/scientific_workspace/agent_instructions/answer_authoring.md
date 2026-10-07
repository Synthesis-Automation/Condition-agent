Saved answer contract and handoff:

For synthesis routes, use w.route_answer with compact input. Supply only the target,
ordered reaction SMILES, conditions and inspected support. The application assembles
molecule/source/step IDs, exact shared intermediates, citations and the existing
scientific_answer.v2 view. The browser generates the route SVG and exposes reaction
SMILES. Do not write SVG, a custom answer builder or a second prose route.

```python
draft = w.route_answer({
    'target_smiles': target_smiles,
    'routes': [{'title': 'Proposed route', 'steps': [
        {'reaction_smiles': reactants + '>>' + product,
         'conditions': 'Conditions to develop', 'conditions_basis': 'unknown',
         'support': [{'ref': captured_source_ref, 'locator': 'Example 1',
                      'note': 'Analogue only; explain the relevant substrate difference.'}],
         'limitations': ['Retain any material unresolved chemistry.']}
    ]}]
})
receipt = w.finalize_answer(draft_path, draft)
print(receipt)
```

Steps and conditions default to proposed. Use basis='reported' or
conditions_basis='reported' only for the actual inspected experiment, with support.
Missing conditions and support may be omitted; their absence remains visible.
Optional reagents lists concise reagent names for route arrows. Do not invent
conditions, yields, atom contributors, source identities or reaction structures.
Do not transfer an analogue's yield to the proposed reaction. Keep conditions for
the proposal separate from the reported source experiment.

Support entries need only ref, locator and an optional note about relevance or
limitations. External refs must be saved literature_source/literature_excerpt
artifacts: use an exact captured passage with w.capture_source or import an existing
browser export with w.capture_source_file. A bare URL or remembered citation is not
inspected support. Title and URL are read from the saved source. A captured passage
is enough; do not copy entire pages just to produce the answer. Local refs may cite
completed scientific calls. Saved inspect_step_precedents, inspect_condition_precedents
and assess_proposed_recipe refs are attached to their corresponding step evidence
fields automatically. These references must match the actual final structures and
stereochemistry. If local support exists, publication still requires a matching
nonempty inspect_step_precedents inspection. Inspect and attach the records named
by the validation error; no new search or chemistry check is required for formatting.

Put steps in synthetic order. Unique exact intermediate links are assembled only
within each route, preserving stereochemistry and component multiplicity. For an
ambiguous producer, supply after_steps as one-based earlier step positions in that
route; [] declares an independent step. Missing steps remain gaps. No feasibility,
atom mapping, source interpretation or experimental validation is inferred.

No separate rationale, five-area self-review, reconstructed literature drawing or
full experimental procedure is required. Retain material counterevidence, uncertain
conditions and missing starting-material preparation in limitations. Detailed source
drawings and explicit review remain optional when requested or decision-relevant;
use w.help for their existing helpers. Do not run tools just to fill answer panels.

For other answers, w.answer_template(answer_markdown) creates the same canonical
answer with empty optional collections. Supply only relevant attributed objects;
reported/computed objects require source_ids. w.help('finalize_answer') exposes claim
shapes. Do not read or reproduce the full schema unless a validation error requires it.

Finish once with w.finalize_answer using this attempt's exact runtime-provided
answer-draft.json path. It validates evidence and saves the answer. Return only its
small handoff receipt; do not print the full draft. On failure, correct the indicated
fields using saved evidence and retry. w.answer_preflight is optional for debugging;
its missing-check warnings do not require new scientific calls merely to publish a
clearly qualified proposal. Do not weaken validators or rewrite saved evidence.
