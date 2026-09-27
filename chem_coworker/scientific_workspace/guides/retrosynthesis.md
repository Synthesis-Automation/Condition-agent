# Retrosynthesis investigation adviser v1

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

The agent owns multi-step planning. Select concrete single-step realizations and
choose the next precursor yourself. Avoid cycles and repeated expansions; record
branch choices, alternatives, constraints and stopping reasons. An exact literature
route may make further searches unnecessary. Do not invoke the built-in multistep
planner, including through custom scripts.

Assemble steps with `external_step_id`, `target_smiles` and dot-separated
`precursor_smiles`. Use `assess_route_step` or `assess_route_proposal`; extend or
replace explicit branches with `revise_route_branch`, then compare alternatives
when useful. Unsupported reconstruction can reflect library coverage, not chemical
impossibility. Never fabricate atom maps, atom donors, or structural correspondence.
For bromination, inspect the reported bromine source and represent known contributing
reactants accurately; if the donor is missing or ambiguous, retain that limitation
instead of inventing a reagent to satisfy an atom-source check.

If direct retrieval fails, inspect the primary source in the browser and preserve
actual visible text with `w.capture_source(text, url=url, locator=locator)`. Record
an exact passage with `w.record_source_excerpt`. Pass its artifact reference in
`evidence_refs` as provenance; citations do not override graph or topology gates.
Failed downloads remain debugging records, not supporting literature.

Stop with a supported route or a clearly bounded partial proposal. Identify unresolved
steps/leaves, availability assumptions, condition gaps and the next useful check.
Structural admission is not experimental feasibility; a claimed published yield or
stereoisomer needs an inspected source that actually supports it.
