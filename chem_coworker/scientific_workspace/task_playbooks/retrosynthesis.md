# Retrosynthesis investigation adviser v8

This is an optional menu. Choose the questions that matter to the user's decision;
skip, reorder, repeat or replace suggestions. Shared scientific contracts still apply.

## Purpose and context

Help the user assess a route, a useful next disconnection, or a bounded partial proposal.
Establish the target's connectivity, chemical form, protecting groups and specified
stereochemistry. Consider material, condition and scope constraints. Define what makes
a terminal material acceptable; a vendor listing does not confirm stock.

## Scientific questions

- Is there an exact inspected route, a transferable analogue, or only a hypothesis?
- Which construction or unresolved step is the current synthesis bottleneck?
- Which concrete precursor branch deserves expansion, and could it change the plan?
- Which steps, terminal materials or stereochemical requirements lack support?

Use `analyze_molecule`, `disconnect_target`, saved step inspections and
`assess_route_step` or `assess_route_proposal` as useful. The agent selects realizations,
subsequent targets and stopping decisions under the catalogue's execution policies.
Choose an initial investigation budget suited to the question; broaden it for a
decision-changing gap. Reuse explicit source-supported steps and unchanged results.
Keep branch choices, alternatives and stopping reasons in recorded notes; avoid cycles.

For a decision-changing bond-construction question, `disconnect_target` accepts
`required_disconnection_bond=[atom_a, atom_b]` with `focus_target_smiles` from a
saved `inspect_reactive_sites` or `compare_molecules` canonical molecule. IDs are
zero-based in that canonical SMILES, not in the input string or atom-map labels.
Inspect the atom inventory before selecting a bond. The bond must be absent in
the proposed precursors and verified as formed in the forward transformation;
ring closure is allowed, but bond-order changes alone do not qualify. Focus is
checked before validation-budget selection and again against final observations.
Inspect `bond_focus`, candidate `bond_focus_check`, and `search_diagnostics` for
scope, rejection counts and budget exclusions. An empty focused search must not
be presented as impossibility or silently broadened. Broaden explicitly only
when another question could change the plan. Template retrieval is unchanged;
the ladder keeps per-level budgets and may try more levels when focus is narrow.

Use `assess_starting_material` when deciding whether to expand a route leaf. It
checks exact registry identity, then exact product occurrence in the optional
fragment index, then a configurable molecular-weight cutoff (default <200 g/mol).
Pass the user's `unavailable_starting_materials`; exclusions override all stops.
Registry membership permits assumed obtainability, literature occurrence needs
preparation inspection, and the MW fallback is only a planning heuristic. Each
stage can be disabled. Preserve the call ref, stop reason, failed lookups and
availability assumptions in notes and the final answer. A stopped branch is not
verified supply; literature/MW leaves leave the route partially resolved. This
operation does not upgrade route admission or automatically expand any branch.

Read supporting advice only when relevant:

- `w.task_guide('retrosynthesis_fragments')`: construction precedents for a chosen core.
- `w.task_guide('retrosynthesis_revision')`: change an explicit route branch or constraints.
- `w.task_guide('retrosynthesis_forward')`: challenge a consequential product-competition question.

## Evidence distinctions

Keep exact reported routes, analogue transfers and proposed steps distinct. Compare
actual source structures, forms and stereochemistry. Templates, similarity and structural
admission do not establish experimental feasibility. A fragment construction witness
may support only part of a core. Cite inspected records for the actual chosen steps;
the shared answer contract owns their required evidence links.

Distinguish uninspected steps from bounded searches with no supporting reactions.
Missing atom contributors, ambiguous correspondence, condition gaps and material
availability remain limitations. Unsupported reconstruction can reflect library coverage.
Never invent donors or mappings to pass a check. Literature-only support retains its
captured citations; optional computation need not accompany an inspected source route.

## Stopping criteria

Stop with a supported route or a clearly bounded partial proposal. Identify unresolved
steps and leaves, availability assumptions, condition gaps and the next useful check.
Further expansion is useful when it can change the user's decision. Shared answer
and presentation instructions own publication requirements and formatting.
