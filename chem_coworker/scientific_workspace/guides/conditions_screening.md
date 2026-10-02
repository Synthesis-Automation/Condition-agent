# Condition screening adviser v1

Optional supporting advice for `conditions`. Use it when the user needs diverse
experiments to screen rather than a single best-supported recipe.

## Purpose and context

Establish the intended transformation, requested panel size and experimental constraints.
Decide whether weak-label evidence is suitable for the user's screening decision.
`generate_weak_label_screening_array` supplies distinct intact recipes; inspect its
catalogue entry for arguments and dataset prerequisites.

## Scientific questions

- Which compatible recipes provide useful diversity within the user's constraints?
- Are there enough recipes, and what evidence or identity uncertainties limit the panel?
- Does the user need a summary, the complete array, or an export for experiment planning?

## Evidence distinctions

Source reaction structures in this dataset are unverified. Graph-derived checks on
the query and recipe compatibility do not verify those source reactions. Historical
yields do not predict the query's yield. Keep complete recipes, recipe IDs, source row
numbers, ambiguous source types, missing values and warnings. A source-type hint may
narrow retrieval only when consistent with structural evidence.

Weak-label row numbers belong to this dataset; `get_precedents` and `get_procedures`
query the structural corpus. Do not pass those row numbers to these operations.
For full arrays or JSON/CSV exports, use the complete saved result, not the summary
preview. Save requested exports inside the investigation and attach them with the
generating call as evidence, preserving provenance and missing values.

## Stopping criteria

Return the available compatible panel with its evidence scope and unresolved risks.
Report a shortfall rather than inventing experiments to reach the requested size.
Screening suggestions remain experiments to evaluate, not predicted successful conditions.
