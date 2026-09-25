# Condition investigation and recipe assessment migration

Date: 2026-09-25

This development change improves direct condition investigations without changing
dataset admission, index formats, retrieval ordering or compatibility definitions.
Independent chemist review and the roadmap's release gates remain pending.

## Reaction-level assessment contract

`assess_reaction_recipe()` now returns the immutable `ReactionRecipeAssessment`,
with `schema_version="reaction_recipe_assessment.v2"`. Existing compatibility
fields remain available, with additional `analysis_status` and `analysis_warnings`.

Previously, missing structural signatures produced `status="conflict"` and the
hard-conflict code `UNRESOLVED_REACTION_FOR_RECIPE_ASSESSMENT`. They now produce
`status="unknown"`, an empty hard-conflict list and
`unresolved_requirements=["VERIFIED_REACTION_SIGNATURE_REQUIRED"]`. Invalid input
is separately labeled `invalid_input`. `compatible=False` in these cases means
assessment cannot establish admission; consumers must inspect status and must
not equate this Boolean with demonstrated chemical incompatibility.

Supported reactions retain the canonical rule engine's `conflict`, `unknown`, or
`no_known_conflict` results. Direct assessment now supplies unchanged spectator
groups from the full reaction analysis to that engine: the identity signature
alone does not retain all of this context. This fixes missed spectator-condition
conflicts. Regression cases cover supported reactions, unresolved transformations,
invalid notation, conflicting mapping and an unchanged aldehyde with an oxidant.
No new compatibility rule or definition version is introduced; the existing
`compatibility.v1` definition schema remains `1.3`.

## Inspection and proposal boundaries

`condition_recommender.compare_condition_evidence()` exposes selected indexed
observations, full indexed recipes, bond-edit/environment multiset differences,
canonical compatibility results, missing operating fields and distinct-reference
counts. It reuses one query analysis and does not rank recipes or certify transfer.
Counts are for the selected observations, not the whole corpus; distinct citation
identities do not establish independent experiments or patent-family independence.

The scientific workspace's `inspect_condition_precedents` composes this function
with existing pagination and procedure access. Exact observation IDs link records;
reaction-only procedures remain explicitly unassigned. All matching procedure
records and their source metadata remain inspectable, including missing text.
No operating values are borrowed from a different observation. Missing indexed
data does not establish absence from the original publication.

`propose_condition_adaptation` records a complete proposed recipe normalized by
`condition_registry`, its original inspected recipe, every changed field and
reason, cited evidence, assumptions, risks and canonical compatibility. It always
remains an unreviewed agent proposal with `transfer_status="not_established"`.
The service verifies evidence existence, not whether it scientifically supports
the change. Staged protocols and declared absences are currently rejected for
adaptation rather than silently dropped. Nothing is admitted to source datasets.

## Execution and answer correction

`ScientificWorkspace.run_python()` snapshots a local script and JSON inputs,
including cited artifacts, and records its interpreter, baseline identity, output,
exit status and log hashes. Only successful recorded executions may support the
answer's `computed` label. A mere file attachment still does not establish execution.
The runner has a deadline and output-size checks. It is a trusted local process,
not OS isolation: the script has the investigator's permissions and must declare
external inputs and respect the fixed baseline. Arbitrary scripts are never replayed
automatically. Correctness, complete dependency declaration and chemical validity
are not inferred from exit status.

After a schema/evidence-reference failure, the conversation service permits one
correction using the same runtime thread. Both rejected drafts and correction
traces remain saved, and usage is retained for each attempt. A second invalid
answer fails visibly. Runtime failures, cancellation and baseline drift are not
repaired. The configured runtime deadline applies per attempt, so a corrected
turn can use up to two runtime deadlines plus preparation overhead.

No API route changes are required. Start a new investigation after updating code;
historical fixed-baseline records remain readable. See the
[workspace quickstart](../AI-native/Scientific_Workspace_Quickstart.md) and
[implementation status](../AI-native/Scientific_Workspace_Implementation_Status.md).
