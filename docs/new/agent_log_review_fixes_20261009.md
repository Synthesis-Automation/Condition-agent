# Workspace route evidence and tool-call corrections

The October 9 workspace review compared two fresh investigations of the same
stereodefined fused-ring target. Their baselines and five initial ranked
disconnections were identical. The agents explored different branches. These
changes preserve that freedom; they introduce no route-choice, search-order,
ranking, or stopping-policy requirement.

## Scaffold query stereochemistry

`precedent_discovery.v1@1.1` preserves hydrogen counts at specified tetrahedral
centers. The virtual hydrogen participates in neighbor ordering: omitting its
SMARTS predicate could invert a ring-junction query during serialization, causing
it to reject its own target. Other hydrogen counts remain unconstrained.
Generated queries still must match the exact target with chirality enabled.
Query construction copies cached patterns before modifying atoms.

This changes generated query semantics and versioned discovery results, not the
stored fragment index. Existing indexes do not need rebuilding. Tests include
the reviewed target, enantiomer/stereocenter inversions, a quaternary center,
serialization invariance, explicit atom correspondence, and cache integrity.
This is a chemistry correction, not an independent review or release-gate result.

## Saved assessments in route presentation

Exact matching structural assessments are recovered from the investigation even
when the answer cites only a precedent inspection. Both complete reaction sides,
stereochemistry and chemical form must match. Conflicting and ambiguous receipts
remain visible alongside successful checks. Recovered receipts are labeled
`saved_exact_step_check`; they do not become claimed experimental support.

Historical conversations are bounded by their recorded answer event, so later
checks do not alter the evidence shown for an earlier answer. Existing saved
answers are not rewritten. Recipe assessments still require an explicit link:
matching molecular graphs alone cannot identify the proposed recipe.

## Optional route-input source searches

`inspect_route_inputs` contract v2 accepts omitted `leaf_queries`, omitted `terms`,
and `terms=[]`. Every route leaf is still assessed. A leaf without terms retains
`source_search_status=terms_not_supplied` and no source-search result; this is
not a negative literature result or verified availability.

Malformed terms and unknown fields are rejected before starting-material checks.
Errors identify `leaf_queries[index]`. `w.help('inspect_route_inputs')` exposes the
actual nested model, constraints, and a valid example. Nonempty searches remain
limited to 1–10 nonblank literal terms of at most 500 characters.

## Carrying selected reconstruction into assessment

Step proposals for `assess_route_step`, `assess_route_proposal` and
`revise_route_branch` may include an explicit saved candidate selection:

```python
step = {
    "external_step_id": "epoxidation",  # Route steps only.
    "target_smiles": candidate["target_smiles"],
    "precursor_smiles": candidate["precursor_smiles"],
    "saved_candidate": {
        "source_ref": disconnection_ref,
        "strategy_id": strategy["strategy_id"],
        "realization_id": candidate["realization_id"],
    },
}
```

The adapter retrieves the selected reconstruction and delegates validation to the
existing canonical proposal assessor. Changed graphs or stereochemistry fail
selection. Missing mappings remain missing; ambiguous selections and conflicting
saved mappings are rejected. Do not also supply `mapped_reaction_smiles`.

The result and event retain the source reference and `mapping_evidence`, identifying
the origin as `saved_operator_reconstruction`, with `observed_reaction=False`.
An invalid saved map fails the existing structural gates even if its saved
candidate has a successful validation label. Route revision retains provenance.
`assess_retro_validity` also retains the reconstruction when selecting a saved
disconnection, rather than discarding it before assessment.

Operation contracts are v2 for step/route assessment and route revision, and v3
for retro validity. Public domain reaction and answer schemas are unchanged.
The epoxidation from the reviewed log has a regression showing that its unmapped
proposal has unresolved correspondence while the saved reconstruction passes
that gate; neither result establishes experimental feasibility.

## Rollout

Restart the scientific server and refresh the browser to recover diagnostic
notices for existing answers. Start a new investigation for new scientific calls
because existing investigations pin their code and definition baseline. The
scientific UI requires no frontend build. There are no HTTP route changes,
dataset snapshot updates, corpus conversions or index rebuilds.

## Validation

`pytest -q`: 2,913 passed, 4 skipped (781.44 seconds). The saved-answer projection
was also checked against both original investigations without modifying their
artifacts: the earlier imine step and later epoxidation step now expose their
recorded unresolved correspondence. All four later-run steps recover their saved
assessments. These diagnostics do not resolve the original chemistry gaps.
