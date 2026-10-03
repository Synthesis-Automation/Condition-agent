# Precedent-backed retrosynthesis validity tool

Status: development implementation, pending independent chemistry review and
untouched evaluation. Existing release gates remain open.

`assess_retro_validity` assesses one concrete precursor-to-target realization.
It combines existing external-proposal structural gates, graph-qualified source
precedents, compatibility assessments, and an optional separately recorded forward
audit. It supplies advisory actions to the agent without changing canonical
disconnection ranking, route admission, or public actionability.

## Evidence grades

The versioned policy is `core_retrosynthesis/definitions/retro_validity.v1.json`.
The rank is ordinal evidence strength, not reaction success probability or yield.

| Grade | Rank | Meaning |
| --- | ---: | --- |
| `exact_reaction` | 4 | Canonical whole-reaction identity and qualified edit graph agree |
| `detailed_local_match` | 3 | Detailed local graph agrees; condition retrieval L0 |
| `close_transformation_analogue` | 2 | Shared transformation with generalized source/departing ports; L1 |
| `broad_core_analogue` | 1 | Broader qualified core with generalized local environment; L2 |
| `operator_only` | 0 | Exact admitted operator support without a qualified precedent |
| `unresolved` | null | No sufficient evidence for an ordinal grade |

Whole-reaction identity retains specified stereochemistry and component
multiplicity while ignoring atom-map numbers and component serialization order.
Fingerprint equality, family names and product identity alone cannot establish
an exact grade. Retro-template specificity levels use a separate namespace.

`condition_recommender.reaction_precedent_support` reads the canonical condition
index and its bound shared-core projections. It uses the same
`compare_reaction_cores` qualification as condition recommendation, without
recipe aggregation, recipe ranking, or a new conversion/index path. Conditions
and outcomes remain attributed to individual source observations. A precedent
recipe is not silently applied to the proposal.

The operator-library match set comes from the existing canonical external-step
assessor. Operator support and optional corpus support retain separate scoped
counts; they are not added together as independent evidence. Repeated templates
do not become additional publications. Distinct references are not proof of
independent replication. Missing reference IDs contribute no reference count.
Both lookup and returned matches are bounded and disclose truncation.

When no mapping is supplied, corpus lookup uses the original structure-only
query, as condition retrieval does. Operator support uses the assessor's validated
materialized mapping. Partial internal mapping can yield different projections;
both queries are retained and the report requires review when they differ.
Supplied mappings remain authoritative and are never discarded to improve support.

## Decision contract

`RetroValidityAssessment` is immutable and serializes as
`retro_validity_assessment.v1`. The workspace wrapper is
`retro_validity_investigation.v1`; source indexes and existing assessment schemas
are unchanged and require no rebuild.

Inspect these axes separately:

- `structural_status`: verified, unresolved, contradicted, or invalid input.
- `precedent_grade` and `evidence_rank`: strongest qualified support.
- `operator_precedent_support` and `corpus_precedent_support`: source reactions,
  differences, available outcomes/recipes, counts, exclusions and coverage.
- `recipe_assessment`: canonical supplied-recipe compatibility, or null.
- `forward_status` and `forward_execution_status`: a separately bound challenge,
  missing check, or incomplete execution. Check traces and policy versions remain.
- `cautions`, `unresolved_checks`, `warnings` and `suggested_action`.

Invalid structures require correction. Verified contradictions, hard recipe
conflicts or a contradicted forward audit suggest rejecting or revising the
realization. Analogue differences and compatibility warnings suggest investigating
cautions. Operator-only and unresolved support suggest seeking evidence. Strong
support suggests inspecting the supporting source evidence.

`supported` means precedent-backed structural support. Missing conditions,
unfinished forward checks, incomplete scope rules and experimental feasibility
remain explicit even with rank 4. Evidence rank is preserved when another axis
contradicts the proposal; the overall status/action carries that contradiction.
Missing forward recovery or missing libraries never become chemical rejection.

## Agent and CLI usage

In an initialized workspace with a pinned `retro_library`:

```python
event = workspace.run("assess_retro_validity", {
    "proposal": {"target_smiles": "CCN", "precursor_smiles": "CC=O.N"},
})
print(workspace.call_summary(event))
```

Alternatively supply `source_ref` from `disconnect_target` plus `realization_id`;
from a route assessment/revision plus `step_id`; or from `assess_route_step` or a
previous validity call with neither selector. Supply exactly one of proposal or
source_ref. Each alternative realization requires its own assessment.

Pinned `condition_index` and `shared_core_index` add corpus support automatically
when both are present. An incomplete pair is disclosed; a stale or mismatched
present pair fails explicitly. No indexes or forward libraries are built by the tool.

`forward_ref` can join an existing `assess_route_step_forward` call on an eligible
saved route step. Both molecular sides, specified stereo and the supplied recipe
must match. Completed audits retain competition evidence; timeout/error/cancelled
checks retain their incomplete state. The tool does not run forward prediction.

The existing recorded CLI accepts a JSON arguments file:

```powershell
python -m chem_coworker.scientific_workspace run results/ai_native/MY_RUN assess_retro_validity --input arguments.json
```

Use `inspect_step_precedents` on the saved source or validity result to inspect
operator precedents. Use `inspect_condition_precedents` with the retrieved corpus
reaction IDs to inspect indexed source observations and procedures. Conditions,
stereo, substrate differences and outcome reporting require source inspection.

The task adviser now recommends this assessment during branch selection. It
remains the agent's decision whether to expand, revise, seek evidence or abstain.
No autonomous planner or additional model is introduced.

Restart the scientific workspace server and start a new investigation after this
scientific code change. Saved investigations remain readable. Development tests
do not establish calibrated feasibility, independent chemist acceptance, or an
untouched-evaluation pass.
