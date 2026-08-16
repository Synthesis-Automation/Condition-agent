"""Retrosynthesis route integration for independent forward audits."""

from __future__ import annotations

from core_retrosynthesis import (
    MoleculeOccurrenceNode,
    PlannedRouteAction,
    ReactionRouteTree,
    RouteReactionNode,
    RouteStepEvidence,
    assess_route_tree_forward,
    build_forward_library_from_generic,
    build_generic_library,
)


def _libraries():
    generic = build_generic_library(
        (
            {
                "reaction_id": "reaction-1",
                "reference_id": "reference-1",
                "reaction_smiles": "CCBr.N>>CCN",
            },
        ),
        levels=("L1", "L2"),
        admission_mode="data_driven",
    )
    return generic, build_forward_library_from_generic(generic)


def test_planned_route_step_receives_advisory_forward_assessment() -> None:
    generic, forward = _libraries()
    operator_id = generic.templates[0].operator_id
    tree = ReactionRouteTree(
        tree_id="tree-1",
        route_kind="planned",
        target_smiles="CCN",
        root=MoleculeOccurrenceNode(
            occurrence_id="molecule-product",
            smiles="CCN",
            depth=0,
            terminal=False,
            terminal_evidence="none",
            unresolved_reason=None,
            reaction=RouteReactionNode(
                reaction_node_id="reaction-node-1",
                step_id="step-1",
                depth=1,
                reaction_smiles="CCBr.N>>CCN",
                evidence=RouteStepEvidence(evidence_kind="predicted"),
                planned_action=PlannedRouteAction(
                    operator_id=operator_id,
                    disconnection_site_key="SITE1:test",
                ),
                children=(
                    MoleculeOccurrenceNode(
                        occurrence_id="molecule-1",
                        smiles="CCBr",
                        depth=1,
                        terminal=True,
                        terminal_evidence="stock",
                        unresolved_reason=None,
                    ),
                    MoleculeOccurrenceNode(
                        occurrence_id="molecule-2",
                        smiles="N",
                        depth=1,
                        terminal=True,
                        terminal_evidence="stock",
                        unresolved_reason=None,
                    ),
                ),
            ),
        ),
        reaction_count=1,
        maximum_depth=1,
        fingerprint_tokens=(),
    )

    report = assess_route_tree_forward(tree, forward)

    assert report.route_ranking_impact == "none_advisory_only"
    assert report.disposition_counts == {"clear": 1}
    assert report.high_risk_step_ids == ()
    step = report.step_assessments[0]
    assert step.step_id == "step-1"
    assert step.assessment.targeted_replay_status == "structurally_reproduced"
    assert step.assessment.intended_product_rank == 1
    assert step.assessment.validity == "structurally_supported"
    checks = {check.check_id: check.status for check in step.assessment.checks}
    assert checks["targeted_operator_replay"] == "pass"
    assert checks["reaction_signature"] == "pass"
    assert checks["reverse_precursor_recovery"] == "pass"
    assert checks["operator_edit_agreement"] == "pass"
    assert checks["condition_compatibility"] == "not_evaluated"
