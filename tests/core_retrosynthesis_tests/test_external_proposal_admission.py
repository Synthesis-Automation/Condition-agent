"""External proposal admission and complete-route topology regressions."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from reactive_taxonomy import featurize_reaction

from core_retrosynthesis.external_proposal_assessment import (
    ExternalRetrosynthesisProposal,
    assess_external_retrosynthesis_proposal,
    load_external_proposal_admission_policy,
    validate_external_proposal_admission_policy,
)
from core_retrosynthesis.external_proposal_review import (
    render_external_route_assessment_html,
)
from core_retrosynthesis.external_route_admission import (
    ExternalRouteProposal,
    assess_external_route_proposal,
    external_route_proposal_from_tree,
)
from core_retrosynthesis.forward_assessment import (
    build_forward_library_from_generic,
)
from core_retrosynthesis.generic_compiler import analyze_generic_reaction
from core_retrosynthesis.generic_library import build_generic_library
from core_retrosynthesis.generic_search import disconnect_generic_target
from core_retrosynthesis.mapping import materialize_atom_mapping
from core_retrosynthesis.multistep import (
    RetrosynthesisRouteStep,
    StartingMaterialAssessment,
)
from core_retrosynthesis.route_contract import validate_route_tree
from core_retrosynthesis.route_tree import build_canonical_route_tree


FIRST_REACTION = "CC=O.N>>CCN"
SECOND_REACTION = "O=C(O)c1ccccc1.NCC>>O=C(NCC)c1ccccc1"
TARGET = "O=C(NCC)c1ccccc1"
SUPPLIED_FIRST_MAPPING = (
    "[CH3:1][CH:2]=[O:3].[NH3:4]>>[CH3:1][CH2:2][NH2:4]"
)


def _row(reaction_smiles: str, ordinal: int) -> dict:
    analysis = featurize_reaction(reaction_smiles)
    assert analysis.valid and analysis.reaction_signature is not None
    value = analysis.to_dict()
    return {
        "reaction_id": f"reaction-{ordinal}",
        "reference_id": f"reference-{ordinal}",
        "reaction_smiles": reaction_smiles,
        "reaction_core": value["reaction_core"],
        "reaction_observation": value["observation"],
        "reaction_completeness": value["reaction_completeness"],
    }


@pytest.fixture(scope="module")
def operator_library():
    return build_generic_library(
        (_row(FIRST_REACTION, 1), _row(SECOND_REACTION, 2)),
        levels=("L0", "L1", "L2"),
    )


def _route_value() -> dict:
    return {
        "target_smiles": TARGET,
        "external_route_id": "web-route-claimed-id",
        "sources": [
            {
                "source_id": "chemist-sketch-1",
                "source_kind": "chemist",
                "quoted_claim": "two-step illustrative route",
            }
        ],
        "steps": [
            {
                "external_step_id": "step-2",
                "target_smiles": TARGET,
                "precursor_smiles": "O=C(O)c1ccccc1.NCC",
                "claimed_transformation": "an externally claimed named reaction",
            },
            {
                "external_step_id": "step-1",
                "target_smiles": "NCC",
                "precursor_smiles": "N.CC=O",
            },
        ],
    }


def test_policy_fixes_gate_order_and_review_only_boundary() -> None:
    policy = load_external_proposal_admission_policy()

    assert len(policy.required_gate_order) == 14
    assert policy.required_gate_order[0] == "input_structure"
    assert policy.required_gate_order[-1] == "condition_support"
    assert policy.ranking_influence == "none_review_only"
    assert policy.public_actionability == "disabled_until_release_gate"

    invalid = json.loads(
        Path(
            "core_retrosynthesis/definitions/external_proposal_admission_policy.v1.json"
        ).read_text(encoding="utf-8")
    )
    invalid["admission"]["ranking_influence"] = "rerank"
    with pytest.raises(ValueError, match="review-only"):
        validate_external_proposal_admission_policy(invalid)


def test_unmapped_step_recovers_internal_operator_and_precedent(
    operator_library,
) -> None:
    proposal = ExternalRetrosynthesisProposal(
        target_smiles="CCN",
        precursor_smiles="N.CC=O",
        claimed_transformation="fabricated name that must not route admission",
        external_proposal_id="untrusted-operator-like-id",
    )

    assessment = assess_external_retrosynthesis_proposal(
        proposal, operator_library
    )
    materialized = materialize_atom_mapping(FIRST_REACTION)
    assert materialized is not None
    internal_identity = analyze_generic_reaction(materialized.reaction_smiles)

    assert assessment.status == "precedent_supported"
    assert assessment.strongest_evidence_tier == "precedent_supported"
    assert assessment.admission_eligible
    assert not assessment.actionable
    assert assessment.internal_confidence_parity
    assert assessment.reaction_identity is not None
    assert internal_identity is not None
    assert assessment.reaction_identity.operator_signature == (
        internal_identity.operator_signature
    )
    assert len(assessment.precedent_matches) == 1
    assert assessment.precedent_matches[0].reaction_id == "reaction-1"
    assert tuple(gate.gate_id for gate in assessment.gates) == (
        load_external_proposal_admission_policy().required_gate_order
    )


def test_source_wording_and_external_id_do_not_change_chemical_identity(
    operator_library,
) -> None:
    first = ExternalRetrosynthesisProposal(
        target_smiles="CCN",
        precursor_smiles="CC=O.N",
        external_proposal_id="first-external-id",
        claimed_transformation="claim one",
    )
    second = ExternalRetrosynthesisProposal(
        target_smiles="CCN",
        precursor_smiles="N.CC=O",
        external_proposal_id="second-external-id",
        claimed_transformation="contradictory claim two",
    )

    first_result = assess_external_retrosynthesis_proposal(first, operator_library)
    second_result = assess_external_retrosynthesis_proposal(second, operator_library)

    assert first_result.proposal_id == second_result.proposal_id
    assert first_result.assessment_id == second_result.assessment_id
    assert first_result.selected_operator_match_id == (
        second_result.selected_operator_match_id
    )


def test_valid_supplied_mapping_is_independently_checked(operator_library) -> None:
    result = assess_external_retrosynthesis_proposal(
        ExternalRetrosynthesisProposal(
            target_smiles="CCN",
            precursor_smiles="CC=O.N",
            mapped_reaction_smiles=SUPPLIED_FIRST_MAPPING,
        ),
        operator_library,
    )

    gates = {item.gate_id: item for item in result.gates}
    assert result.admission_eligible
    assert gates["reaction_side_consistency"].status == "pass"
    assert gates["atom_correspondence"].status == "pass"
    assert "EXTERNAL_MAPPING_REQUIRES_INDEPENDENT_VALIDATION" in result.warnings


def test_mapped_display_contradiction_is_invalid(operator_library) -> None:
    result = assess_external_retrosynthesis_proposal(
        ExternalRetrosynthesisProposal(
            target_smiles="CCO",
            precursor_smiles="CC=O.N",
            mapped_reaction_smiles=SUPPLIED_FIRST_MAPPING,
        ),
        operator_library,
    )

    gates = {item.gate_id: item for item in result.gates}
    assert result.status == "invalid"
    assert not result.admission_eligible
    assert result.reaction_signature is None
    assert gates["reaction_side_consistency"].status == "fail"
    assert gates["operator_support"].status == "not_run"


def test_unresolved_correspondence_is_preserved_without_operator_claim(
    operator_library,
) -> None:
    result = assess_external_retrosynthesis_proposal(
        ExternalRetrosynthesisProposal(
            target_smiles="CCN",
            precursor_smiles="CC.CN",
        ),
        operator_library,
    )

    assert result.status == "ambiguous"
    assert not result.admission_eligible
    assert not result.operator_matches
    assert any(
        item.gate_id == "atom_correspondence" and item.status == "unresolved"
        for item in result.gates
    )


def test_forward_challenge_is_separate_additional_evidence(operator_library) -> None:
    forward_library = build_forward_library_from_generic(operator_library)
    result = assess_external_retrosynthesis_proposal(
        ExternalRetrosynthesisProposal(
            target_smiles="CCN",
            precursor_smiles="CC=O.N",
        ),
        operator_library,
        forward_library=forward_library,
    )

    gates = {item.gate_id: item for item in result.gates}
    assert result.status == "forward_reproduced"
    assert result.strongest_evidence_tier == "forward_reproduced"
    assert gates["forward_reproduction"].status == "pass"
    assert gates["forward_competition"].status == "pass"


def test_complete_supported_route_builds_valid_review_tree(operator_library) -> None:
    result = assess_external_route_proposal(
        ExternalRouteProposal.from_dict(_route_value()), operator_library
    )

    assert result.status == "admitted_review_only"
    assert result.admission_eligible
    assert not result.actionable
    assert not result.unresolved_step_ids
    assert result.disconnected_step_ids == ()
    assert result.leaf_smiles == ("CC=O", "N", "O=C(O)c1ccccc1")
    assert result.admitted_route_tree is not None
    assert result.admitted_route_tree.reaction_count == 2
    assert result.admitted_route_tree.maximum_depth == 2
    assert validate_route_tree(result.admitted_route_tree).valid
    assert "LEAF_STOCK_STATUS_NOT_ASSESSED" in result.admitted_route_tree.warnings
    serialized = result.to_dict()
    assert serialized["admitted_route_tree"]["root"]["reaction"] is not None

    html = render_external_route_assessment_html(result)
    assert html.count("<svg") >= 4
    assert "External provenance (not validation)" in html
    assert "none_review_only" not in html


def test_planned_tree_retains_mapping_for_identity_stripped_parity(
    operator_library,
) -> None:
    candidate = disconnect_generic_target("CCN", operator_library, top_k=5)[0]
    step = RetrosynthesisRouteStep(
        step_id="planned-step",
        depth=1,
        product_smiles="CCN",
        precursor_smiles=("CC=O", "N"),
        step_cost=0.1,
        step_cost_components=(),
        candidate=candidate,
        product_node_id="root",
        precursor_node_ids=("leaf-1", "leaf-2"),
    )
    leaves = tuple(
        StartingMaterialAssessment(
            smiles=smiles,
            canonical_smiles=smiles,
            depth=1,
            molecular_weight=weight,
            terminal=True,
            terminal_reasons=("fixture",),
            terminal_evidence="fixture",
            catalog_role_status="reactant_supported",
            route_node_id=node_id,
        )
        for smiles, weight, node_id in (
            ("CC=O", 44.0, "leaf-1"),
            ("N", 17.0, "leaf-2"),
        )
    )
    tree = build_canonical_route_tree("CCN", "root", (step,), leaves)

    assert tree.root.reaction is not None
    evidence = tree.root.reaction.evidence
    assert evidence.reactants_smiles
    assert evidence.product_smiles_mapped
    assert "VALIDATED_MAPPED_REACTION_RETAINED" in evidence.warnings

    stripped = external_route_proposal_from_tree(tree)
    assert stripped.steps[0].proposal.mapped_reaction_smiles
    recovered = assess_external_route_proposal(stripped, operator_library)
    assert recovered.status == "admitted_review_only"
    assert recovered.admission_eligible


def test_disconnected_step_invalidates_route_without_erasing_step_evidence(
    operator_library,
) -> None:
    value = _route_value()
    value["steps"].append(
        {
            "external_step_id": "disconnected",
            "target_smiles": "CCO",
            "precursor_smiles": "CC=O",
        }
    )
    result = assess_external_route_proposal(
        ExternalRouteProposal.from_dict(value), operator_library
    )

    assert result.status == "invalid"
    assert result.disconnected_step_ids == ("disconnected",)
    assert result.admitted_route_tree is None
    assert len(result.step_assessments) == 3
    assert result.step_assessments[0].assessment.admission_eligible


def test_route_dependency_cycle_is_rejected(operator_library) -> None:
    route = ExternalRouteProposal.from_dict(
        {
            "target_smiles": "CC",
            "steps": [
                {
                    "external_step_id": "to-cc",
                    "target_smiles": "CC",
                    "precursor_smiles": "CO",
                },
                {
                    "external_step_id": "to-co",
                    "target_smiles": "CO",
                    "precursor_smiles": "CC",
                },
            ],
        }
    )
    result = assess_external_route_proposal(route, operator_library)

    connectivity = next(
        item for item in result.topology_gates if item.gate_id == "route_connectivity"
    )
    assert result.status == "invalid"
    assert connectivity.status == "fail"
    assert "cycle" in connectivity.summary


def test_external_input_cannot_supply_trusted_internal_ids() -> None:
    value = _route_value()["steps"][0]
    value["operator_id"] = "fabricated-operator"

    with pytest.raises(ValueError, match="cannot supply trusted fields"):
        ExternalRetrosynthesisProposal.from_dict(value)
