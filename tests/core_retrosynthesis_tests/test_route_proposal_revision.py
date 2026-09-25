"""Graph-preserving branch edits and declared starting-material constraints."""

from dataclasses import replace

import pytest

from core_retrosynthesis.external_route_admission import ExternalRouteProposal, ExternalRouteStepProposal, assess_external_route_proposal
from core_retrosynthesis.route_proposal_revision import (
    assess_declared_route_materials, external_route_step_neighbors, revise_external_route_proposal,
)
from core_retrosynthesis.generic_library import build_generic_library
from .test_external_proposal_admission import _route_value, _row, FIRST_REACTION, SECOND_REACTION


@pytest.fixture(scope="module")
def operator_library():
    return build_generic_library((_row(FIRST_REACTION, 1), _row(SECOND_REACTION, 2)), levels=("L0", "L1", "L2"))


def test_revision_replaces_only_selected_branch_and_retains_target() -> None:
    source = ExternalRouteProposal.from_dict(_route_value())
    before = source.to_dict()
    new = ExternalRouteStepProposal.from_dict({"external_step_id": "step-1", "target_smiles": "CCN", "precursor_smiles": "CCBr.N"})
    result = revise_external_route_proposal(source, remove_step_ids=("step-1",), replacement_steps=(new,))
    assert result.proposal.target_smiles == source.target_smiles
    assert result.proposal.steps[0] is source.steps[0]
    assert result.replaced_step_ids == ("step-1",)
    assert result.preserved_step_ids == ("step-2",)
    assert source.to_dict() == before
    assert external_route_step_neighbors(result.proposal, "step-1")["downstream_step_ids"] == ["step-2"]
    assert external_route_step_neighbors(result.proposal, "step-2")["upstream_step_ids"] == ["step-1"]


@pytest.mark.parametrize("removals,replacements,message", [
    (("missing",), (), "unknown"), (("step-1", "step-1"), (), "unique"),
    ((), ("step-1",), "explicit removal"),
    (("step-1",), ("step-1", "step-1"), "unique"),
    (("step-1",), ("step-1",), "change at least"),
])
def test_revision_rejects_ambiguous_edits_and_noops(removals, replacements, message) -> None:
    source = ExternalRouteProposal.from_dict(_route_value())
    by_id = {step.external_step_id: step for step in source.steps}
    with pytest.raises(ValueError, match=message):
        revise_external_route_proposal(source, remove_step_ids=removals,
                                       replacement_steps=tuple(by_id[key] for key in replacements))


def test_unavailable_intermediate_can_be_made_but_new_stock_is_not_assumed(operator_library) -> None:
    route = ExternalRouteProposal.from_dict(_route_value())
    terminal_amine = replace(route, steps=route.steps[:1])
    before = assess_external_route_proposal(terminal_amine, operator_library)
    declared = assess_declared_route_materials(before, ("NCC",))
    assert declared.status == "violated"
    assert declared.blocked_leaf_smiles == ("CCN",)
    revision = revise_external_route_proposal(terminal_amine, remove_step_ids=(), replacement_steps=route.steps[1:])
    after = assess_external_route_proposal(revision.proposal, operator_library)
    resolved = assess_declared_route_materials(after, ("NCC",))
    assert resolved.status == "satisfied_for_declared_constraints"
    assert resolved.stock_availability == "not_assessed"
    assert revision.added_step_ids == ("step-1",)
    assert not after.actionable


def test_broken_topology_cannot_clear_material_constraints(operator_library) -> None:
    route = ExternalRouteProposal.from_dict(_route_value())
    # Removing the target step leaves a disconnected upstream branch.
    edited = revise_external_route_proposal(route, remove_step_ids=("step-2",), replacement_steps=()).proposal
    result = assess_external_route_proposal(edited, operator_library)
    assert result.status == "invalid"
    assert assess_declared_route_materials(result, ("CCN",)).status == "unknown"


def test_declared_material_identity_uses_graphs_and_rejects_invalid_input(operator_library) -> None:
    result = assess_external_route_proposal(ExternalRouteProposal.from_dict(_route_value()), operator_library)
    assert assess_declared_route_materials(result, ("C(C)=O", "CC=O")).blocked_leaf_smiles == ("CC=O",)
    assert assess_declared_route_materials(result, ()).status == "not_requested"
    for value in ("invalid", "CCO.N", ""):
        with pytest.raises(ValueError, match="valid molecular graph"):
            assess_declared_route_materials(result, (value,))
