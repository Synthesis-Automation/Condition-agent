"""Executable two-step strategy replay regressions."""

from __future__ import annotations

import pytest

from core_retrosynthesis.cli import _parser
from core_retrosynthesis.coupled_strategy_search import (
    compare_coupled_strategy_policies,
    search_coupled_strategy_replays,
    write_coupled_strategy_policy_comparison,
)


def _occurrence(
    occurrence_id: str,
    *,
    first_reaction: str,
    second_reaction: str,
    intermediate: str,
    target: str,
    relationship: str,
    dependency: str,
    admission: str,
    score: float = 1.0,
) -> dict[str, object]:
    return {
        "occurrence_id": occurrence_id,
        "exact_strategy_id": f"exact:{occurrence_id}",
        "typed_strategy_id": f"typed:{occurrence_id}",
        "relationship_class": relationship,
        "dependency_class": dependency,
        "admission_class": admission,
        "coupling_score": score,
        "first_reaction_smiles": first_reaction,
        "second_reaction_smiles": second_reaction,
        "intermediate_smiles": intermediate,
        "final_product_smiles": target,
        "overall_reaction_smiles": "",
        "first_source_reaction_id": f"{occurrence_id}:1",
        "second_source_reaction_id": f"{occurrence_id}:2",
        "patent_id": f"PATENT:{occurrence_id}",
        "source_route_id": f"ROUTE:{occurrence_id}",
    }


def _nitration_reduction() -> dict[str, object]:
    return _occurrence(
        "nitration-reduction",
        first_reaction=(
            "c1ccccc1.O=[N+]([O-])O>>O=[N+]([O-])c1ccccc1"
        ),
        second_reaction="O=[N+]([O-])c1ccccc1>>Nc1ccccc1",
        intermediate="O=[N+]([O-])c1ccccc1",
        target="Nc1ccccc1",
        relationship="handle_progression",
        dependency="created_handle_consumed",
        admission="strict",
    )


def _activation_substitution() -> dict[str, object]:
    return _occurrence(
        "activation-substitution",
        first_reaction="CCO.CS(=O)(=O)Cl>>CCOS(C)(=O)=O",
        second_reaction="CCOS(C)(=O)=O.N>>CCN",
        intermediate="CCOS(C)(=O)=O",
        target="CCN",
        relationship="handle_progression",
        dependency="activation_then_conversion",
        admission="strict",
        score=0.95,
    )


def _temporary_group() -> dict[str, object]:
    return _occurrence(
        "temporary-group",
        first_reaction="Oc1ccccc1.CC(=O)Cl>>CC(=O)Oc1ccccc1",
        second_reaction="CC(=O)Oc1ccccc1>>Oc1ccccc1",
        intermediate="CC(=O)Oc1ccccc1",
        target="Oc1ccccc1",
        relationship="handle_progression",
        dependency="temporary_group_removed",
        admission="review",
        score=0.65,
    )


def _independent_sites() -> dict[str, object]:
    return _occurrence(
        "independent-sites",
        first_reaction=(
            "CC(=O)Cl.NCC(=O)OC>>CC(=O)NCC(=O)OC"
        ),
        second_reaction="CC(=O)NCC(=O)OC>>CC(=O)NCC(=O)O",
        intermediate="CC(=O)NCC(=O)OC",
        target="CC(=O)NCC(=O)O",
        relationship="independent_sites",
        dependency="independent_sites",
        admission="rejected",
        score=0.0,
    )


def test_nitration_reduction_replays_as_two_physical_steps() -> None:
    occurrence = _nitration_reduction()

    comparison = compare_coupled_strategy_policies(
        "c1ccc(cc1)N", [occurrence]
    )

    assert len(comparison.v1.actions) == 1
    assert len(comparison.v2.actions) == 1
    action = comparison.v2.actions[0]
    assert action.dependency_class == "created_handle_consumed"
    assert action.intermediate_smiles == "O=[N+]([O-])c1ccccc1"
    assert len(action.physical_steps) == 2
    assert [
        step.forward_step_number for step in action.retrosynthetic_steps
    ] == [2, 1]
    assert action.validation_status == "validated_exact_replay"
    assert action.logical_reaction_smiles.endswith(">>Nc1ccccc1")


def test_activation_then_conversion_is_admitted_by_v2() -> None:
    result = search_coupled_strategy_replays(
        "CCN",
        [_activation_substitution()],
        policy_version="v2",
    )

    assert len(result.actions) == 1
    assert result.actions[0].dependency_class == "activation_then_conversion"
    assert result.diagnostics.valid_action_count == 1
    assert result.actions[0].physical_steps[0].product_smiles == "CCOS(C)(=O)=O"


def test_v2_suppresses_temporary_group_macro_admitted_by_v1() -> None:
    comparison = compare_coupled_strategy_policies(
        "Oc1ccccc1", [_temporary_group()]
    )

    assert len(comparison.v1.actions) == 1
    assert len(comparison.v2.actions) == 0
    assert comparison.v1.diagnostics.policy_admitted_count == 1
    assert comparison.v2.diagnostics.policy_rejected_count == 1
    assert (
        comparison.v1.actions[0].dependency_class
        == "temporary_group_removed"
    )


def test_independent_site_pair_is_never_a_macro_action() -> None:
    comparison = compare_coupled_strategy_policies(
        "CC(=O)NCC(=O)O", [_independent_sites()]
    )

    assert not comparison.v1.actions
    assert not comparison.v2.actions
    assert comparison.v1.diagnostics.policy_rejected_count == 1
    assert comparison.v2.diagnostics.policy_rejected_count == 1


def test_broken_intermediate_chain_is_rejected_after_policy_admission() -> None:
    broken = dict(_nitration_reduction())
    broken["intermediate_smiles"] = "CCO"

    result = search_coupled_strategy_replays(
        "Nc1ccccc1", [broken], policy_version="v2"
    )

    assert not result.actions
    assert result.diagnostics.policy_admitted_count == 1
    assert result.diagnostics.invalid_chain_count == 1


def test_report_sample_and_serialization_contract() -> None:
    report = {"sample": [_nitration_reduction(), _activation_substitution()]}

    result = search_coupled_strategy_replays(
        "Nc1ccccc1", report, policy_version="v2", top_k=1
    )
    serialized = result.to_dict()

    assert serialized["schema_version"] == "1.0"
    assert serialized["diagnostics"]["source_occurrence_count"] == 2
    assert serialized["actions"][0]["retrosynthetic_step_numbers"] == [2, 1]


def test_comparison_writer_and_cli_contract(tmp_path) -> None:
    comparison = compare_coupled_strategy_policies(
        "Nc1ccccc1", [_nitration_reduction()]
    )
    destination = tmp_path / "comparison.json"

    summary = write_coupled_strategy_policy_comparison(
        comparison, destination
    )
    arguments = _parser().parse_args(
        [
            "replay-coupled-route-strategies",
            "report.json",
            "Nc1ccccc1",
            "comparison.json",
            "--top-k",
            "4",
        ]
    )

    assert summary["v1_action_count"] == 1
    assert summary["v2_action_count"] == 1
    assert destination.read_text(encoding="utf-8").endswith("\n")
    assert arguments.command == "replay-coupled-route-strategies"
    assert arguments.top_k == 4


def test_invalid_arguments_are_rejected() -> None:
    with pytest.raises(ValueError, match="policy version"):
        search_coupled_strategy_replays(
            "CC", [], policy_version="v3"  # type: ignore[arg-type]
        )
    with pytest.raises(ValueError, match="target SMILES"):
        search_coupled_strategy_replays(
            "not-a-smiles", [], policy_version="v2"
        )
    with pytest.raises(ValueError, match="sample sequence"):
        search_coupled_strategy_replays(
            "CC", {"sample": None}, policy_version="v2"
        )
