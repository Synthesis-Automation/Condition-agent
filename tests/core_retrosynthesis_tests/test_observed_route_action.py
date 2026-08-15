"""Evidence-faceted observed route-action label regressions."""

from __future__ import annotations

from core_retrosynthesis.observed_route_action import (
    ObservedRouteActionLabel,
    build_observed_route_action_label,
    load_observed_route_action_label_policy,
)


SULFONAMIDE_FORMATION = (
    "[NH2:1][c:4]1[cH:5][cH:6][cH:7][cH:8][cH:9]1."
    "[S:2]([CH3:10])(=[O:11])(=[O:12])[Cl:3]>>"
    "[NH:1]([c:4]1[cH:5][cH:6][cH:7][cH:8][cH:9]1)"
    "[S:2]([CH3:10])(=[O:11])=[O:12]"
)


def test_policy_is_versioned_and_rejects_only_declared_review_reasons() -> None:
    policy = load_observed_route_action_label_policy()

    assert policy.definition_id == "observed_route_action_label.v1@1.1"
    assert policy.allowed_review_reasons == ("not_all_edits_graph_checked",)
    assert "active_atom_mapping_complete" in policy.required_passed_checks


def test_departing_leaving_group_can_label_strategy_without_admitting_operator() -> None:
    label = build_observed_route_action_label(
        SULFONAMIDE_FORMATION,
        route_product_smiles="CS(=O)(=O)Nc1ccccc1",
        reaction_id="observed-sulfonamide",
        reference_id="US00000001A1",
    )

    assert label.core_quality_status == "review"
    assert label.core_review_reasons == ("not_all_edits_graph_checked",)
    assert label.departing_edit_descriptors == ("broken:S-Cl:SINGLE:NONE",)
    assert label.product_site_verified
    assert label.retained_edits_verified
    assert label.synthon_partition_verified
    assert label.exact_precursors_verified
    assert label.strategy_verified
    assert label.realization_verified
    assert not label.operator_roundtrip_verified
    assert label.search_eligible
    assert label.strategy_id and label.strategy_id.startswith("STRAT1:")
    assert ObservedRouteActionLabel.from_dict(label.to_dict()) == label


def test_conflicting_product_identity_disables_all_learning_facets() -> None:
    label = build_observed_route_action_label(
        SULFONAMIDE_FORMATION,
        route_product_smiles="CC",
    )

    assert not label.product_site_verified
    assert not label.retained_edits_verified
    assert not label.exact_precursors_verified
    assert not label.search_eligible
    assert "target_identity_mismatch" in label.limitations
