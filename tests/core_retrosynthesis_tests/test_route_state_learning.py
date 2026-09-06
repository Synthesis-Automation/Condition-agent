"""Literature-derived latent-state and route-order guidance regressions."""

from __future__ import annotations

from dataclasses import replace

from core_retrosynthesis import (
    GenericTemplateLibrary,
    LiteratureRouteActionSelector,
    ROUTE_STATE_LEARNING_VERSION,
    RouteStateActionObservation,
    RouteStateLearningCatalog,
    load_route_state_learning_catalog,
    load_route_state_learning_policy,
    mine_route_state_learning_catalog,
    save_route_state_learning_catalog,
)
from core_retrosynthesis.route_conversion import build_observed_route_tree
from core_retrosynthesis.route_core import build_route_core_projection

from .test_latent_state_search import CONVERGENT, SINGLE_CARRIER


def _record(split: str, patent_id: str) -> dict:
    return {
        "schema_version": "1.0",
        "route_id": f"route-{split}",
        "patent_id": patent_id,
        "split": split,
        "target_smiles": "CNC",
        "original_reaction_count": 1,
        "higher_level_reaction_count": 1,
        "higher_level_depth": 1,
        "abstraction_reduction": 0,
        "steps": [
            {
                "source_reaction_id": "release",
                "product_smiles": "CNC",
                "precursor_smiles": ["CC(=O)N(C)C"],
                "reactants_smiles": "[CH3:1][N:2]([CH3:3])[C:4](=[O:5])[CH3:6]",
                "reagents_smiles": "",
                "product_smiles_mapped": "[CH3:1][NH:2][CH3:3]",
                "reaction_smiles": (
                    "[CH3:1][N:2]([CH3:3])[C:4](=[O:5])[CH3:6]"
                    ">>[CH3:1][NH:2][CH3:3]"
                ),
                "abstracted_reaction_smiles": "",
            },
        ],
    }


def _empty_library() -> GenericTemplateLibrary:
    return GenericTemplateLibrary(
        templates=(),
        source_row_count=0,
        accepted_observation_count=0,
        rejection_counts={},
        definition={},
    )


def _catalog(operator_id: str) -> RouteStateLearningCatalog:
    return RouteStateLearningCatalog(
        definition_id=ROUTE_STATE_LEARNING_VERSION,
        schema_version="1.0",
        source_path="fixture",
        source_sha256="abc",
        operator_library_path="fixture",
        route_counts={"train": 1},
        action_counts={"train:protection_state_change": 1},
        dependency_counts={},
        train_state_operator_support={operator_id: 1},
        train_state_operator_patent_support={operator_id: 1},
        train_reverse_operator_transition_support={},
        train_reverse_operator_transition_patent_support={},
        train_action_transition_support={},
        heldout_metrics={},
        action_samples=(),
        ordering_samples=(),
    )


def test_policy_is_versioned_and_non_authoritative() -> None:
    policy = load_route_state_learning_policy()

    assert policy.definition_id == ROUTE_STATE_LEARNING_VERSION
    assert policy.ranking_influence == (
        "opt_in_validated_state_ordering_and_portfolio_reservation_only"
    )
    assert "amine_deprotection_like" in policy.state_release_pattern_ids


def test_mining_uses_train_only_for_state_operator_support(
    monkeypatch, tmp_path
) -> None:
    train = build_route_core_projection(
        build_observed_route_tree(_record("train", "PAT-TRAIN"))
    )
    test = build_route_core_projection(
        build_observed_route_tree(_record("test", "PAT-TEST"))
    )
    source = tmp_path / "source.jsonl"
    source.write_text("fixture", encoding="utf-8")
    monkeypatch.setattr(
        "core_retrosynthesis.route_state_learning.iter_route_core_projections",
        lambda _: iter((train, test)),
    )
    monkeypatch.setattr(
        "core_retrosynthesis.route_state_learning._operator_indexes",
        lambda _: (
            {
                "PAT-TRAIN:release": ("OP:STATE",),
                "PAT-TEST:release": ("OP:HELDOUT",),
            },
            {},
        ),
    )

    catalog = mine_route_state_learning_catalog(source, _empty_library())

    assert catalog.route_counts == {"test": 1, "train": 1}
    assert catalog.train_state_operator_support == {"OP:STATE": 1}
    assert "OP:HELDOUT" not in catalog.train_state_operator_support
    assert catalog.heldout_metrics["test"]["mapped_state_action_count"] == 1
    assert catalog.heldout_metrics["test"]["mapped_state_action_coverage"] == 0.0
    assert any(
        item.action_class == "protection_state_change"
        for item in catalog.action_samples
    )


def test_catalog_round_trip_preserves_typed_samples(tmp_path) -> None:
    observation = RouteStateActionObservation(
        observation_id="obs",
        route_core_id="route",
        patent_id="patent",
        split="train",
        source_reaction_id="reaction",
        retrosynthetic_depth=1,
        action_class="protection_state_change",
        pattern_ids=("amine_deprotection_like",),
        operator_ids=("OP:STATE",),
        exact_core_key="exact",
        typed_core_key="typed",
        shape_core_key="shape",
        reaction_smiles="CN(C)C=O>>CNC",
    )
    catalog = replace(_catalog("OP:STATE"), action_samples=(observation,))
    destination = tmp_path / "catalog.json.gz"

    save_route_state_learning_catalog(catalog, destination)
    first_bytes = destination.read_bytes()
    save_route_state_learning_catalog(catalog, destination)

    assert load_route_state_learning_catalog(destination) == catalog
    assert destination.read_bytes() == first_bytes


def test_selector_reserves_only_supported_validated_state_action() -> None:
    supported = replace(
        SINGLE_CARRIER,
        template_id="supported-state",
        operator_id="OP:STATE",
        strategy_id="",
        strategic_class="complexity_increasing",
    )
    unsupported = replace(
        SINGLE_CARRIER,
        template_id="unsupported-state",
        operator_id="OP:OTHER",
        strategy_id="",
        strategic_class="complexity_increasing",
    )
    candidates = (
        CONVERGENT,
        replace(CONVERGENT, template_id="second"),
        replace(CONVERGENT, template_id="third"),
        unsupported,
        supported,
    )
    selector = LiteratureRouteActionSelector(_catalog("OP:STATE"))

    selected = selector("CNC", candidates, top_k=3)

    assert selected[0] is CONVERGENT
    assert supported in selected
    assert unsupported not in selected
    assert all(item in candidates for item in selected)
