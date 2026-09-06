"""Baseline-preserving latent-route portfolio regressions."""

from __future__ import annotations

from dataclasses import replace

from core_retrosynthesis import (
    BaselinePreservingLatentPortfolio,
    LiteratureRouteActionSelector,
    RouteStateLearningCatalog,
    build_baseline_preserving_latent_portfolio,
    build_dataset_latent_portfolio_review,
    load_latent_route_portfolio_policy,
    plan_baseline_preserving_latent_portfolio,
    plan_multistep_routes,
    write_latent_route_portfolio,
    write_dataset_latent_portfolio_review,
)
from core_retrosynthesis.cli import _parser

from .test_multistep import _LiteratureIndex, _candidate, _expander


def _catalog() -> RouteStateLearningCatalog:
    return RouteStateLearningCatalog(
        definition_id="route_state_learning.v1",
        schema_version="1.0",
        source_path="fixture",
        source_sha256="abc",
        operator_library_path="fixture",
        route_counts={"train": 1},
        action_counts={"train:protection_state_change": 1},
        dependency_counts={},
        train_state_operator_support={"OP:STATE": 3},
        train_state_operator_patent_support={"OP:STATE": 2},
        train_reverse_operator_transition_support={},
        train_reverse_operator_transition_patent_support={},
        train_action_transition_support={},
        heldout_metrics={},
        action_samples=(),
        ordering_samples=(),
    )


def _results():
    target = "CCNCC"
    baseline_candidate = _candidate(target, "CCN.CC")
    state_candidate = replace(
        _candidate(target, "CCN(CC)C=O"),
        operator_id="OP:STATE",
        strategy_id="",
        strategic_class="complexity_increasing",
    )
    literature = _LiteratureIndex(known=("CCN", "CC", "CCN(CC)C=O"))
    baseline = plan_multistep_routes(
        target,
        object(),
        literature,
        max_depth=1,
        molecular_weight_threshold=20.0,
        expander=_expander({target: (baseline_candidate,)}),
    )
    exploratory = plan_multistep_routes(
        target,
        object(),
        literature,
        max_depth=1,
        molecular_weight_threshold=20.0,
        expander=_expander({target: (state_candidate,)}),
    )
    return baseline, exploratory


def test_policy_requires_a_separate_baseline_preserving_lane() -> None:
    policy = load_latent_route_portfolio_policy()

    assert policy.preserve_baseline_routes is True
    assert policy.maximum_exploratory_routes == 1
    assert policy.ranking_influence == "separate_opt_in_exploratory_lane"


def test_latent_alternative_is_appended_without_reordering_baseline() -> None:
    baseline, exploratory = _results()

    portfolio = build_baseline_preserving_latent_portfolio(
        baseline, exploratory, _catalog()
    )

    assert portfolio.baseline is baseline
    assert portfolio.baseline_preserved is True
    assert portfolio.augmented_route_ids[: len(portfolio.baseline_route_ids)] == (
        portfolio.baseline_route_ids
    )
    assert len(portfolio.alternatives) == 1
    assert portfolio.alternatives[0].supported_state_operator_ids == (
        "OP:STATE",
    )


def test_duplicate_or_unsupported_routes_are_not_added() -> None:
    baseline, exploratory = _results()

    duplicate = build_baseline_preserving_latent_portfolio(
        baseline, baseline, _catalog()
    )
    unsupported = build_baseline_preserving_latent_portfolio(
        baseline,
        exploratory,
        replace(
            _catalog(),
            train_state_operator_support={},
            train_state_operator_patent_support={},
        ),
    )

    assert duplicate.alternatives == ()
    assert unsupported.alternatives == ()


def test_two_search_wrapper_keeps_hooks_out_of_baseline() -> None:
    baseline, exploratory = _results()
    calls = []

    def planner(target, library, literature, **kwargs):
        del target, library, literature
        calls.append(kwargs)
        return exploratory if kwargs.get("route_action_selector") else baseline

    result = plan_baseline_preserving_latent_portfolio(
        baseline.target_smiles,
        object(),
        object(),
        _catalog(),
        planner=planner,
        max_depth=1,
    )

    assert isinstance(result, BaselinePreservingLatentPortfolio)
    assert "route_action_selector" not in calls[0]
    assert isinstance(calls[1]["route_action_selector"], LiteratureRouteActionSelector)
    assert calls[1]["search_guidance"] is None
    assert result.baseline is baseline


def test_portfolio_report_is_structure_rich(tmp_path) -> None:
    baseline, exploratory = _results()
    portfolio = build_baseline_preserving_latent_portfolio(
        baseline, exploratory, _catalog()
    )

    summary = write_latent_route_portfolio(
        portfolio,
        tmp_path / "portfolio.json",
        tmp_path / "portfolio.html",
    )

    document = (tmp_path / "portfolio.html").read_text(encoding="utf-8")
    assert summary["baseline_preserved"] is True
    assert summary["novel_latent_route_count"] == 1
    assert "Novel latent-state lane" in document
    assert "Baseline preserved" in document
    assert "<svg" in document


def test_cli_exposes_baseline_preserving_latent_portfolio() -> None:
    arguments = _parser().parse_args(
        [
            "plan-latent-portfolio",
            "library.json.gz",
            "stock.sqlite",
            "catalog.json.gz",
            "CCNCC",
            "portfolio.json",
            "portfolio.html",
        ]
    )

    assert arguments.command == "plan-latent-portfolio"
    assert arguments.route_state_catalog == "catalog.json.gz"


def test_dataset_review_keeps_baseline_and_adds_only_novel_state_route(
    tmp_path,
) -> None:
    baseline, exploratory = _results()
    target = {
        "target_id": "case-1",
        "name": "Case 1",
        "category": "fixture",
        "smiles": baseline.target_smiles,
    }
    base_artifact = {
        "evaluation_id": "base",
        "summary": {"planner_solved_count": 1},
        "cases": [
            {
                "case_id": "case-1",
                "selection_rank": 1,
                "panel_case": {
                    "target": target,
                    "baseline": baseline.to_dict(),
                    "reference_reactions": [],
                },
            }
        ],
    }
    exploratory_artifact = {
        "evaluation_id": "state",
        "cases": [
            {
                "case_id": "case-1",
                "selection_rank": 1,
                "panel_case": {
                    "target": target,
                    "baseline": exploratory.to_dict(),
                },
            }
        ],
    }

    review = build_dataset_latent_portfolio_review(
        base_artifact, exploratory_artifact, _catalog()
    )
    summary = write_dataset_latent_portfolio_review(
        review,
        tmp_path / "dataset.json",
        tmp_path / "dataset.html",
    )

    assert review["summary"]["baseline_preserved_count"] == 1
    assert review["summary"]["targets_with_novel_latent_route"] == 1
    assert review["cases"][0]["baseline_route_ids"] == [
        baseline.routes[0].route_id
    ]
    assert summary["summary"]["targets_with_novel_latent_route"] == 1
    document = (tmp_path / "dataset.html").read_text(encoding="utf-8")
    assert "Planner baseline — preserved" in document
    assert "Additional latent-state route" in document
