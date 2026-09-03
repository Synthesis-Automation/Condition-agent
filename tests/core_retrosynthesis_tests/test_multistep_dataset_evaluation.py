"""Bounded dataset multistep evaluation regressions."""

from __future__ import annotations

from pathlib import Path

from core_retrosynthesis import GenericTemplateLibrary
from core_retrosynthesis.cli import _parser
from core_retrosynthesis.multistep import plan_multistep_routes
from core_retrosynthesis.multistep_dataset_evaluation import (
    MultistepDatasetEvaluationConfig,
    evaluate_partition_review_routes,
    write_multistep_dataset_evaluation,
)

from .test_multistep import _LiteratureIndex, _candidate, _expander
from .test_operator_coverage_comparison import _tree


def test_evaluates_exact_observed_route_recovery(tmp_path: Path) -> None:
    tree = _tree()
    review = {
        "cases": [
            {
                "case_id": "case",
                "selection_rank": 1,
                "route_tree": tree.to_dict(),
            }
        ]
    }
    library = GenericTemplateLibrary(
        templates=(),
        source_row_count=1,
        accepted_observation_count=1,
        rejection_counts={},
        definition={"core_admission_policy": "validated_departures"},
    )

    def planner(target, loaded_library, literature_index, **kwargs):
        del loaded_library, kwargs
        return plan_multistep_routes(
            target,
            object(),
            literature_index,
            max_depth=1,
            molecular_weight_threshold=150.0,
            expander=_expander(
                {target: (_candidate(target, "CC(=O)O.CCN"),)}
            ),
        )

    evaluation = evaluate_partition_review_routes(
        review,
        library,
        _LiteratureIndex(),
        config=MultistepDatasetEvaluationConfig(max_depth=1),
        planner=planner,
    )

    assert evaluation.summary == {
        "target_count": 1,
        "planner_solved_count": 1,
        "planner_partial_only_count": 0,
        "exact_root_action_recovery_count": 1,
        "any_observed_action_recovery_count": 1,
        "full_observed_route_recovery_count": 1,
        "verified_observed_action_count": 1,
        "matched_observed_action_count": 1,
    }
    case = evaluation.cases[0]
    assert case.exact_root_action_recovered is True
    assert case.full_observed_route_recovered is True

    report = write_multistep_dataset_evaluation(
        evaluation,
        tmp_path / "evaluation.json",
        tmp_path / "evaluation.html",
    )
    assert report["summary"]["planner_solved_count"] == 1
    document = (tmp_path / "evaluation.html").read_text(encoding="utf-8")
    assert "Exact root actions" in document
    assert "BOUNDED_SEARCH_IS_NOT_A_COMPLETE_RETROSYNTHESIS_PROOF" in document
    assert "Observed precedent" in document
    assert "<svg" in document


def test_cli_exposes_partition_review_route_evaluation() -> None:
    arguments = _parser().parse_args(
        [
            "evaluate-partition-review-routes",
            "review.json",
            "library.json",
            "stock.sqlite",
            "evaluation.json",
            "evaluation.html",
            "--max-depth",
            "4",
            "--allow-untyped-literature-terminals",
        ]
    )

    assert arguments.command == "evaluate-partition-review-routes"
    assert arguments.max_depth == 4
    assert arguments.allow_untyped_literature_terminals is True
