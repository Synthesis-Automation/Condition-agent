"""Data-driven route-action policy and multistep integration regressions."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from core_retrosynthesis.route_action_evaluation import (
    RouteActionEvaluationConfig,
    evaluate_route_actions,
)
from core_retrosynthesis.route_action_policy import (
    RouteActionPolicyModel,
    RoutePolicyCandidate,
    RoutePolicyExample,
    build_route_policy_examples,
    load_route_action_policy,
    load_route_action_policy_definition,
    save_route_action_policy,
)
from core_retrosynthesis.route_conversion import build_observed_route_tree

from .test_multistep import _LiteratureIndex, _candidate, _expander
from .test_route_action_evaluation import (
    PRECURSOR,
    PRODUCT,
    _candidate as _replay_candidate,
    _empty_library,
    _record,
)
from core_retrosynthesis.multistep import plan_multistep_routes


def _policy_candidate(
    candidate_id: str,
    *,
    rank: int,
    strategic_class: str,
) -> RoutePolicyCandidate:
    return RoutePolicyCandidate(
        candidate_id=candidate_id,
        candidate_rank=rank,
        precursor_smiles="C.C" if strategic_class == "strategic" else "N.O",
        abstraction_level="L2",
        template_id=f"template-{candidate_id}",
        operator_id="operator",
        synthon_signature="synthon",
        score=0.9,
        structural_score_band=4,
        strategic_complexity_score=(0.8 if strategic_class == "strategic" else 0.0),
        strategic_class=strategic_class,
        independent_reference_support=2,
        precursor_compatibility_disposition="pass",
        supervision_label=(
            "observed_exact" if strategic_class == "strategic" else "unchosen_alternative"
        ),
        source_patent_precedent_overlap=False,
    )


def _training_examples() -> tuple[RoutePolicyExample, ...]:
    values = []
    for index in range(4):
        tactical = _policy_candidate(
            f"tactical-{index}", rank=1, strategic_class="tactical"
        )
        strategic = _policy_candidate(
            f"strategic-{index}", rank=2, strategic_class="strategic"
        )
        values.append(
            RoutePolicyExample(
                example_id=f"example-{index}",
                route_id=f"route-{index}",
                patent_id=f"US{index}A1",
                split="train",
                target_product_smiles="CCCCCCCC",
                retrosynthetic_depth=0,
                observed_remaining_steps=2,
                maximum_route_depth=3,
                selected_candidate_ids=(strategic.candidate_id,),
                label_source="observed_exact",
                label_strength=1.0,
                candidates=(tactical, strategic),
            )
        )
    return tuple(values)


def _trained_policy() -> RouteActionPolicyModel:
    model = RouteActionPolicyModel()
    report = model.fit(_training_examples())
    assert report.final_loss < report.initial_loss
    model.residual_scale = 1.0
    return model


def test_policy_definition_and_serialization_are_deterministic(tmp_path: Path) -> None:
    definition = load_route_action_policy_definition()
    assert definition.definition_id == "route_action_policy.v1@1.2"
    model = _trained_policy()
    first = tmp_path / "first.json.gz"
    second = tmp_path / "second.json.gz"

    save_route_action_policy(model, first)
    save_route_action_policy(model, second)

    assert first.read_bytes() == second.read_bytes()
    restored = load_route_action_policy(first)
    assert restored.model_id == model.model_id
    assert restored.to_dict() == model.to_dict()


def test_replay_choices_become_listwise_route_policy_examples() -> None:
    def searcher(*args: object, **kwargs: object):
        return (
            (
                _replay_candidate(
                    "CC",
                    operator="other",
                    site="other",
                    synthon="other",
                    precedent="US9A1:1",
                ),
                _replay_candidate(
                    PRECURSOR,
                    operator="expected",
                    site="expected",
                    synthon="expected",
                    precedent="US8A1:1",
                ),
            ),
            {},
        )

    tree = build_observed_route_tree(_record())
    observed = evaluate_route_actions(
        tree,
        _empty_library(),
        config=RouteActionEvaluationConfig(use_hierarchical_ranking=False),
        searcher=searcher,
    )
    exact = replace(
        observed.steps[0].candidates[1],
        exact_precursor_match=True,
        supervision_label="observed_exact",
    )
    evaluation = replace(
        observed,
        steps=(
            replace(
                observed.steps[0],
                candidates=(observed.steps[0].candidates[0], exact),
                candidate_count=2,
            ),
        ),
    )

    examples = build_route_policy_examples((evaluation,))

    assert len(examples) == 1
    assert examples[0].label_source == "observed_exact"
    assert len(examples[0].selected_candidate_ids) == 1
    assert examples[0].target_product_smiles == PRODUCT


def test_policy_without_sufficient_validation_has_zero_planner_influence() -> None:
    model = RouteActionPolicyModel()
    report = model.fit(_training_examples())

    assert report.policy_active is False
    assert report.selected_residual_scale == 0.0
    assert report.activation_reason == "insufficient_validation_examples:0<20"
    assert model.planner_probability_deficit_weight == 0.0


def test_learned_policy_reorders_validated_multistep_actions() -> None:
    model = _trained_policy()
    tactical = replace(
        _candidate("CCCCCCCC", "N.O", score=0.9),
        strategic_class="tactical",
        strategic_complexity_score=0.0,
    )
    strategic = replace(
        _candidate("CCCCCCCC", "C.C", score=0.9),
        strategic_class="strategic",
        strategic_complexity_score=0.8,
    )

    result = plan_multistep_routes(
        "CCCCCCCC",
        object(),
        _LiteratureIndex(),
        max_depth=2,
        molecular_weight_threshold=100.0,
        per_step_top_k=2,
        route_action_policy=model,
        expander=_expander({"CCCCCCCC": (tactical, strategic)}),
    )

    assert result.route_action_policy_model_id == model.model_id
    assert result.route_action_policy_active is True
    assert result.route_action_policy_residual_scale == 1.0
    assert result.diagnostics.route_policy_scored_actions == 2
    assert result.diagnostics.route_policy_reordered_expansions == 1
    assert result.routes[0].steps[0].candidate.precursor_smiles == "C.C"
    assert result.routes[0].steps[0].route_policy_rank == 1
    assert result.routes[0].steps[0].original_candidate_rank == 2
    assert "route_policy_probability_deficit" in dict(
        result.routes[0].steps[0].step_cost_components
    )
