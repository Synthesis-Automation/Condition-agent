"""Chemist-facing promoted route-action review regressions."""

from __future__ import annotations

from dataclasses import replace

from core_retrosynthesis.generic_models import GenericTemplateLibrary
from core_retrosynthesis.route_action_evaluation import (
    RouteActionEvaluationConfig,
    evaluate_route_actions,
)
from core_retrosynthesis.route_action_review import (
    is_promoted_route_action,
    render_route_action_review_html,
    sample_promoted_route_actions,
)
from core_retrosynthesis.route_conversion import build_observed_route_tree

from .test_observed_route_action import SULFONAMIDE_FORMATION


def _evaluation():
    product = "CS(=O)(=O)Nc1ccccc1"
    record = {
        "schema_version": "1.0",
        "route_id": "US00000002A1_0",
        "patent_id": "US00000002A1",
        "split": "test",
        "target_smiles": product,
        "original_reaction_count": 1,
        "higher_level_reaction_count": 1,
        "higher_level_depth": 1,
        "abstraction_reduction": 0,
        "steps": [
            {
                "source_reaction_id": "sulfonamide-1",
                "product_smiles": product,
                "precursor_smiles": SULFONAMIDE_FORMATION.split(">>")[0].split("."),
                "reactants_smiles": SULFONAMIDE_FORMATION.split(">>")[0],
                "reagents_smiles": "",
                "product_smiles_mapped": SULFONAMIDE_FORMATION.split(">>")[1],
                "reaction_smiles": SULFONAMIDE_FORMATION,
                "abstracted_reaction_smiles": "",
            }
        ],
    }
    library = GenericTemplateLibrary(
        templates=(),
        source_row_count=0,
        accepted_observation_count=0,
        rejection_counts={},
        definition={},
    )
    return evaluate_route_actions(
        build_observed_route_tree(record),
        library,
        config=RouteActionEvaluationConfig(run_search=False),
    )


def test_promoted_review_sampling_and_html_are_deterministic() -> None:
    evaluation = _evaluation()

    assert is_promoted_route_action(evaluation.steps[0])
    first = sample_promoted_route_actions((evaluation,), sample_size=1)
    second = sample_promoted_route_actions((evaluation,), sample_size=1)
    assert first[0][1].evaluation_id == second[0][1].evaluation_id

    document = render_route_action_review_html((evaluation,), sample_size=1)
    assert "Promoted route-action label review" in document
    assert "broken:S-Cl:SINGLE:NONE" in document
    assert "Accept label" in document
    assert "does not admit an executable template" in document


def test_strictly_admitted_operator_is_not_promoted_for_review() -> None:
    evaluation = _evaluation()
    observed = replace(
        evaluation.steps[0].observed_action,
        operator_roundtrip_verified=True,
        operator_admission_stage="accepted",
        operator_admission_reason=None,
    )
    step = replace(evaluation.steps[0], observed_action=observed)

    assert not is_promoted_route_action(step)
