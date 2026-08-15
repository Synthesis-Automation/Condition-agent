"""Observed route-action replay schema and matching regressions."""

from __future__ import annotations

from core_retrosynthesis.generic_compiler import (
    analyze_generic_reaction,
    compile_generic_templates,
)
from core_retrosynthesis.generic_models import (
    GenericDisconnectionCandidate,
    GenericTemplateLibrary,
)
from core_retrosynthesis.route_action_evaluation import (
    RouteActionEvaluation,
    RouteActionEvaluationConfig,
    evaluate_route_actions,
    normalize_observed_reaction,
)
from core_retrosynthesis.route_conversion import build_observed_route_tree


REACTION = "O=Cc1ccccc1>>OCc1ccccc1"
PRODUCT = "OCc1ccccc1"
PRECURSOR = "O=Cc1ccccc1"


def _record(reaction: str = REACTION) -> dict:
    return {
        "schema_version": "1.0",
        "route_id": "US00000001A1_0",
        "patent_id": "US00000001A1",
        "split": "test",
        "target_smiles": PRODUCT,
        "original_reaction_count": 1,
        "higher_level_reaction_count": 1,
        "higher_level_depth": 1,
        "abstraction_reduction": 0,
        "steps": [
            {
                "source_reaction_id": "reaction-1",
                "product_smiles": PRODUCT,
                "precursor_smiles": [PRECURSOR],
                "reactants_smiles": reaction.split(">", 1)[0],
                "reagents_smiles": "",
                "product_smiles_mapped": reaction.rsplit(">", 1)[-1],
                "reaction_smiles": reaction,
                "abstracted_reaction_smiles": "",
            }
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


def _candidate(
    precursor: str,
    *,
    operator: str,
    site: str,
    synthon: str,
    precedent: str,
) -> GenericDisconnectionCandidate:
    return GenericDisconnectionCandidate(
        target_smiles=PRODUCT,
        precursor_smiles=precursor,
        proposed_reaction_smiles=f"{precursor}>>{PRODUCT}",
        transformation_kind="carbonyl_reduction",
        abstraction_level="L2",
        compiler_engine="test",
        template_id=f"template-{precursor}",
        score=0.8,
        context_similarity=0.8,
        product_similarity=0.8,
        precursor_similarity=0.8,
        template_specificity=0.8,
        independent_reference_support=1,
        forward_validation_status="verified_signature",
        center_transition_key="center",
        disconnection_site_key=site,
        precedent_reaction_ids=(precedent,),
        operator_id=operator,
        realization_id=f"realization-{precursor}",
        operator_signature=f"signature-{operator}",
        synthon_signature=synthon,
    )


def test_normalizes_three_part_reaction_without_promoting_reagents() -> None:
    assert normalize_observed_reaction("CCO>O.[Na+]>CC=O") == "CCO>>CC=O"
    assert normalize_observed_reaction("CCO>>CC=O") == "CCO>>CC=O"
    assert normalize_observed_reaction("CCO") is None


def test_route_replay_ranks_distinct_strategy_and_exact_realization() -> None:
    compiled = compile_generic_templates(
        {"reaction_smiles": REACTION, "reaction_id": "expected"},
        levels=("L0", "L1", "L2"),
        admission_mode="data_driven",
    )
    template = compiled.templates[0]
    identity = analyze_generic_reaction(template.precedents[0].mapped_reaction_smiles)
    assert identity is not None

    def searcher(*args: object, **kwargs: object) -> tuple[tuple[object, ...], dict]:
        return (
            (
                _candidate(
                    "CC",
                    operator="OP1:other",
                    site="SITE1:other",
                    synthon="SYN1:other",
                    precedent="US99999999A1:row-1",
                ),
                _candidate(
                    "CCO",
                    operator=template.operator_id,
                    site=identity.disconnection_site_key,
                    synthon=template.synthon_signature,
                    precedent="US00000001A1:row-2",
                ),
                _candidate(
                    PRECURSOR,
                    operator=template.operator_id,
                    site=identity.disconnection_site_key,
                    synthon=template.synthon_signature,
                    precedent="US00000001A1:row-3",
                ),
            ),
            {"proposed_action_count": 3, "validation_attempt_count": 3},
        )

    evaluation = evaluate_route_actions(
        build_observed_route_tree(_record()),
        _empty_library(),
        config=RouteActionEvaluationConfig(use_hierarchical_ranking=False),
        searcher=searcher,
    )

    step = evaluation.steps[0]
    assert step.observed_action.strategy_verified
    assert step.observed_action.operator_roundtrip_verified
    assert step.search_status == "searched"
    assert step.exact_precursor_rank == 3
    assert step.site_rank == 2
    assert step.operator_rank == 2
    assert step.synthon_rank == 2
    assert step.strategy_rank == 2
    assert step.outcome == "exact_lower_rank"
    assert step.source_patent_precedent_overlap
    assert step.candidates[1].strategy_match
    assert not step.candidates[1].exact_precursor_match
    assert step.candidates[1].supervision_label == "strategy_equivalent"
    assert step.candidates[2].supervision_label == "observed_exact"
    assert RouteActionEvaluation.from_dict(evaluation.to_dict()) == evaluation


def test_unreconstructable_observation_is_retained_as_ineligible() -> None:
    record = _record("[CH3:1][Cl:3]>>[CH3:1][NH2:3]")
    record["target_smiles"] = "CN"
    record["steps"][0]["product_smiles"] = "CN"
    tree = build_observed_route_tree(record)

    evaluation = evaluate_route_actions(
        tree,
        _empty_library(),
        config=RouteActionEvaluationConfig(use_hierarchical_ranking=False),
    )

    step = evaluation.steps[0]
    assert not step.observed_action.search_eligible
    assert "reaction_core_unavailable" in step.observed_action.limitations
    assert step.outcome == "not_search_eligible"
    assert step.candidates == ()
