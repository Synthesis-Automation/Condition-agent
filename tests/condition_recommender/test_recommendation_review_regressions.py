"""Reproduced ranking, constraint and condition-coverage failures."""

from copy import deepcopy
from dataclasses import asdict
from functools import lru_cache


from condition_registry import ConditionConstraintSet, normalize_condition_constraint
from condition_recommender.compatibility import assess_recipe_compatibility
from condition_recommender.conversion.generic import convert_record
from condition_recommender.conversion.input_schema import adapt_row
from condition_recommender.generic_api import GenericConditionRecommender
from condition_recommender.generic_indexing import build_generic_index
from condition_recommender.generic_retrieval import retrieve_generic_pool_with_trace
from condition_recommender.recipe_ranking import rank_condition_recipes
from reactive_taxonomy import featurize_reaction


REACTION = "[CH3:1][Br:2].[NH2:3][CH3:4]>>[CH3:1][NH:3][CH3:4]"


@lru_cache(maxsize=1)
def _base_record():
    return convert_record(
        adapt_row(
            {
                "reaction_id": "condition-review",
                "reaction_smiles": REACTION,
                "reagent_cas": "584-08-7",
                "solvent_cas": "64-17-5",
                "yield_pct": "80",
                "reference": "development-review",
            },
            source_dataset="development",
            source_path="development.csv",
            source_row_number=2,
        )
    ).to_dict()


def _record(
    number, temperature=25, yield_pct=70, core="same-ingredients", reference=None
):
    value = deepcopy(_base_record())
    recipe = value["resolved_recipe"]
    recipe.update(
        recipe_id=f"{core}:{temperature}",
        recipe_core_id=core,
        temperature_c=temperature,
    )
    value.update(
        reaction_id=f"r-{number:03}",
        observation_id=f"obs-{number:03}",
        reference_id=reference or f"ref-{number}",
        resolved_recipe_id=recipe["recipe_id"],
        resolved_recipe_core_id=core,
        yield_pct=yield_pct,
    )
    return value


def _rank(records, top_k=5):
    index = build_generic_index(records)
    signature = _base_record()["reaction_signature"]
    pool = tuple(
        (row, assess_recipe_compatibility(signature, row.resolved_recipe))
        for row in index.rows
    )
    return rank_condition_recipes(
        signature,
        pool,
        retrieval_level="exact_signature",
        top_k=top_k,
        query_reaction_smiles=REACTION,
    )


def test_constraints_select_surviving_variant_before_aggregation():
    index = build_generic_index([_record(1, 120, 90), _record(2, 25, 40)])
    constraint = normalize_condition_constraint(
        "maximum_temperature_c", 50, provenance="explicit_user"
    ).constraint
    result = GenericConditionRecommender(index).recommend(
        REACTION,
        top_k=1,
        minimum_pool_size=1,
        condition_constraints=ConditionConstraintSet((constraint,)),
    )
    assert result.valid and len(result.recommendations) == 1
    assert result.recommendations[0].resolved_recipe["temperature_c"] == 25
    assert result.recommendations[0].historical_yield_pct == 40
    assert result.recommendations[0].observation_support == 1


def test_constraints_do_not_stop_at_twenty_one_forbidden_recipes():
    records = [_record(i, 120, 90, core=f"hot-{i}") for i in range(25)] + [
        _record(99, 25, 40, core="cold")
    ]
    index = build_generic_index(records)
    constraint = normalize_condition_constraint(
        "maximum_temperature_c", 50, provenance="explicit_user"
    ).constraint
    result = GenericConditionRecommender(index).recommend(
        REACTION,
        top_k=1,
        minimum_pool_size=1,
        condition_constraints=ConditionConstraintSet((constraint,)),
    )
    assert result.valid and result.recommendations[0].recipe_core_id == "cold"


def test_direct_pool_retrieval_filters_constraints_before_support():
    index = build_generic_index([_record(1, 120), _record(2, 25)])
    constraint = normalize_condition_constraint(
        "maximum_temperature_c", 50, provenance="explicit_user"
    ).constraint
    _, rows, trace = retrieve_generic_pool_with_trace(
        _base_record()["reaction_signature"],
        index,
        minimum_pool_size=1,
        condition_constraints=ConditionConstraintSet((constraint,)),
    )
    assert len(rows) == 1 and rows[0].resolved_recipe["temperature_c"] == 25
    assert any(item.excluded_candidate_count == 1 for item in trace)


def test_neighbor_budget_is_per_recipe_and_preserves_total_support():
    records = [_record(i, core="A") for i in range(55)] + [_record(99, core="B")]
    recommendations = _rank(records)
    assert {value.recipe_core_id for value in recommendations} == {"A", "B"}
    assert (
        next(
            value for value in recommendations if value.recipe_core_id == "A"
        ).reference_support
        == 55
    )


def test_historical_yield_includes_low_outcomes_and_balances_references():
    values = [
        _record(1, yield_pct=10, reference="same"),
        _record(2, yield_pct=90, reference="same"),
    ]
    recommendation = _rank(values)[0]
    assert recommendation.historical_yield_pct == 50
    assert recommendation.historical_yield_summary["observation_count"] == 2
    assert recommendation.historical_yield_summary["minimum_pct"] == 10
    assert recommendation.historical_yield_summary["maximum_pct"] == 90
    assert recommendation.historical_yield_summary["is_prediction"] is False
    assert (
        _rank(values + [_record(3, yield_pct=80, reference="independent")])[
            0
        ].historical_yield_pct
        == 65
    )


def test_selected_variant_does_not_borrow_other_temperature_yields():
    recommendation = _rank([_record(1, 25, 20), _record(2, 120, 90)])[0]
    expected = {25: 20, 120: 90}[recommendation.resolved_recipe["temperature_c"]]
    assert recommendation.historical_yield_pct == expected
    assert recommendation.historical_yield_summary["observation_count"] == 1


def test_missing_conditions_have_unknown_coverage_without_becoming_hard_conflict():
    assessment = assess_recipe_compatibility(_base_record()["reaction_signature"], {})
    assert assessment.compatible and assessment.status == "unknown"
    assert assessment.score == 0 and assessment.evidence


def test_hydrogen_source_check_is_structure_backed_and_preserves_ambiguity():
    signature = asdict(
        featurize_reaction("[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]").reaction_signature
    )
    assert signature["named_family"] is None
    recipe = deepcopy(_base_record()["resolved_recipe"])
    absent = assess_recipe_compatibility(signature, recipe)
    assert absent.status == "unknown" and absent.unresolved_requirements
    supplied = assess_recipe_compatibility(
        signature, {**recipe, "atmosphere": "hydrogen"}
    )
    assert supplied.compatible and not supplied.unresolved_requirements
    assert supplied.checked_requirements
    assert assess_recipe_compatibility(
        signature, {**recipe, "atmosphere": "hydrogen-free"}
    ).unresolved_requirements
    unrelated = {**signature, "order_changes": ["C-O:SINGLE>DOUBLE"]}
    assert not assess_recipe_compatibility(unrelated, recipe).checked_requirements
