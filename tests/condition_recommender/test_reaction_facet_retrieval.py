"""Regressions for graph-derived progressive reaction-facet retrieval."""

from dataclasses import asdict

from reactive_taxonomy import featurize_reaction

from condition_recommender.generic_api import recommend_indexed_signature
from condition_recommender.generic_indexing import build_generic_index
from condition_recommender.reaction_facets import (
    load_reaction_facet_rules,
    reaction_facet_keys,
)


QUERY = "Clc1ccncc1.NC(=O)c1ccccc1>>O=C(Nc1ccncc1)c1ccccc1"
ARYL_PRIMARY_AMIDE = (
    "Clc1ccccc1.NC(=O)c1ccccc1>>O=C(Nc1ccccc1)c1ccccc1"
)
HETARYL_SECONDARY_AMIDE = (
    "Clc1ccncc1.CNC(=O)c1ccccc1>>CN(c1ccncc1)C(=O)c1ccccc1"
)
HETARYL_CARBAMATE = (
    "Clc1ccncc1.NC(=O)OC(C)(C)C>>CC(C)(C)OC(=O)Nc1ccncc1"
)
SUZUKI = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"


def _analysis(reaction_smiles: str):
    analysis = featurize_reaction(reaction_smiles)
    assert analysis.valid
    assert analysis.reaction_signature is not None
    assert analysis.reaction_core is not None
    assert analysis.fallback_descriptor is not None
    return analysis


def _keys(reaction_smiles: str) -> dict[str, str]:
    analysis = _analysis(reaction_smiles)
    return reaction_facet_keys(
        asdict(analysis.reaction_signature),
        asdict(analysis.reaction_core),
        asdict(analysis.fallback_descriptor),
    )


def _record(index: int, reaction_smiles: str, recipe_core: str) -> dict:
    analysis = _analysis(reaction_smiles)
    recipe_id = recipe_core.replace("RCORE1", "RCR1")
    return {
        "schema_version": "10.0",
        "converter_definition_version": "generic_conversion.v10.0",
        "admission_tier": "verified",
        "index_eligibility": "eligible",
        "precedent_tier": "trusted",
        "core_eligibility": "trusted_core",
        "core_eligibility_definition_version": "core_eligibility.v1@1.0",
        "chemistry_status": "verified",
        "condition_status": "resolved_complete",
        "condition_stage_status": "single_stage",
        "outcome_status": "usable",
        "reaction_id": f"reaction-{index}",
        "observation_id": f"observation-{index}",
        "canonical_reaction_id": f"CRX1:{index}",
        "reaction_smiles": reaction_smiles,
        "reaction_label": asdict(analysis.reaction_label),
        "yield_pct": 70.0 + index,
        "source_dataset": "facet-test",
        "reference_id": f"REF1:{index}",
        "reference_condition_series_id": f"RCS1:{index}",
        "reaction_signature": asdict(analysis.reaction_signature),
        "reaction_core": asdict(analysis.reaction_core),
        "fallback_descriptor": asdict(analysis.fallback_descriptor),
        "resolved_recipe_id": recipe_id,
        "resolved_recipe_core_id": recipe_core,
        "resolved_recipe": {
            "recipe_id": recipe_id,
            "recipe_core_id": recipe_core,
        },
        "condition_resolution": {"has_uncertainty": False},
    }


def test_facet_definition_and_projection_are_type_agnostic() -> None:
    rules = load_reaction_facet_rules()
    assert rules["retrieval_ladder"] == (
        "reaction_facet_exact",
        "reaction_facet_attachment_relaxed",
    )
    assert set(_keys(QUERY)) == set(rules["retrieval_ladder"])
    assert set(_keys(SUZUKI)) == set(rules["retrieval_ladder"])


def test_attachment_relaxation_preserves_active_atom_chemistry() -> None:
    query = _keys(QUERY)
    aryl = _keys(ARYL_PRIMARY_AMIDE)
    secondary = _keys(HETARYL_SECONDARY_AMIDE)
    carbamate = _keys(HETARYL_CARBAMATE)

    assert query["reaction_facet_exact"] != aryl["reaction_facet_exact"]
    assert (
        query["reaction_facet_attachment_relaxed"]
        == aryl["reaction_facet_attachment_relaxed"]
    )
    assert (
        query["reaction_facet_attachment_relaxed"]
        != secondary["reaction_facet_attachment_relaxed"]
    )
    assert (
        query["reaction_facet_attachment_relaxed"]
        != carbamate["reaction_facet_attachment_relaxed"]
    )


def test_progressive_retrieval_fills_distinct_recipes_by_tier() -> None:
    query = _analysis(QUERY)
    index = build_generic_index(
        [
            _record(1, QUERY, "RCORE1:exact-a"),
            _record(2, QUERY, "RCORE1:exact-b"),
            _record(3, ARYL_PRIMARY_AMIDE, "RCORE1:relaxed"),
            _record(4, HETARYL_CARBAMATE, "RCORE1:carbamate"),
        ]
    )

    result = recommend_indexed_signature(
        asdict(query.reaction_signature),
        index,
        reaction_core=asdict(query.reaction_core),
        fallback_descriptor=asdict(query.fallback_descriptor),
        top_k=3,
        minimum_pool_size=1,
    )

    assert result.valid
    assert result.retrieval_level == "progressive_reaction_facets"
    assert result.retrieval_definition_version.startswith(
        "reaction_facet_retrieval.v1@1.0;"
    )
    assert [value.recipe_core_id for value in result.recommendations] == [
        "RCORE1:exact-b",
        "RCORE1:exact-a",
        "RCORE1:relaxed",
    ]
    assert [value.retrieval_level for value in result.recommendations] == [
        "reaction_facet_exact",
        "reaction_facet_exact",
        "reaction_facet_attachment_relaxed",
    ]
    assert "RCORE1:carbamate" not in {
        value.recipe_core_id for value in result.recommendations
    }
    assert tuple(value.level for value in result.retrieval_trace) == (
        "reaction_facet_exact",
        "reaction_facet_attachment_relaxed",
    )
    assert result.retrieval_trace[-1].status == "selected_target_reached"
    assert "REACTION_FACET_RELAXATION_USED" in result.warnings


def test_top_k_stops_before_relaxation_when_exact_recipes_are_sufficient() -> None:
    query = _analysis(QUERY)
    index = build_generic_index(
        [
            _record(1, QUERY, "RCORE1:exact-a"),
            _record(2, QUERY, "RCORE1:exact-b"),
            _record(3, ARYL_PRIMARY_AMIDE, "RCORE1:relaxed"),
        ]
    )

    result = recommend_indexed_signature(
        asdict(query.reaction_signature),
        index,
        reaction_core=asdict(query.reaction_core),
        fallback_descriptor=asdict(query.fallback_descriptor),
        top_k=2,
        minimum_pool_size=1,
    )

    assert result.retrieval_level == "reaction_facet_exact"
    assert len(result.recommendations) == 2
    assert tuple(value.level for value in result.retrieval_trace) == (
        "reaction_facet_exact",
    )
    assert "REACTION_FACET_RELAXATION_USED" not in result.warnings
