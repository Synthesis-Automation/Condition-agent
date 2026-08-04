from dataclasses import asdict

import pytest

from reactive_taxonomy import featurize_reaction

from condition_recommender.models import ChemistRankingPreferences
from condition_recommender.preference_scoring import (
    assess_functional_group_tolerance,
    assess_partner_category_similarity,
    reacting_center_categories,
)
from condition_recommender.ranking_preferences import (
    RANKING_COMPONENTS,
    available_ranking_profiles,
    resolve_ranking_preferences,
)


def _core(reaction_smiles: str) -> dict:
    analysis = featurize_reaction(reaction_smiles)
    assert analysis.valid
    assert analysis.reaction_core is not None
    return asdict(analysis.reaction_core)


def test_named_chemist_profiles_are_complete_and_normalized() -> None:
    profiles = available_ranking_profiles()
    assert {value["profile_id"] for value in profiles} >= {
        "default",
        "reactant_category",
        "functional_group_tolerance",
        "evidence",
        "yield",
        "procedure_completeness",
    }
    for value in profiles:
        resolved = resolve_ranking_preferences(
            ChemistRankingPreferences(profile_id=value["profile_id"])
        )
        assert tuple(resolved.weights) == RANKING_COMPONENTS
        assert sum(resolved.weights.values()) == pytest.approx(1.0)


def test_custom_priorities_are_normalized_without_changing_components() -> None:
    weights = {name: 1.0 for name in RANKING_COMPONENTS}
    weights["partner_category"] = 6.0
    resolved = resolve_ranking_preferences(
        ChemistRankingPreferences(profile_id="chemist-custom", weights=weights)
    )
    assert resolved.customized
    assert sum(resolved.weights.values()) == pytest.approx(1.0)
    assert resolved.weights["partner_category"] > resolved.weights["yield"]


def test_reacting_center_categories_distinguish_amine_classes() -> None:
    aliphatic = _core("Brc1ccccc1.CCNCC>>CCN(CC)c1ccccc1")
    aryl = _core("Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1")

    aliphatic_labels = {
        value["label"] for value in reacting_center_categories(aliphatic)
    }
    aryl_labels = {value["label"] for value in reacting_center_categories(aryl)}
    assert "secondary aliphatic amine" in aliphatic_labels
    assert "primary aryl amine" in aryl_labels

    exact_score, exact_evidence = assess_partner_category_similarity(
        aliphatic, aliphatic
    )
    mismatch_score, mismatch_evidence = assess_partner_category_similarity(
        aliphatic, aryl
    )
    assert exact_score == 1.0
    assert mismatch_score is not None and mismatch_score < exact_score
    assert exact_evidence["status"] == "available"
    assert any(
        value["status"] == "category_mismatch"
        for value in mismatch_evidence["matches"]
    )


def test_functional_group_tolerance_requires_direct_precedent_evidence() -> None:
    query = {
        "spectator_groups": (
            {
                "group_id": "ester",
                "chemist_label": "ester",
                "graph_distance": 2,
            },
            {
                "group_id": "nitrile",
                "chemist_label": "nitrile",
                "graph_distance": 4,
            },
        )
    }
    precedents = (
        (
            "reference:one",
            {"spectator_groups": ({"group_id": "ester"},)},
        ),
        ("reference:two", {"spectator_groups": ()}),
    )

    score, evidence = assess_functional_group_tolerance(query, precedents)

    assert score is not None and 0.0 < score < 1.0
    by_group = {value["group_id"]: value for value in evidence["groups"]}
    assert by_group["ester"]["status"] == "directly_demonstrated"
    assert by_group["nitrile"]["status"] == "unknown_not_tolerant"
    assert by_group["nitrile"]["matched_independent_count"] == 0


def test_no_query_spectators_makes_tolerance_not_applicable() -> None:
    score, evidence = assess_functional_group_tolerance(
        {"spectator_groups": ()},
        (("reference:one", {"spectator_groups": ()}),),
    )
    assert score is None
    assert evidence["status"] == "not_applicable"
