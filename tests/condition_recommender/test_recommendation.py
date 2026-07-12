from pathlib import Path

from condition_recommender import recommend_conditions
from condition_recommender.indexing import IndexedReaction
from condition_recommender.models import ConditionIdentity
from condition_recommender.retrieval import load_retrieval_rules, retrieve_pool, structured_similarity


def _indexed(reaction_id: str, e_context: str = "Ar", t_context: str = "Ar") -> IndexedReaction:
    return IndexedReaction(
        reaction_id=reaction_id,
        reaction_smiles="",
        yield_pct=80.0,
        conditions=ConditionIdentity(("14221-01-3",), ("584-08-7",), ("108-88-3",)),
        recipe_id="recipe-1",
        electrophile={"handle_token": "Br", "anchor_context": e_context, "steric": {"class": "open"}, "electronic": {"class": "neutral"}},
        transfer_partner={"handle_token": "B(OH)2", "anchor_context": t_context, "steric": {"class": "open"}, "electronic": {"class": "neutral"}},
        spectator_group_ids=(),
        family_flags=(),
        product_connection_label="Ar–Ar",
    )


def _query() -> dict:
    return {
        "electrophile": {"handle_token": "Br", "anchor_context": "Ar", "steric": {"class": "open"}, "electronic": {"class": "neutral"}},
        "transfer_partner": {"handle_token": "B(OH)2", "anchor_context": "Ar", "steric": {"class": "open"}, "electronic": {"class": "neutral"}},
        "spectator_group_ids": (),
        "family_flags": (),
    }


def test_structured_similarity_rewards_exact_partner_signature() -> None:
    exact, _ = structured_similarity(_query(), _indexed("exact"))
    different, _ = structured_similarity(_query(), _indexed("different", "HeteroAr", "HeteroAr"))
    assert exact == 1.0
    assert exact > different


def test_retrieval_weights_are_normalized() -> None:
    rules = load_retrieval_rules()
    assert round(sum(rules["similarity_weights"].values()), 10) == 1.0
    assert round(sum(rules["ranking_weights"].values()), 10) == 1.0


def test_retrieval_uses_exact_signature_when_supported() -> None:
    rows = tuple(_indexed(str(index)) for index in range(12)) + (_indexed("other", "HeteroAr"),)
    level, pool = retrieve_pool(_query(), rows)
    assert level == "exact_partner_signature"
    assert len(pool) == 12


def test_pilot_recommendation_returns_ranked_recipes() -> None:
    path = Path("results/suzuki_recommendation_pilot/verified.csv")
    result = recommend_conditions(
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        records_path=path,
        top_k=3,
    )
    assert result.valid
    assert result.product_connection_label == "Ar–Ar"
    assert len(result.recommendations) == 3
    assert [item.rank for item in result.recommendations] == [1, 2, 3]
    assert all(item.conditions.complete for item in result.recommendations)
    assert all(0.0 <= item.score <= 1.0 for item in result.recommendations)
    assert all(0.0 <= item.expected_yield_pct <= 100.0 for item in result.recommendations)


def test_unverified_query_is_rejected() -> None:
    result = recommend_conditions(
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1",
        top_k=3,
    )
    assert not result.valid
    assert result.error == "QUERY_REACTION_NOT_EXACTLY_VERIFIED"
