"""R-group profiles replace overlapping environment features when available."""

from condition_recommender.signature_features import substituent_profile_similarity
from condition_recommender.similarity import assess_signature_similarity
from reactive_taxonomy import featurize_reaction


def _payload(reaction: str):
    analysis = featurize_reaction(reaction)
    assert analysis.reaction_signature is not None
    assert analysis.reaction_core is not None
    return analysis.to_dict()["reaction_signature"], analysis.to_dict()["reaction_core"]


def test_signature_similarity_uses_port_profiles_when_both_cores_have_them() -> None:
    query, query_core = _payload("CCBr.CN>>CCNC")
    exact, exact_core = _payload("CCBr.CN>>CCNC")
    different, different_core = _payload("c1ccccc1Br.CN>>CNc1ccccc1")

    exact_score = assess_signature_similarity(
        query,
        exact,
        query_reaction_core=query_core,
        precedent_reaction_core=exact_core,
    )
    different_score = assess_signature_similarity(
        query,
        different,
        query_reaction_core=query_core,
        precedent_reaction_core=different_core,
    )

    assert exact_score.components["environment"] == 1.0
    assert different_score.components["environment"] < 1.0
    assert exact_score.score > different_score.score


def test_aromatic_position_changes_reduce_profile_similarity() -> None:
    _, para_core = _payload("Brc1ccc(F)cc1.CN>>CNc1ccc(F)cc1")
    _, meta_core = _payload("Brc1cccc(F)c1.CN>>CNc1cccc(F)c1")

    assert substituent_profile_similarity(para_core, para_core) == 1.0
    positional_score = substituent_profile_similarity(para_core, meta_core)
    assert positional_score is not None
    assert 0.0 < positional_score < 1.0
