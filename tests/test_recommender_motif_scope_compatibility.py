from __future__ import annotations

from chemtools.recommend import recommender as hte


def test_heteroaryl_biaryl_compatibility_map_links_to_aryl_biaryl() -> None:
    compat = hte._load_motif_compatibility_map()
    assert "HeteroAr-Ar" in (compat.get("Ar-Ar") or set())
    assert "Ar-Ar" in (compat.get("HeteroAr-Ar") or set())


def test_motif_token_compatibility_treats_heteroaryl_biaryl_as_aryl_biaryl() -> None:
    assert hte._motif_tokens_compatible("HeteroAr-Ar", "Ar-Ar")
    assert hte._motif_tokens_compatible("Ar-Ar", "HeteroAr-Ar")


def test_transformation_scoring_accepts_heteroaryl_biaryl_query_for_aryl_biaryl_key() -> None:
    recommender = hte.HTERecommender.__new__(hte.HTERecommender)
    db_key = "CRK-v1 |Ar-B(OH)2|HeteroAr-Br -> Ar-Ar | bond_formed: C(ar)-C(ar)"
    score = recommender._score_transformation_match(
        query_motifs={"Ar-B(OH)2", "HeteroAr-Br", "HeteroAr-Ar"},
        db_key=db_key,
        query_reacted={"Ar-B(OH)2", "HeteroAr-Br"},
        query_formed={"HeteroAr-Ar"},
        query_spectators=set(),
    )
    assert score > 0.0
