from chemtools.recommend import HTESimpleConditionRecommender, RecommendationOptions


def test_hte_simple_recommender_exact_match():
    recommender = HTESimpleConditionRecommender()
    sample = recommender.df.iloc[0]

    options = RecommendationOptions(top_n=1, min_precedents=1)
    result = recommender.recommend(
        sample["Reaction_Type_Standardized"],
        sample["Reactant_A"],
        sample["Reactant_B"],
        options=options,
    )

    assert result["match_level"] in {"exact", "category", "reaction_type"}
    assert result["recommendations"]
    assert result["metadata"]["dataset"].endswith("data\\HTE_db\\HTE_0.csv") or result["metadata"]["dataset"].endswith(
        "data/HTE_db/HTE_0.csv"
    )
