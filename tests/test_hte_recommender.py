import pandas as pd

from chemtools.HTE import recommender as hte
from chemtools.HTE.recommender import HTERecommender


def _make_min_hte_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Reaction_Type_Standardized": ["suzuki_miyaura"],
            "Reactant_A_Type": ["Ar-X"],
            "Reactant_B_Type": ["Ar-X"],
            "Reactant_C_Type": [""],
            "Catalyst": ["Pd"],
            "Ligand": ["PPh3"],
            "Base": ["K2CO3"],
            "Solvent": ["DMF"],
            "Secondary Solvent": [""],
            "Additive": [""],
            "Coupling Reagent": [""],
            "AREA_TOTAL_REDUCED": [80.0],
            "z-Score": [2.0],
            "Reactant_A_Category": [""],
            "Reactant_B_Category": [""],
            "Reaction_Category": [""],
            "Is_Intramolecular": [False],
            "Source_File": ["tests"],
            "Source_Group": ["datasets"],
            "spectator_groups": [""],
        }
    )


def test_fallback_tiers_order() -> None:
    tiers = hte._build_fallback_tiers(["Ar-Br", "Ar-X", "Ar-R"], set(), set())
    assert tiers
    assert "Ar-Br" in tiers[0]
    assert "Ar-R" not in tiers[0]
    if len(tiers) > 1:
        assert "Ar-R" in tiers[1]


def test_recommend_falls_back_when_direct_key_missing(monkeypatch) -> None:
    df = _make_min_hte_df()
    indexed_data = {"Ar-X": df}
    reaction_type_patterns = {}
    transformation_indices = {}

    def fake_load_db(path: str):
        return df, indexed_data, reaction_type_patterns, transformation_indices

    def fake_detect(self, smiles: str):
        if "Br" in smiles:
            return ["Ar-Br", "Ar-X"], "Aryl Halide"
        return ["Ar-Cl", "Ar-X"], "Aryl Halide"

    monkeypatch.setattr(hte, "_load_hte_database_cached", fake_load_db)
    monkeypatch.setattr(HTERecommender, "_detect_reactant_types", fake_detect)

    recommender = HTERecommender(hte_db_path="data/HTE_db")
    result = recommender.recommend(
        reactant_a_smiles="Brc1ccccc1",
        reactant_b_smiles="Clc1ccccc1",
        top_k=1,
        min_experiments=1,
    )

    assert result.is_fallback_match is True
    assert len(result.recommendations) == 1
    assert "datasets" in result.recommendations_by_source
    assert result.total_matching_experiments >= 1
    assert "Ar-X" in (result.matched_motifs or ())
