import pandas as pd

from chemtools.recommend import recommender as hte
from chemtools.recommend.recommender import HTERecommender


def _make_min_hte_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Reaction_Type_Standardized": ["Suzuki_miyaura"],
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
            "Source_Group": ["literature"],
            "spectator_groups": [""],
        }
    )


def _make_source_group_df() -> pd.DataFrame:
    base = _make_min_hte_df()
    rules = base.copy()
    rules["Catalyst"] = ["Ni"]
    rules["Source_Group"] = ["rules"]
    return pd.concat([base, rules], ignore_index=True)


def _make_summary_df() -> pd.DataFrame:
    base = _make_min_hte_df()
    alt = base.copy()
    alt["Reactant_A_Type"] = ["Ar-Cl"]
    alt["Catalyst"] = ["Ni"]
    return pd.concat([base, alt], ignore_index=True)


def _make_weighted_df() -> pd.DataFrame:
    base = _make_min_hte_df()
    base["Reactant_A_Type"] = ["Ar-Cl"]
    base["Reactant_B_Type"] = ["Ar-B(OH)2"]
    base["Catalyst"] = ["Pd"]
    base["reaction_smiles"] = ["rxn_good"]
    alt = base.copy()
    alt["Catalyst"] = ["Ni"]
    alt["reaction_smiles"] = ["rxn_bad"]
    return pd.concat([base, alt], ignore_index=True)


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
    assert "literature" in result.recommendations_by_source
    assert result.total_matching_experiments >= 1
    assert "Ar-X" in (result.matched_motifs or ())


def test_recommend_source_group_filter(monkeypatch) -> None:
    df = _make_source_group_df()
    indexed_data = {"Ar-X": df}
    reaction_type_patterns = {}
    transformation_indices = {}

    def fake_load_db(path: str):
        return df, indexed_data, reaction_type_patterns, transformation_indices

    def fake_detect(self, smiles: str):
        return ["Ar-X"], "Aryl Halide"

    monkeypatch.setattr(hte, "_load_hte_database_cached", fake_load_db)
    monkeypatch.setattr(HTERecommender, "_detect_reactant_types", fake_detect)

    recommender = HTERecommender(hte_db_path="data/HTE_db")
    result = recommender.recommend(
        reactant_a_smiles="Brc1ccccc1",
        top_k=1,
        min_experiments=1,
        source_group="datasets",
    )

    assert result.total_matching_experiments == 1
    assert len(result.recommendations) == 1
    assert result.recommendations[0].catalyst == "Pd"


def test_recommend_top_k_zero_returns_all_matches(monkeypatch) -> None:
    base = _make_min_hte_df()
    alt_1 = base.copy()
    alt_1["Catalyst"] = ["Ni"]
    alt_2 = base.copy()
    alt_2["Catalyst"] = ["Cu"]
    df = pd.concat([base, alt_1, alt_2], ignore_index=True)

    indexed_data = {"Ar-X": df}
    reaction_type_patterns = {}
    transformation_indices = {}

    def fake_load_db(path: str):
        return df, indexed_data, reaction_type_patterns, transformation_indices

    def fake_detect(self, smiles: str):
        return ["Ar-X"], "Aryl Halide"

    def fake_precedent(self, *args, **kwargs):
        return []

    monkeypatch.setattr(hte, "_load_hte_database_cached", fake_load_db)
    monkeypatch.setattr(HTERecommender, "_detect_reactant_types", fake_detect)
    monkeypatch.setattr(HTERecommender, "_build_precedent_recommendations", fake_precedent)

    recommender = HTERecommender(hte_db_path="data/HTE_db")
    result = recommender.recommend(
        reactant_a_smiles="Brc1ccccc1",
        top_k=0,
        min_experiments=1,
    )

    assert result.total_matching_experiments == 3
    assert len(result.recommendations) == 3
    assert set(rec.catalyst for rec in result.recommendations) == {"Pd", "Ni", "Cu"}


def test_summarize_conditions_filters(monkeypatch) -> None:
    df = _make_summary_df()
    indexed_data = {"Ar-X": df}
    reaction_type_patterns = {}
    transformation_indices = {}

    def fake_load_db(path: str):
        return df, indexed_data, reaction_type_patterns, transformation_indices

    monkeypatch.setattr(hte, "_load_hte_database_cached", fake_load_db)
    recommender = HTERecommender(hte_db_path="data/HTE_db")
    payload = recommender.summarize_conditions(
        reaction_type_filter="suzuki_miyaura",
        reactant_type_filters=["Ar-Cl"],
        source_group="literature",
        top_k=3,
        min_experiments=1,
    )

    assert payload["total_matching_experiments"] == 1
    recs = payload["recommendations"]
    assert recs
    assert recs[0].catalyst == "Ni"


def test_aryl_weighting_adjusts_match_score(monkeypatch) -> None:
    df = _make_weighted_df()
    key = hte._reactant_key(["Ar-B(OH)2", "Ar-Cl"])
    indexed_data = {key: df}
    reaction_type_patterns = {}
    transformation_indices = {}

    def fake_load_db(path: str):
        return df, indexed_data, reaction_type_patterns, transformation_indices

    def fake_detect(self, smiles: str):
        if "Cl" in smiles:
            return ["Ar-Cl"], "Aryl Halide"
        return ["Ar-B(OH)2"], "Aryl Boronate"

    def fake_normalize_reaction(rsmi: str):
        mapping = {
            "rxn_good": {
                "reactants": [
                    {"smiles_norm": "ArCl_good"},
                    {"smiles_norm": "ArB_good"},
                ]
            },
            "rxn_bad": {
                "reactants": [
                    {"smiles_norm": "ArCl_bad"},
                    {"smiles_norm": "ArB_good"},
                ]
            },
        }
        return mapping.get(rsmi, {"reactants": []})

    def fake_featurize(smiles: str):
        payload = {
            "motifs": [],
            "steric": {"aryl": [], "alkyl": []},
            "electronics": {"aryl": []},
        }
        if "ArCl" in smiles or "Ar-Cl" in smiles:
            payload["motifs"] = [{"compound_id": "Ar-Cl"}]
            score = 2.0 if "good" in smiles or "query" in smiles else 8.0
            payload["steric"]["aryl"] = [{"result": {"score_0_10": score}}]
            payload["electronics"]["aryl"] = [{"result": {"score_0_10": 6.0 if "good" in smiles or "query" in smiles else 2.0}}]
        elif "ArB" in smiles or "Ar-B" in smiles:
            payload["motifs"] = [{"compound_id": "Ar-B(OH)2"}]
            payload["steric"]["aryl"] = [{"result": {"score_0_10": 1.0}}]
            payload["electronics"]["aryl"] = [{"result": {"score_0_10": 5.0}}]
        return payload

    monkeypatch.setattr(hte, "_load_hte_database_cached", fake_load_db)
    monkeypatch.setattr(HTERecommender, "_detect_reactant_types", fake_detect)
    monkeypatch.setattr(hte, "normalize_reaction", fake_normalize_reaction)
    monkeypatch.setattr(hte, "featurize_molecule", fake_featurize)

    recommender = HTERecommender(hte_db_path="data/HTE_db")
    result = recommender.recommend(
        reactant_a_smiles="ArCl_query",
        reactant_b_smiles="ArB_query",
        top_k=2,
        min_experiments=1,
        use_aryl_steric_electronic_weighting=True,
    )

    assert len(result.recommendations) == 2
    assert result.recommendations[0].catalyst == "Pd"
    assert result.recommendations[0].match_score > result.recommendations[1].match_score


def test_multi_motif_reactant_key_matches(monkeypatch) -> None:
    df = _make_min_hte_df()
    df["Reactant_A_Type"] = ["Ar-Br|Ar-SH"]
    df["Reactant_B_Type"] = ["RCH2-NHR"]

    key = hte._reactant_key([hte._collapse_motif_field("Ar-Br|Ar-SH"), "RCH2-NHR"])
    indexed_data = {key: df}
    reaction_type_patterns = {}
    transformation_indices = {}

    def fake_load_db(path: str):
        return df, indexed_data, reaction_type_patterns, transformation_indices

    def fake_detect(self, smiles: str):
        if "NHR" in smiles:
            return ["RCH2-NHR"], "Amine"
        return ["Ar-Br", "Ar-SH"], "Aryl Halide"

    monkeypatch.setattr(hte, "_load_hte_database_cached", fake_load_db)
    monkeypatch.setattr(HTERecommender, "_detect_reactant_types", fake_detect)

    recommender = HTERecommender(hte_db_path="data/HTE_db")
    result = recommender.recommend(
        reactant_a_smiles="Sc1ccccc1Br",
        reactant_b_smiles="NHR_substrate",
        top_k=1,
        min_experiments=1,
    )

    assert result.total_matching_experiments == 1
    assert result.is_fallback_match is False
    assert result.recommendations[0].catalyst == "Pd"


def test_scaffold_alias_expands_to_aromn_h() -> None:
    motif_sets = hte._load_motif_sets()
    scope_map = hte._load_scope_map()
    expanded = hte._expand_reactant_field("Pyrazole", motif_sets, scope_map)
    assert "Pyrazole" in expanded
    assert "AromN-H" in expanded


def test_exact_match_priority_beats_fallback_alias(monkeypatch) -> None:
    df = _make_min_hte_df()
    df["Reactant_A_Type"] = ["Ar-Br"]
    df["Reactant_B_Type"] = ["Ar-B(OH)2"]
    df["Source_Group"] = ["rules"]
    df["Reactant_Signature_Core"] = ["Ar-B(OH)2|Ar-Br"]
    df["Reactant_Signature_Ext"] = ["Ar-B(OH)2|Ar-Br"]

    direct_key = hte._reactant_key(["Ar-Br", "Ar-B(OH)2"])
    fallback_key = hte._reactant_key(["Ar-X", "Ar-B(OH)2"])
    indexed_data = {
        direct_key: df.iloc[[0]].copy(),
        fallback_key: df.iloc[[0]].copy(),
    }
    reaction_type_patterns = {direct_key: hte.Counter({"Suzuki_miyaura": 1})}
    transformation_indices = {}

    def fake_load_db(path: str):
        return df, indexed_data, reaction_type_patterns, transformation_indices

    def fake_detect(self, smiles: str):
        if "B(" in smiles:
            return ["Ar-B(OH)2"], "Aryl Boronate"
        return ["Ar-Br", "Ar-X"], "Aryl Halide"

    monkeypatch.setattr(hte, "_load_hte_database_cached", fake_load_db)
    monkeypatch.setattr(HTERecommender, "_detect_reactant_types", fake_detect)

    recommender = HTERecommender(hte_db_path="data/HTE_db")
    result = recommender.recommend(
        reactant_a_smiles="Brc1ccccc1",
        reactant_b_smiles="OB(O)c1ccccc1",
        top_k=1,
        min_experiments=2,
    )

    assert result.recommendations
    assert result.recommendations[0].match_score == 1.0


def test_key_match_backfills_rules_and_experiments_with_aromatic_alias(monkeypatch) -> None:
    df = pd.DataFrame(
        {
            "Reaction_Type_Standardized": ["Suzuki_miyaura", "Suzuki_miyaura", "Suzuki_miyaura"],
            "Reactant_A_Type": ["HeteroAr-Br", "Ar-Br", "Ar-Br"],
            "Reactant_B_Type": ["Ar-B(OH)2", "Ar-B(OH)2", "Ar-B(OH)2"],
            "Reactant_C_Type": ["", "", ""],
            "Catalyst": ["Pd", "Ni", "Pd"],
            "Ligand": ["PPh3", "dtbpy", "SPhos"],
            "Base": ["K2CO3", "Cs2CO3", "K3PO4"],
            "Solvent": ["DMF", "DMAc", "THF"],
            "Secondary Solvent": ["", "", ""],
            "Additive": ["", "", ""],
            "Coupling Reagent": ["", "", ""],
            "AREA_TOTAL_REDUCED": [85.0, 70.0, 78.0],
            "z-Score": [2.2, 1.3, 1.6],
            "Reactant_A_Category": ["", "", ""],
            "Reactant_B_Category": ["", "", ""],
            "Reaction_Category": ["", "", ""],
            "Is_Intramolecular": [False, False, False],
            "Source_File": ["tests", "tests", "tests"],
            "Source_Group": ["literature", "rules", "experiments"],
            "spectator_groups": ["", "", ""],
            "Reactant_Signature_Core": [
                "Ar-B(OH)2|HeteroAr-Br",
                "Ar-B(OH)2|Ar-Br",
                "Ar-B(OH)2|Ar-Br",
            ],
            "Reactant_Signature_Ext": [
                "Ar-B(OH)2|HeteroAr-Br",
                "Ar-B(OH)2|Ar-Br",
                "Ar-B(OH)2|Ar-Br",
            ],
            "Reaction_Key": [
                "CRK-v1 |Ar-B(OH)2|HeteroAr-Br -> Ar-Ar | bond_formed: C(ar)-C(ar)",
                "",
                "",
            ],
        }
    )

    transformation_indices = {
        "Ar-B(OH)2|HeteroAr-Br": df.iloc[[0]].copy(),
    }
    indexed_data = {
        hte._reactant_key(["Ar-B(OH)2", "Ar-Br"]): df.iloc[[1, 2]].copy(),
    }
    reaction_type_patterns = {
        hte._reactant_key(["Ar-B(OH)2", "Ar-Br"]): hte.Counter({"Suzuki_miyaura": 2}),
    }

    def fake_load_db(path: str):
        return df, indexed_data, reaction_type_patterns, transformation_indices

    def fake_detect(self, smiles: str):
        return [], "Unknown"

    def fake_featurize_reaction(smiles: str, options=None):
        return {
            "reaction_type": "Suzuki_miyaura",
            "reaction_key": "CRK-v1 |Ar-B(OH)2|HeteroAr-Br -> Ar-Ar | bond_formed: C(ar)-C(ar)",
            "aggregates": {
                "reacted_motifs": ["Ar-B(OH)2", "HeteroAr-Br"],
                "formed_motifs": ["Ar-Ar"],
                "spectator_motifs": [],
            },
        }

    monkeypatch.setattr(hte, "_load_hte_database_cached", fake_load_db)
    monkeypatch.setattr(HTERecommender, "_detect_reactant_types", fake_detect)
    monkeypatch.setattr(hte, "featurize_reaction", fake_featurize_reaction)

    recommender = HTERecommender(hte_db_path="data/HTE_db")
    result = recommender.recommend(
        reactant_a_smiles="Brc1cccnc1",
        reactant_b_smiles="Nc1ccccc1B(O)O",
        product_smiles="Nc1ccccc1-c1cccnc1",
        top_k=3,
        min_experiments=1,
    )

    assert result.total_matching_experiments >= 3
    assert "rules" in result.recommendations_by_source
    assert "experiments" in result.recommendations_by_source
