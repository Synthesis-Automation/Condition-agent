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


def test_recommend_can_match_by_reaction_events_key(monkeypatch) -> None:
    df = _make_min_hte_df()
    df["Reaction_Type_Standardized"] = ["C_N_Coupling"]
    df["Reaction_Key"] = [""]
    df["Reaction_Events"] = ["sig:LGDisp+C-N | form:C(ar)-N | break:Br-C(ar)"]

    event_key = hte._reaction_events_to_match_key(df.loc[0, "Reaction_Events"])
    indexed_data = {}
    reaction_type_patterns = {}
    transformation_indices = {event_key: df.copy()}

    def fake_load_db(path: str):
        return df, indexed_data, reaction_type_patterns, transformation_indices

    def fake_detect(self, smiles: str):
        return [], "Unknown"

    def fake_featurize_reaction(smiles: str, options=None):
        return {
            "reaction_type": {"reaction_type": "C_N_Coupling", "confidence": 0.95},
            "reaction_key": "",
            "reaction_events": {
                "events": [
                    {"kind": "c_n_bond_formation", "confidence": 0.9},
                    {"kind": "leaving_group_displacement", "confidence": 0.92},
                ],
                "bond_pairs": {"formed": [("C", "N")], "broken": [("Br", "C")]},
            },
            "aggregates": {
                "reacted_motifs": ["Ar-Br", "HeteroAr-NH2"],
                "formed_motifs": ["HeteroAr-NHR"],
                "spectator_motifs": [],
            },
        }

    monkeypatch.setattr(hte, "_load_hte_database_cached", fake_load_db)
    monkeypatch.setattr(HTERecommender, "_detect_reactant_types", fake_detect)
    monkeypatch.setattr(hte, "featurize_reaction", fake_featurize_reaction)

    recommender = HTERecommender(hte_db_path="data/HTE_db")
    result = recommender.recommend(
        reactant_a_smiles="Brc1ccccc1",
        reactant_b_smiles="Nc1ncccn1",
        product_smiles="c1ccc(Nc2ncccn2)cc1",
        top_k=1,
        min_experiments=1,
    )

    assert result.total_matching_experiments == 1
    assert result.recommendations
    assert result.query_reaction_events_key == event_key


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


def _make_mixed_reaction_type_key_df() -> tuple[pd.DataFrame, str]:
    query_key = "CRK-v1 |Ar-B(OH)2|Ar-Cl -> Ar-Ar | bond_formed: C(ar)-C(ar)"
    suzuki = _make_min_hte_df()
    suzuki["Reaction_Type_Standardized"] = ["Suzuki_miyaura"]
    suzuki["Catalyst"] = ["Pd"]
    suzuki["z-Score"] = [1.5]
    suzuki["AREA_TOTAL_REDUCED"] = [78.0]
    suzuki["Reaction_Key"] = [query_key]
    suzuki["Reactant_Signature_Core"] = ["Ar-B(OH)2|Ar-Cl"]
    suzuki["Reactant_Signature_Ext"] = ["Ar-B(OH)2|Ar-Cl"]

    chan_lam = suzuki.copy()
    chan_lam["Reaction_Type_Standardized"] = ["ChanLam_dataset_converted (2)"]
    chan_lam["Catalyst"] = ["Cu"]
    chan_lam["z-Score"] = [4.2]
    chan_lam["AREA_TOTAL_REDUCED"] = [89.0]

    df = pd.concat([chan_lam, suzuki], ignore_index=True)
    return df, query_key


def test_reaction_type_filter_blocks_cross_family_key_matches(monkeypatch) -> None:
    df, query_key = _make_mixed_reaction_type_key_df()
    indexed_data = {}
    reaction_type_patterns = {}
    transformation_indices = {query_key: df.copy()}

    def fake_load_db(path: str):
        return df, indexed_data, reaction_type_patterns, transformation_indices

    def fake_detect(self, smiles: str):
        return [], "Unknown"

    def fake_featurize_reaction(smiles: str, options=None):
        return {
            "reaction_type": {"reaction_type": "Suzuki_miyaura", "confidence": 0.95},
            "reaction_key": query_key,
            "aggregates": {
                "reacted_motifs": ["Ar-B(OH)2", "Ar-Cl"],
                "formed_motifs": ["Ar-Ar"],
                "spectator_motifs": [],
            },
        }

    def fake_precedent(self, *args, **kwargs):
        return []

    monkeypatch.setattr(hte, "_load_hte_database_cached", fake_load_db)
    monkeypatch.setattr(HTERecommender, "_detect_reactant_types", fake_detect)
    monkeypatch.setattr(hte, "featurize_reaction", fake_featurize_reaction)
    monkeypatch.setattr(HTERecommender, "_build_precedent_recommendations", fake_precedent)

    recommender = HTERecommender(hte_db_path="data/HTE_db")
    result = recommender.recommend(
        reactant_a_smiles="Clc1ccccc1",
        reactant_b_smiles="B(O)Oc1ccccc1",
        product_smiles="c1ccccc1-c1ccccc1",
        top_k=5,
        min_experiments=1,
        reaction_type_filter="Suzuki",
        source_group="literature",
    )

    assert result.recommendations
    assert result.total_matching_experiments == 1
    assert all(rec.reaction_type == "Suzuki_miyaura" for rec in result.recommendations)


def test_detected_type_filter_applies_even_when_query_key_matches(monkeypatch) -> None:
    df, query_key = _make_mixed_reaction_type_key_df()
    indexed_data = {}
    reaction_type_patterns = {}
    transformation_indices = {query_key: df.copy()}

    def fake_load_db(path: str):
        return df, indexed_data, reaction_type_patterns, transformation_indices

    def fake_detect(self, smiles: str):
        return [], "Unknown"

    def fake_featurize_reaction(smiles: str, options=None):
        return {
            "reaction_type": {"reaction_type": "Suzuki_miyaura", "confidence": 0.95},
            "reaction_key": query_key,
            "aggregates": {
                "reacted_motifs": ["Ar-B(OH)2", "Ar-Cl"],
                "formed_motifs": ["Ar-Ar"],
                "spectator_motifs": [],
            },
        }

    def fake_precedent(self, *args, **kwargs):
        return []

    monkeypatch.setattr(hte, "_load_hte_database_cached", fake_load_db)
    monkeypatch.setattr(HTERecommender, "_detect_reactant_types", fake_detect)
    monkeypatch.setattr(hte, "featurize_reaction", fake_featurize_reaction)
    monkeypatch.setattr(HTERecommender, "_build_precedent_recommendations", fake_precedent)

    recommender = HTERecommender(hte_db_path="data/HTE_db")
    result = recommender.recommend(
        reactant_a_smiles="Clc1ccccc1",
        reactant_b_smiles="B(O)Oc1ccccc1",
        product_smiles="c1ccccc1-c1ccccc1",
        top_k=5,
        min_experiments=1,
        source_group="literature",
    )

    assert result.recommendations
    assert result.is_filtered_by_detected_type is True
    assert all(rec.reaction_type == "Suzuki_miyaura" for rec in result.recommendations)


def test_resolve_reaction_type_label_handles_sample_suffix() -> None:
    assert hte._resolve_reaction_type_label("suzuki_miyaura_sample500") == "Suzuki_miyaura"


def test_best_seller_scoring_beats_single_outlier(monkeypatch) -> None:
    outlier = pd.DataFrame(
        {
            "Reaction_Type_Standardized": ["Suzuki_miyaura"] * 6,
            "Reactant_A_Type": ["Ar-Cl"] * 6,
            "Reactant_B_Type": ["Ar-B(OH)2"] * 6,
            "Reactant_C_Type": [""] * 6,
            "Catalyst": ["Pd_outlier"] * 6,
            "Ligand": ["L1"] * 6,
            "Base": ["K2CO3"] * 6,
            "Solvent": ["DMF"] * 6,
            "Secondary Solvent": [""] * 6,
            "Additive": [""] * 6,
            "Coupling Reagent": [""] * 6,
            "AREA_TOTAL_REDUCED": [10.0, 8.0, 7.0, 9.0, 6.0, 5.0],
            "z-Score": [5.0, -2.0, -2.0, -2.0, -2.0, -2.0],
            "Reactant_A_Category": [""] * 6,
            "Reactant_B_Category": [""] * 6,
            "Reaction_Category": [""] * 6,
            "Is_Intramolecular": [False] * 6,
            "Source_File": ["tests"] * 6,
            "Source_Group": ["experiments"] * 6,
            "spectator_groups": [""] * 6,
            "Reactant_Signature_Core": ["Ar-B(OH)2|Ar-Cl"] * 6,
            "Reactant_Signature_Ext": ["Ar-B(OH)2|Ar-Cl"] * 6,
        }
    )
    robust = outlier.copy()
    robust["Catalyst"] = ["Pd_robust"] * 6
    robust["z-Score"] = [1.4, 1.6, 1.7, 1.8, 1.9, 2.0]
    robust["AREA_TOTAL_REDUCED"] = [70.0, 72.0, 75.0, 77.0, 80.0, 83.0]

    df = pd.concat([outlier, robust], ignore_index=True)
    key = hte._reactant_key(["Ar-B(OH)2", "Ar-Cl"])
    indexed_data = {key: df}
    reaction_type_patterns = {}
    transformation_indices = {}

    def fake_load_db(path: str):
        return df, indexed_data, reaction_type_patterns, transformation_indices

    def fake_detect(self, smiles: str):
        if "B(" in smiles:
            return ["Ar-B(OH)2"], "Aryl Boronate"
        return ["Ar-Cl"], "Aryl Halide"

    monkeypatch.setattr(hte, "_load_hte_database_cached", fake_load_db)
    monkeypatch.setattr(HTERecommender, "_detect_reactant_types", fake_detect)

    recommender = HTERecommender(hte_db_path="data/HTE_db")
    result = recommender.recommend(
        reactant_a_smiles="Clc1ccccc1",
        reactant_b_smiles="OB(O)c1ccccc1",
        top_k=2,
        min_experiments=1,
        source_group="experiments",
        reaction_type_filter="Suzuki_miyaura",
    )

    assert len(result.recommendations) == 2
    assert result.recommendations[0].catalyst == "Pd_robust"
    assert result.recommendations[0].avg_z_score > result.recommendations[1].avg_z_score


def test_aggregate_deduplicates_identical_rows_for_support_gate(monkeypatch) -> None:
    template = _make_min_hte_df()
    template["Reactant_A_Type"] = ["Ar-Cl"]
    template["Reactant_B_Type"] = ["Ar-B(OH)2"]
    template["Source_Group"] = ["experiments"]
    template["Reactant_Signature_Core"] = ["Ar-B(OH)2|Ar-Cl"]
    template["Reactant_Signature_Ext"] = ["Ar-B(OH)2|Ar-Cl"]

    duplicate_a = pd.concat([template.copy(), template.copy()], ignore_index=True)

    b1 = template.copy()
    b1["Catalyst"] = ["Pd_B"]
    b1["z-Score"] = [1.2]
    b1["AREA_TOTAL_REDUCED"] = [65.0]
    b2 = b1.copy()
    b2["z-Score"] = [1.8]
    b2["AREA_TOTAL_REDUCED"] = [75.0]

    df = pd.concat([duplicate_a, b1, b2], ignore_index=True)
    key = hte._reactant_key(["Ar-B(OH)2", "Ar-Cl"])
    indexed_data = {key: df}
    reaction_type_patterns = {}
    transformation_indices = {}

    def fake_load_db(path: str):
        return df, indexed_data, reaction_type_patterns, transformation_indices

    def fake_detect(self, smiles: str):
        if "B(" in smiles:
            return ["Ar-B(OH)2"], "Aryl Boronate"
        return ["Ar-Cl"], "Aryl Halide"

    monkeypatch.setattr(hte, "_load_hte_database_cached", fake_load_db)
    monkeypatch.setattr(HTERecommender, "_detect_reactant_types", fake_detect)

    recommender = HTERecommender(hte_db_path="data/HTE_db")
    result = recommender.recommend(
        reactant_a_smiles="Clc1ccccc1",
        reactant_b_smiles="OB(O)c1ccccc1",
        top_k=5,
        min_experiments=2,
        source_group="experiments",
        reaction_type_filter="Suzuki_miyaura",
    )

    assert result.recommendations
    catalysts = [rec.catalyst for rec in result.recommendations]
    assert "Pd_B" in catalysts
    assert "Pd" not in catalysts


def test_detected_type_prefers_family_specific_halide_fallback(monkeypatch) -> None:
    chan_lam = _make_min_hte_df()
    chan_lam["Reaction_Type_Standardized"] = ["ChanLam_dataset_converted (2)"]
    chan_lam["Reactant_A_Type"] = ["Ar-B(OH)2"]
    chan_lam["Reactant_B_Type"] = ["Ar-Cl"]
    chan_lam["Catalyst"] = ["Cu"]
    chan_lam["z-Score"] = [4.0]
    chan_lam["AREA_TOTAL_REDUCED"] = [86.0]

    suzuki = chan_lam.copy()
    suzuki["Reaction_Type_Standardized"] = ["suzuki_miyaura_sample500"]
    suzuki["Reactant_B_Type"] = ["Ar-Br"]
    suzuki["Catalyst"] = ["Pd"]
    suzuki["z-Score"] = [1.8]
    suzuki["AREA_TOTAL_REDUCED"] = [74.0]

    df = pd.concat([chan_lam, suzuki], ignore_index=True)
    indexed_data = {
        hte._reactant_key(["Ar-B(OH)2", "Ar-Cl"]): chan_lam.copy(),
        hte._reactant_key(["Ar-B(OH)2", "Ar-Br"]): suzuki.copy(),
    }
    reaction_type_patterns = {}
    transformation_indices = {}

    def fake_load_db(path: str):
        return df, indexed_data, reaction_type_patterns, transformation_indices

    def fake_detect(self, smiles: str):
        return [], "Unknown"

    def fake_featurize_reaction(smiles: str, options=None):
        return {
            "reaction_type": {"reaction_type": "Suzuki_miyaura", "confidence": 0.95},
            "reaction_key": "CRK-v1 |Ar-B(OH)2|Ar-Cl -> Ar-Ar | bond_formed: C(ar)-C(ar)",
            "aggregates": {
                "reacted_motifs": ["Ar-Cl", "Ar-B(OH)2"],
                "formed_motifs": ["Ar-Ar"],
                "spectator_motifs": [],
            },
        }

    def fake_precedent(self, *args, **kwargs):
        return []

    monkeypatch.setattr(hte, "_load_hte_database_cached", fake_load_db)
    monkeypatch.setattr(HTERecommender, "_detect_reactant_types", fake_detect)
    monkeypatch.setattr(hte, "featurize_reaction", fake_featurize_reaction)
    monkeypatch.setattr(HTERecommender, "_build_precedent_recommendations", fake_precedent)

    recommender = HTERecommender(hte_db_path="data/HTE_db")
    result = recommender.recommend(
        reactant_a_smiles="Clc1ccccc1",
        reactant_b_smiles="B(O)Oc1ccccc1",
        product_smiles="c1ccccc1-c1ccccc1",
        top_k=5,
        min_experiments=1,
        source_group="literature",
    )

    assert result.recommendations
    assert result.is_filtered_by_detected_type is True
    assert all(rec.reaction_type == "suzuki_miyaura_sample500" for rec in result.recommendations)


def test_rules_required_core_motif_path_ignores_non_core_query_token(monkeypatch) -> None:
    query_key = "CRK-v1 |Ar-B(OH)2|Ar-OMs|R_acidic-H -> Ar-Ar | bond_formed: C(ar)-C(ar)"
    df = _make_min_hte_df()
    df["Source_Group"] = ["rules"]
    df["Reaction_Type_Standardized"] = ["Suzuki_miyaura"]
    df["Reactant_A_Type"] = ["@sp2_electrophiles"]
    df["Reactant_B_Type"] = ["@organoboron"]
    df["Reactant_Signature_Core"] = ["Ar-B(OH)2|Ar-OMs"]
    df["Reactant_Signature_Ext"] = ["Ar-B(OH)2|Ar-OMs"]
    df["Catalyst"] = ["Pd(OAc)2"]
    df["Base"] = ["K3PO4"]
    df["Solvent"] = ["1,4-dioxane"]

    indexed_data = {}
    reaction_type_patterns = {}
    transformation_indices = {query_key: df.copy()}

    def fake_load_db(path: str):
        return df, indexed_data, reaction_type_patterns, transformation_indices

    def fake_detect(self, smiles: str):
        return [], "Unknown"

    def fake_featurize_reaction(smiles: str, options=None):
        return {
            "reaction_type": {"reaction_type": "Suzuki_miyaura", "confidence": 0.95},
            "reaction_key": query_key,
            "aggregates": {
                "reacted_motifs": ["Ar-B(OH)2", "Ar-OMs", "R_acidic-H"],
                "formed_motifs": ["Ar-Ar"],
                "spectator_motifs": [],
            },
        }

    def fake_precedent(self, *args, **kwargs):
        return []

    monkeypatch.setattr(hte, "_load_hte_database_cached", fake_load_db)
    monkeypatch.setattr(HTERecommender, "_detect_reactant_types", fake_detect)
    monkeypatch.setattr(hte, "featurize_reaction", fake_featurize_reaction)
    monkeypatch.setattr(HTERecommender, "_build_precedent_recommendations", fake_precedent)

    recommender = HTERecommender(hte_db_path="data/HTE_db")
    result = recommender.recommend(
        reactant_a_smiles="O=S(Oc1ccccc1)(C)=O",
        reactant_b_smiles="OB(c2ccccc2)O",
        product_smiles="c3(c4ccccc4)ccccc3",
        source_group="rules",
        top_k=5,
        min_experiments=1,
    )

    assert result.total_matching_experiments == 1
    assert "rules" in result.recommendations_by_source
    assert result.recommendations
    assert result.recommendations[0].reaction_type == "Suzuki_miyaura"


def test_precedent_recommendations_respect_source_group_filter(monkeypatch) -> None:
    df = _make_min_hte_df()
    indexed_data = {"Ar-X": df}
    reaction_type_patterns = {}
    transformation_indices = {}

    def fake_load_db(path: str):
        return df, indexed_data, reaction_type_patterns, transformation_indices

    def fake_knn(*args, **kwargs):
        return {
            "precedents": [
                {
                    "conditions": {
                        "catalyst": ["Pd(OAc)2"],
                        "ligand": ["SPhos"],
                        "base": ["K3PO4"],
                        "solvent": ["THF"],
                    },
                    "similarity": 0.91,
                    "yield": 82.0,
                    "dataset_reaction_id": "Suzuki_miyaura",
                    "rxn_type": "Suzuki",
                    "reaction_id": "lit:1",
                    "source_group": "literature",
                },
                {
                    "conditions": {
                        "catalyst": ["PdCl2"],
                        "ligand": ["PPh3"],
                        "base": ["K2CO3"],
                        "solvent": ["DMF"],
                    },
                    "similarity": 0.87,
                    "yield": 74.0,
                    "dataset_reaction_id": "Suzuki_miyaura",
                    "rxn_type": "Suzuki",
                    "reaction_id": "rules:1",
                    "source_group": "rules",
                },
            ]
        }

    monkeypatch.setattr(hte, "_load_hte_database_cached", fake_load_db)
    monkeypatch.setattr("chemtools.precedent.knn", fake_knn)

    recommender = HTERecommender(hte_db_path="data/HTE_db")
    recs = recommender._build_precedent_recommendations(
        reactant_a_smiles="Clc1ccccc1",
        reactant_b_smiles="B(O)Oc1ccccc1",
        product_smiles="c1ccccc1-c1ccccc1",
        reaction_type="Suzuki_miyaura",
        top_k=10,
        source_group="literature",
    )

    assert recs
    assert len(recs) == 1
    assert recs[0].reaction_id == "lit:1"


def test_recommend_returns_precedent_when_structured_match_missing(monkeypatch) -> None:
    df = _make_min_hte_df()
    indexed_data = {}
    reaction_type_patterns = {}
    transformation_indices = {}

    def fake_load_db(path: str):
        return df, indexed_data, reaction_type_patterns, transformation_indices

    def fake_detect(self, smiles: str):
        return ["Ar-X"], "Aryl Halide"

    def fake_precedent(self, *args, **kwargs):
        return [
            hte.ConditionRecommendation(
                catalyst="Pd",
                ligand="PPh3",
                base="K2CO3",
                solvent="DMF",
                match_score=0.95,
                reaction_type="Suzuki_miyaura",
                reaction_id="precedent:1",
            )
        ]

    monkeypatch.setattr(hte, "_load_hte_database_cached", fake_load_db)
    monkeypatch.setattr(HTERecommender, "_detect_reactant_types", fake_detect)
    monkeypatch.setattr(HTERecommender, "_build_precedent_recommendations", fake_precedent)

    recommender = HTERecommender(hte_db_path="data/HTE_db")
    result = recommender.recommend(
        reactant_a_smiles="Brc1ccccc1",
        reactant_b_smiles="OB(O)c1ccccc1",
        top_k=5,
        min_experiments=1,
        source_group="literature",
    )

    assert result.total_matching_experiments == 0
    assert "precedent" in result.recommendations_by_source
    assert len(result.recommendations_by_source["precedent"]) == 1
    assert result.recommendations
    assert result.recommendations[0].reaction_id == "precedent:1"


def test_recommend_returns_precedent_when_structured_match_missing_for_protocols(monkeypatch) -> None:
    df = _make_min_hte_df()
    indexed_data = {}
    reaction_type_patterns = {}
    transformation_indices = {}

    def fake_load_db(path: str):
        return df, indexed_data, reaction_type_patterns, transformation_indices

    def fake_detect(self, smiles: str):
        return ["Ar-X"], "Aryl Halide"

    def fake_precedent(self, *args, **kwargs):
        return [
            hte.ConditionRecommendation(
                catalyst="Pd(OAc)2",
                ligand="SPhos",
                base="K3PO4",
                solvent="THF",
                match_score=1.0,
                reaction_type="Suzuki_miyaura",
                reaction_id="protocol_exact:1",
            )
        ]

    monkeypatch.setattr(hte, "_load_hte_database_cached", fake_load_db)
    monkeypatch.setattr(HTERecommender, "_detect_reactant_types", fake_detect)
    monkeypatch.setattr(HTERecommender, "_build_precedent_recommendations", fake_precedent)

    recommender = HTERecommender(hte_db_path="data/HTE_db")
    result = recommender.recommend(
        reactant_a_smiles="Brc1ccccc1",
        reactant_b_smiles="OB(O)c1ccccc1",
        top_k=5,
        min_experiments=1,
        source_group="protocols",
    )

    assert result.total_matching_experiments == 0
    assert "precedent" in result.recommendations_by_source
    assert len(result.recommendations_by_source["precedent"]) == 1
    assert result.recommendations
    assert result.recommendations[0].reaction_id == "protocol_exact:1"


def test_normalize_hte_dataframe_prefers_detected_type_and_backfills_reaction_key(monkeypatch) -> None:
    df = pd.DataFrame(
        {
            "reaction_id": ["palladium_catalyzed_buchwald_hartwig_amination_and"],
            "detected_reaction_type": ["Suzuki_miyaura"],
            "reaction_smiles": [
                "CC(C)(C)c1ccc(OS(=O)(=O)C)cc1.B(O)(O)c2ccccc2>>CC(C)(C)c1ccc(-c2ccccc2)cc1"
            ],
            "Reaction_Key": [""],
            "reactant_1": ["Ar-OMs|R_acidic-H"],
            "reactant_2": ["Ar-B(OH)2"],
            "reactant_3": [""],
            "catalyst": ["Pd(OAc)2"],
            "base": ["Et3N/K3PO4"],
            "solvent": ["DCM/t-BuOH"],
            "yield": [0.0],
            "z_score": [0.0],
        }
    )

    def fake_featurize_reaction(smiles: str, options=None):
        return {
            "reaction_key": (
                "CRK-v1 |Ar-B(OH)2|Ar-OMs|R_acidic-H -> Ar-Ar "
                "| bond_formed: C(ar)-C(ar) "
                "| bond_broken: C(ar)-O "
                "| spectators: Ar-Alkyl "
                "| events: LGDisp+C-C"
            )
        }

    monkeypatch.setattr(hte, "featurize_reaction", fake_featurize_reaction)

    normalized = hte._normalize_hte_dataframe(df)

    assert normalized.loc[0, "Reaction_Type_Standardized"] == "Suzuki_miyaura"
    assert str(normalized.loc[0, "Reaction_Key"]).startswith("|")
    assert "bond_formed:" not in str(normalized.loc[0, "Reaction_Key"])
    assert "events:" not in str(normalized.loc[0, "Reaction_Key"])
    events_text = str(normalized.loc[0, "Reaction_Events"])
    assert "form:" in events_text
    assert "sig:" in events_text


def test_precedent_exact_reaction_rescue_works_across_family_filter(monkeypatch) -> None:
    df = _make_min_hte_df()
    indexed_data = {"Ar-X": df}
    reaction_type_patterns = {}
    transformation_indices = {}

    def fake_load_db(path: str):
        return df, indexed_data, reaction_type_patterns, transformation_indices

    def fake_knn(*args, **kwargs):
        return {"precedents": []}

    query_a = "CO.OB(O)c1ccccc1"
    query_b = "Brc1ccccn1"
    product = "c1ccncc1.c1ccc(-c2ccccn2)cc1.COc1ccccn1"

    def fake_iter_files():
        return ["fake/C_O_Coupling_canonical.csv"]

    def fake_source_group(path: str):
        return "literature"

    def fake_file_family(path: str):
        return "C_O_Coupling"

    def fake_read_records(path: str):
        return [
            {
                "reaction_smiles": f"{query_b}.{query_a}>>{product}",
                "reaction_id": "C_O_Coupling",
                "catalyst": "Pd(OAc)2",
                "base": "KOH",
                "solvent": "MeOH/THF",
                "yield": "3",
            }
        ]

    def fake_make_row(rec, *, row_index, file_family=None, source_group=None):
        return {
            "reaction_id": f"{file_family}:{row_index}",
            "dataset_reaction_id": rec.get("reaction_id"),
            "reaction_smiles": rec.get("reaction_smiles"),
            "yield_value": float(rec.get("yield") or 0.0),
            "base_uid": rec.get("base"),
            "solvent_uid": rec.get("solvent"),
            "rxn_type": "C_O_Coupling",
            "source_file": file_family,
            "source_group": source_group,
            "reference": "",
            "conditions": {
                "catalyst": rec.get("catalyst"),
                "base": rec.get("base"),
                "solvent": rec.get("solvent"),
            },
        }

    monkeypatch.setattr(hte, "_load_hte_database_cached", fake_load_db)
    monkeypatch.setattr("chemtools.precedent.knn", fake_knn)
    monkeypatch.setattr("chemtools.precedent.loader._iter_literature_files", fake_iter_files)
    monkeypatch.setattr("chemtools.precedent.loader._infer_source_group_from_path", fake_source_group)
    monkeypatch.setattr("chemtools.precedent.loader._file_family_from_name", fake_file_family)
    monkeypatch.setattr("chemtools.precedent.loader._read_csv_records", fake_read_records)
    monkeypatch.setattr("chemtools.precedent.loader._make_row_from_csv", fake_make_row)

    recommender = HTERecommender(hte_db_path="data/HTE_db")
    recs = recommender._build_precedent_recommendations(
        reactant_a_smiles=query_a,
        reactant_b_smiles=query_b,
        product_smiles=product,
        reaction_type="Suzuki_miyaura",
        top_k=5,
        source_group="literature",
    )

    assert recs
    assert recs[0].reaction_id == "C_O_Coupling:0"
    assert recs[0].reaction_type == "C_O_Coupling"
    assert recs[0].match_score == 1.0
