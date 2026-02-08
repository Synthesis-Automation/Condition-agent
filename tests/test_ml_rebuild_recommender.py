import pandas as pd

from ml.contracts import FeatureConfig, RecommenderConfig, Stage2Config
from ml.recommender import TwoStageRecommender


def _tiny_training_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "yield": [80.0, 60.0, 30.0, 20.0, 90.0, 70.0],
            "sulfonamide_smiles": ["S1", "S1", "S2", "S2", "S3", "S3"],
            "boronic_smiles": ["B1", "B1", "B2", "B2", "B3", "B3"],
            "reaction_smiles": ["R1", "R1", "R2", "R2", "R3", "R3"],
            "catalyst": ["Cu(OAc)2", "Cu(OAc)2", "Cu(OTf)2", "Cu(OTf)2", "Cu(OAc)2", "Cu(OTf)2"],
            "base": ["K2CO3", "CsF", "K2CO3", "CsF", "K2CO3", "CsF"],
            "solvent": ["DCE", "DCE", "MeCN", "MeCN", "DCE", "MeCN"],
            "formed_motifs_tokens": ["Ar-N", "Ar-N", "Ar-N", "Ar-N", "Ar-N", "Ar-N"],
            "spectator_groups_tokens": ["F OR", "F OR", "Cl", "Cl", "F", "F"],
            "sulf_motif_tokens": ["Ar-SO2NH2"] * 6,
            "bor_motif_tokens": ["Ar-B(OH)2"] * 6,
            "sulf_motif_count": [3, 3, 2, 2, 4, 4],
            "sulf_aryl_steric_max": [1.0, 1.0, 0.8, 0.8, 1.2, 1.2],
            "sulf_alkyl_steric_max": [0.2, 0.2, 0.1, 0.1, 0.3, 0.3],
            "sulf_aryl_electronic_avg": [6.0, 6.0, 5.0, 5.0, 6.5, 6.5],
            "bor_motif_count": [2, 2, 2, 2, 3, 3],
            "bor_aryl_steric_max": [0.4, 0.4, 0.3, 0.3, 0.5, 0.5],
            "bor_alkyl_steric_max": [0.0, 0.0, 0.0, 0.0, 0.1, 0.1],
            "bor_aryl_electronic_avg": [5.0, 5.0, 5.0, 5.0, 5.2, 5.2],
        }
    )


def test_two_stage_recommender_fit_and_score_rows() -> None:
    cfg = RecommenderConfig(
        feature=FeatureConfig(profile="core_full", with_condition_props=False),
        stage2=Stage2Config(model_type="rf_reg", n_estimators=40, random_state=7),
    )
    rec = TwoStageRecommender.create(cfg)
    train_df = _tiny_training_df()
    rec.fit(train_df)
    scored = rec.score_rows(train_df.head(3))
    assert "stage1_score" in scored.columns
    assert "stage2_score" in scored.columns
    assert "final_score" in scored.columns
    assert len(scored) == 3
