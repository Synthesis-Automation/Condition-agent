import pandas as pd

from ml.contracts import BlendConfig, FeatureConfig, RecommenderConfig, Stage2Config
from ml.features import resolve_feature_spec
from ml.models import ConditionPriorRanker


def test_blend_normalized_handles_positive_weights() -> None:
    w1, w2 = BlendConfig(w_stage1=2.0, w_stage2=6.0).normalized()
    assert abs(w1 - 0.25) < 1e-9
    assert abs(w2 - 0.75) < 1e-9


def test_resolve_feature_spec_base_profile() -> None:
    spec = resolve_feature_spec(profile="base_motif_spectator", with_condition_props=False)
    assert spec.categorical_cols == ["base"]
    assert "spectator_groups_tokens" in spec.text_cols


def test_condition_prior_ranker_scores_seen_and_unseen() -> None:
    df = pd.DataFrame(
        {
            "catalyst": ["Cu(OAc)2", "Cu(OAc)2", "Cu(OTf)2"],
            "base": ["K2CO3", "K2CO3", "CsF"],
            "solvent": ["DCE", "DCE", "MeCN"],
            "yield": [90.0, 70.0, 20.0],
        }
    )
    prior = ConditionPriorRanker(condition_cols=["catalyst", "base", "solvent"], shrinkage_m=10.0).fit(df)
    seen = prior.score(df.iloc[[0]])[0]
    unseen = prior.score(
        pd.DataFrame([{"catalyst": "X", "base": "Y", "solvent": "Z"}])
    )[0]
    assert seen > unseen


def test_recommender_config_constructs() -> None:
    cfg = RecommenderConfig(feature=FeatureConfig(profile="core_full"), stage2=Stage2Config(model_type="rf_reg"))
    assert cfg.feature.profile == "core_full"
    assert cfg.stage2.model_type == "rf_reg"
