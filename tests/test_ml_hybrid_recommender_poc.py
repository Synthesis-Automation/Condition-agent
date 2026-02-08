import numpy as np
import pandas as pd

from ml.hybrid_recommender_poc import _condition_key, _topk_hit_metrics


def test_condition_key_uses_catalyst_base_solvent() -> None:
    df = pd.DataFrame(
        {
            "catalyst": ["Cu(OAc)2", None],
            "base": ["K2CO3", "CsF"],
            "solvent": ["DCE", None],
        }
    )
    out = _condition_key(df).tolist()
    assert out[0] == "Cu(OAc)2|K2CO3|DCE"
    assert out[1] == "NA|CsF|NA"


def test_topk_hit_metrics_basic() -> None:
    y = np.array([20.0, 80.0, 30.0, 90.0], dtype=float)
    s = np.array([0.1, 0.9, 0.2, 0.8], dtype=float)
    rxn = np.array(["R1", "R1", "R2", "R2"])
    m = _topk_hit_metrics(y, s, rxn, top_ks=[1], thresholds=[70.0])
    assert "top1_hit_ge_70" in m
    assert abs(m["top1_hit_ge_70"] - 1.0) < 1e-9
