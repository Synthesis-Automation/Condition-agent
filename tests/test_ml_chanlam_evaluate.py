import numpy as np

from ml.evaluate_chanlam import _topk_metrics_for_fold, compute_apyr


def test_topk_metrics_fold() -> None:
    y_true = np.array([20.0, 80.0, 60.0, 10.0], dtype=float)
    y_pred = np.array([0.1, 0.9, 0.7, 0.2], dtype=float)
    out = _topk_metrics_for_fold(
        y_true,
        y_pred,
        top_ks=[1, 3],
        yield_thresholds=[50.0, 70.0],
    )
    assert out["top1_hit_ge_50"] == 1.0
    assert out["top1_hit_ge_70"] == 1.0
    assert out["top3_hit_ge_50"] == 1.0
    assert out["top3_hit_ge_70"] == 1.0


def test_compute_apyr_nonzero() -> None:
    y_true = np.array([10.0, 90.0, 20.0, 80.0], dtype=float)
    y_pred = np.array([15.0, 93.0, 19.0, 81.0], dtype=float)
    groups = np.array(["A", "A", "B", "B"])
    mean_apyr, std_apyr, n_groups = compute_apyr(y_true, y_pred, groups, top_within=5.0)
    assert n_groups == 2
    assert mean_apyr > 70.0
    assert std_apyr >= 0.0

