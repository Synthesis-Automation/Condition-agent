import numpy as np

from ml.chanlam_pipeline import _extract_reaction_substrates, compute_apyr


def test_extract_reaction_substrates_chanlam() -> None:
    rsmi = "NS(=O)(=O)c1cccc(F)n1.COc1ccc(B(O)O)cc1>>COc1ccc(NS(=O)(=O)c2cccc(F)n2)cc1"
    sulfonamide, boronic = _extract_reaction_substrates(rsmi)
    assert sulfonamide == "NS(=O)(=O)c1cccc(F)n1"
    assert boronic == "COc1ccc(B(O)O)cc1"


def test_compute_apyr_prefers_top_predictions() -> None:
    y_true = np.array([10.0, 90.0, 30.0, 80.0], dtype=float)
    y_pred = np.array([15.0, 95.0, 20.0, 82.0], dtype=float)
    groups = np.array(["A", "A", "B", "B"])
    mean_apyr, std_apyr, n_groups = compute_apyr(y_true, y_pred, groups, top_within=5.0)
    assert n_groups == 2
    assert mean_apyr > 70.0
    assert std_apyr >= 0.0

