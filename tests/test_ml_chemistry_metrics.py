import numpy as np
import pandas as pd

from ml.chemistry import extract_reaction_substrates, tokenize_text_field
from ml.metrics import compute_apyr


def test_extract_reaction_substrates_chanlam() -> None:
    rsmi = "NS(=O)(=O)c1cccc(F)n1.COc1ccc(B(O)O)cc1>>COc1ccc(NS(=O)(=O)c2cccc(F)n2)cc1"
    sulf, bor = extract_reaction_substrates(rsmi)
    assert sulf == "NS(=O)(=O)c1cccc(F)n1"
    assert bor == "COc1ccc(B(O)O)cc1"


def test_tokenize_text_field_splits_and_normalizes() -> None:
    text = "Ar-N / HeteroAr, OR; Ar-N"
    tok = tokenize_text_field(text)
    assert tok == "Ar-N HeteroAr OR"


def test_compute_apyr_prefers_top_predictions() -> None:
    y_true = np.array([10.0, 90.0, 30.0, 80.0], dtype=float)
    y_pred = np.array([15.0, 95.0, 20.0, 82.0], dtype=float)
    groups = np.array(["A", "A", "B", "B"])
    mean_apyr, std_apyr, n_groups = compute_apyr(y_true, y_pred, groups, top_within=5.0)
    assert n_groups == 2
    assert mean_apyr > 70.0
    assert std_apyr >= 0.0


def test_compute_apyr_empty_returns_zero() -> None:
    y_true = np.array([], dtype=float)
    y_pred = np.array([], dtype=float)
    groups = pd.Series([], dtype=str)
    m, s, n = compute_apyr(y_true, y_pred, groups)
    assert (m, s, n) == (0.0, 0.0, 0)
