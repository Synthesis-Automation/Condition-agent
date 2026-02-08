"""Evaluation metrics for condition recommendation models."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


@dataclass
class EvaluationResult:
    spearman: float
    apyr_mean: float
    apyr_std: float
    n_groups: int


def percentile_rank(values: np.ndarray, score: float) -> float:
    if values.size == 0:
        return 0.0
    return float(100.0 * np.mean(values <= score))


def compute_apyr(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    groups: Sequence[str],
    *,
    top_within: float = 5.0,
) -> Tuple[float, float, int]:
    if len(y_true) == 0:
        return 0.0, 0.0, 0
    tmp = pd.DataFrame({"y": y_true, "pred": y_pred, "group": groups})
    apyr_values: List[float] = []
    for _, grp in tmp.groupby("group"):
        exp = grp["y"].to_numpy(dtype=float)
        pred = grp["pred"].to_numpy(dtype=float)
        if exp.size == 0:
            continue
        threshold = float(np.max(pred) - top_within)
        chosen_exp = exp[pred >= threshold]
        if chosen_exp.size == 0:
            continue
        ranks = [percentile_rank(exp, val) for val in chosen_exp]
        apyr_values.append(float(np.mean(ranks)))
    if not apyr_values:
        return 0.0, 0.0, 0
    return float(np.mean(apyr_values)), float(np.std(apyr_values)), len(apyr_values)


def evaluate_predictions(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    groups: Sequence[str],
) -> EvaluationResult:
    rho = spearmanr(y_true, y_pred).correlation
    spearman = 0.0 if rho is None or math.isnan(rho) else float(rho)
    apyr_mean, apyr_std, n_groups = compute_apyr(y_true, y_pred, groups)
    return EvaluationResult(
        spearman=spearman,
        apyr_mean=apyr_mean,
        apyr_std=apyr_std,
        n_groups=n_groups,
    )
