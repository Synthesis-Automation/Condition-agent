"""Default-setting tuner for shortlist and blend weights using LOGO predictions."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Dict, Sequence

import numpy as np
import pandas as pd

from ml.metrics import evaluate_predictions


@dataclass(frozen=True)
class TuneResult:
    shortlist_k: int
    w_stage1: float
    w_stage2: float
    apyr_mean: float
    spearman_mean: float
    top1_hit_ge_70: float
    top3_hit_ge_70: float
    objective_score: float = 0.0


def parse_int_list(text: str) -> list[int]:
    return [int(v.strip()) for v in text.split(",") if v.strip()]


def parse_float_list(text: str) -> list[float]:
    return [float(v.strip()) for v in text.split(",") if v.strip()]


def _topk_hit_rate(df: pd.DataFrame, *, score_col: str, k: int, threshold: float) -> float:
    vals = []
    for _, grp in df.groupby("reaction_smiles"):
        top = grp.sort_values(score_col, ascending=False).head(int(k))
        vals.append(float((top["condition_yield"] >= float(threshold)).any()))
    return float(np.mean(vals)) if vals else 0.0


def _score_with_shortlist(
    fold_df: pd.DataFrame,
    *,
    shortlist_k: int,
    w_stage1: float,
    w_stage2: float,
) -> np.ndarray:
    score = w_stage1 * fold_df["stage1_score"].to_numpy(dtype=float) + w_stage2 * fold_df["stage2_score"].to_numpy(
        dtype=float
    )
    if shortlist_k <= 0:
        return score
    rank = fold_df.groupby("reaction_smiles")["stage1_score"].rank(method="first", ascending=False).to_numpy(dtype=float)
    out = score.copy()
    out[rank > float(shortlist_k)] = -1.0e9
    return out


def evaluate_setting(
    pred_df: pd.DataFrame,
    *,
    shortlist_k: int,
    w_stage1: float,
    w_stage2: float,
) -> TuneResult:
    apyr_vals = []
    rho_vals = []
    top1_vals = []
    top3_vals = []
    for _, fold in pred_df.groupby("fold"):
        work = fold.copy()
        work["score"] = _score_with_shortlist(work, shortlist_k=shortlist_k, w_stage1=w_stage1, w_stage2=w_stage2)
        m = evaluate_predictions(
            work["condition_yield"].to_numpy(dtype=float),
            work["score"].to_numpy(dtype=float),
            work["group"].astype(str).to_numpy(),
        )
        apyr_vals.append(float(m.apyr_mean))
        rho_vals.append(float(m.spearman))
        top1_vals.append(_topk_hit_rate(work, score_col="score", k=1, threshold=70.0))
        top3_vals.append(_topk_hit_rate(work, score_col="score", k=3, threshold=70.0))
    return TuneResult(
        shortlist_k=int(shortlist_k),
        w_stage1=float(w_stage1),
        w_stage2=float(w_stage2),
        apyr_mean=float(np.mean(apyr_vals)),
        spearman_mean=float(np.mean(rho_vals)),
        top1_hit_ge_70=float(np.mean(top1_vals)),
        top3_hit_ge_70=float(np.mean(top3_vals)),
    )


def grid_search_defaults(
    pred_df: pd.DataFrame,
    *,
    shortlist_grid: Sequence[int],
    w_stage1_grid: Sequence[float],
) -> tuple[TuneResult, list[Dict[str, Any]]]:
    rows: list[Dict[str, Any]] = []
    best: TuneResult | None = None
    for k in shortlist_grid:
        for w1 in w_stage1_grid:
            w1f = float(w1)
            w2f = float(1.0 - w1f)
            if w2f < 0.0:
                continue
            result = evaluate_setting(pred_df, shortlist_k=int(k), w_stage1=w1f, w_stage2=w2f)
            row = {
                "shortlist_k": result.shortlist_k,
                "w_stage1": result.w_stage1,
                "w_stage2": result.w_stage2,
                "apyr_mean": result.apyr_mean,
                "spearman_mean": result.spearman_mean,
                "top1_hit_ge_70": result.top1_hit_ge_70,
                "top3_hit_ge_70": result.top3_hit_ge_70,
            }
            rows.append(row)
            if best is None:
                best = result
                continue
            lhs = (result.apyr_mean, result.top1_hit_ge_70, result.spearman_mean)
            rhs = (best.apyr_mean, best.top1_hit_ge_70, best.spearman_mean)
            if lhs > rhs:
                best = result
    if best is None:
        raise ValueError("No candidate settings evaluated.")
    return best, rows


def _minmax(values: np.ndarray) -> np.ndarray:
    lo = float(np.min(values))
    hi = float(np.max(values))
    if hi <= lo:
        return np.ones_like(values, dtype=float)
    return (values - lo) / (hi - lo)


def score_composite_objective(
    rows: list[Dict[str, Any]],
    *,
    w_apyr: float,
    w_top1: float,
    w_spearman: float,
    w_top3: float,
    w_stage2_pref: float,
) -> list[Dict[str, Any]]:
    if not rows:
        return rows
    df = pd.DataFrame(rows).copy()
    apyr_n = _minmax(df["apyr_mean"].to_numpy(dtype=float))
    top1_n = _minmax(df["top1_hit_ge_70"].to_numpy(dtype=float))
    spear_n = _minmax(df["spearman_mean"].to_numpy(dtype=float))
    top3_n = _minmax(df["top3_hit_ge_70"].to_numpy(dtype=float))
    stage2_n = df["w_stage2"].to_numpy(dtype=float)
    total_w = float(w_apyr + w_top1 + w_spearman + w_top3 + w_stage2_pref)
    if total_w <= 0.0:
        total_w = 1.0
    score = (
        float(w_apyr) * apyr_n
        + float(w_top1) * top1_n
        + float(w_spearman) * spear_n
        + float(w_top3) * top3_n
        + float(w_stage2_pref) * stage2_n
    ) / total_w
    df["objective_score"] = score
    return df.to_dict(orient="records")


def select_best_by_composite(rows: list[Dict[str, Any]]) -> TuneResult:
    if not rows:
        raise ValueError("No rows available for composite selection.")
    best_row = sorted(
        rows,
        key=lambda r: (
            float(r.get("objective_score", 0.0)),
            float(r.get("apyr_mean", 0.0)),
            float(r.get("top1_hit_ge_70", 0.0)),
            float(r.get("spearman_mean", 0.0)),
        ),
        reverse=True,
    )[0]
    return TuneResult(
        shortlist_k=int(best_row["shortlist_k"]),
        w_stage1=float(best_row["w_stage1"]),
        w_stage2=float(best_row["w_stage2"]),
        apyr_mean=float(best_row["apyr_mean"]),
        spearman_mean=float(best_row["spearman_mean"]),
        top1_hit_ge_70=float(best_row["top1_hit_ge_70"]),
        top3_hit_ge_70=float(best_row["top3_hit_ge_70"]),
        objective_score=float(best_row.get("objective_score", 0.0)),
    )


def filter_rows_by_min_stage2(rows: list[Dict[str, Any]], *, min_stage2_weight: float) -> list[Dict[str, Any]]:
    threshold = float(min_stage2_weight)
    if threshold <= 0.0:
        return list(rows)
    feasible = [r for r in rows if float(r.get("w_stage2", 0.0)) >= threshold]
    return feasible if feasible else list(rows)


def tune_from_prediction_csv(
    prediction_csv: Path,
    *,
    output_json: Path,
    shortlist_grid: Sequence[int],
    w_stage1_grid: Sequence[float],
    objective_mode: str,
    objective_weights: Dict[str, float],
    min_stage2_weight: float,
) -> Dict[str, Any]:
    pred_df = pd.read_csv(prediction_csv)
    best_legacy, rows = grid_search_defaults(pred_df, shortlist_grid=shortlist_grid, w_stage1_grid=w_stage1_grid)
    if objective_mode == "composite":
        rows = score_composite_objective(
            rows,
            w_apyr=float(objective_weights.get("apyr", 0.5)),
            w_top1=float(objective_weights.get("top1", 0.25)),
            w_spearman=float(objective_weights.get("spearman", 0.15)),
            w_top3=float(objective_weights.get("top3", 0.05)),
            w_stage2_pref=float(objective_weights.get("stage2_pref", 0.05)),
        )
        feasible_rows = filter_rows_by_min_stage2(rows, min_stage2_weight=min_stage2_weight)
        best = select_best_by_composite(feasible_rows)
    elif objective_mode == "apyr_lexicographic":
        best = best_legacy
        feasible_rows = rows
    else:
        raise ValueError(f"Unsupported objective_mode: {objective_mode}")
    report = {
        "prediction_csv": str(prediction_csv),
        "n_rows": int(len(pred_df)),
        "n_folds": int(pred_df["fold"].nunique()) if "fold" in pred_df.columns else 0,
        "objective_mode": objective_mode,
        "objective_weights": objective_weights,
        "min_stage2_weight": float(min_stage2_weight),
        "search_space": {
            "shortlist_grid": [int(v) for v in shortlist_grid],
            "w_stage1_grid": [float(v) for v in w_stage1_grid],
            "w_stage2_grid_implied": "1 - w_stage1",
        },
        "best": {
            "shortlist_k": best.shortlist_k,
            "w_stage1": best.w_stage1,
            "w_stage2": best.w_stage2,
            "apyr_mean": best.apyr_mean,
            "spearman_mean": best.spearman_mean,
            "top1_hit_ge_70": best.top1_hit_ge_70,
            "top3_hit_ge_70": best.top3_hit_ge_70,
            "objective_score": best.objective_score,
        },
        "all_results": rows,
        "feasible_results_count": int(len(feasible_rows)),
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    pd.DataFrame(rows).to_csv(output_json.with_suffix(".grid.csv"), index=False)
    return report
