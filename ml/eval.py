"""Evaluation harness for rebuilt two-stage recommender."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Sequence

import numpy as np
import pandas as pd

from ml.contracts import EvalConfig, RecommenderConfig
from ml.data import assert_no_reaction_leakage, load_chanlam_feature_table, logo_splits
from ml.metrics import evaluate_predictions
from ml.recommender import TwoStageRecommender


def _topk_hits(
    df: pd.DataFrame,
    *,
    score_col: str,
    reaction_col: str,
    yield_col: str,
    top_ks: Sequence[int],
    thresholds: Sequence[float],
) -> Dict[str, float]:
    out: Dict[str, float] = {}
    grouped = list(df.groupby(reaction_col))
    for k in top_ks:
        for t in thresholds:
            vals = []
            for _, grp in grouped:
                top = grp.sort_values(score_col, ascending=False).head(int(k))
                vals.append(float((top[yield_col] >= float(t)).any()))
            out[f"top{k}_hit_ge_{int(t)}"] = float(np.mean(vals)) if vals else 0.0
    return out


def _apply_shortlist(
    df: pd.DataFrame,
    *,
    stage1_col: str,
    score_col: str,
    reaction_col: str,
    shortlist_k: int,
) -> np.ndarray:
    if shortlist_k <= 0:
        return df[score_col].to_numpy(dtype=float)
    work = df[[reaction_col, stage1_col, score_col]].copy()
    work["_rank"] = work.groupby(reaction_col)[stage1_col].rank(method="first", ascending=False)
    out = work[score_col].to_numpy(dtype=float)
    mask = work["_rank"].to_numpy(dtype=float) > float(shortlist_k)
    out[mask] = -1.0e9
    return out


def run_logo_backtest(
    rec_cfg: RecommenderConfig,
    eval_cfg: EvalConfig,
) -> Dict[str, Any]:
    df = load_chanlam_feature_table(rec_cfg.dataset, with_condition_props=rec_cfg.feature.with_condition_props)
    fold_rows: list[dict[str, Any]] = []
    pred_rows: list[pd.DataFrame] = []

    max_folds = int(eval_cfg.max_folds)
    for fold_idx, train_idx, test_idx in logo_splits(df, group_col=rec_cfg.dataset.group_col):
        if max_folds > 0 and fold_idx > max_folds:
            break
        train_df = df.iloc[train_idx].copy()
        test_df = df.iloc[test_idx].copy()
        assert_no_reaction_leakage(train_df, test_df, reaction_col=rec_cfg.dataset.reaction_col)

        rec = TwoStageRecommender.create(rec_cfg).fit(train_df)
        scored = rec.score_rows(test_df)
        final_short = _apply_shortlist(
            scored,
            stage1_col="stage1_score",
            score_col="final_score",
            reaction_col=rec_cfg.dataset.reaction_col,
            shortlist_k=rec_cfg.stage1.shortlist_k,
        )
        scored["final_score_shortlist"] = final_short

        y_true = pd.to_numeric(scored[rec_cfg.dataset.yield_col], errors="coerce").fillna(0.0).to_numpy(dtype=float)
        groups = scored[rec_cfg.dataset.group_col].fillna("NA").astype(str).to_numpy()
        fold = {
            "fold": int(fold_idx),
            "group": str(groups[0]) if len(groups) else "NA",
            "rows": int(len(scored)),
        }
        for name, col in (
            ("stage1", "stage1_score"),
            ("stage2", "stage2_score"),
            ("final", "final_score_shortlist"),
        ):
            m = evaluate_predictions(y_true, scored[col].to_numpy(dtype=float), groups)
            fold[f"{name}_spearman"] = float(m.spearman)
            fold[f"{name}_apyr"] = float(m.apyr_mean)
            hits = _topk_hits(
                scored,
                score_col=col,
                reaction_col=rec_cfg.dataset.reaction_col,
                yield_col=rec_cfg.dataset.yield_col,
                top_ks=eval_cfg.top_ks,
                thresholds=eval_cfg.thresholds,
            )
            for k, v in hits.items():
                fold[f"{name}_{k}"] = float(v)
        fold_rows.append(fold)

        pred_rows.append(
            pd.DataFrame(
                {
                    "fold": int(fold_idx),
                    "group": groups,
                    "reaction_smiles": scored[rec_cfg.dataset.reaction_col].astype(str).to_numpy(),
                    "condition_yield": y_true,
                    "stage1_score": scored["stage1_score"].to_numpy(dtype=float),
                    "stage2_score": scored["stage2_score"].to_numpy(dtype=float),
                    "final_score": scored["final_score"].to_numpy(dtype=float),
                    "final_score_shortlist": scored["final_score_shortlist"].to_numpy(dtype=float),
                }
            )
        )

    fold_df = pd.DataFrame(fold_rows)
    pred_df = pd.concat(pred_rows, axis=0, ignore_index=True) if pred_rows else pd.DataFrame()

    def _summ(prefix: str) -> Dict[str, float]:
        out = {
            "spearman_mean": float(fold_df[f"{prefix}_spearman"].mean()),
            "spearman_std": float(fold_df[f"{prefix}_spearman"].std(ddof=0)),
            "apyr_mean": float(fold_df[f"{prefix}_apyr"].mean()),
            "apyr_std": float(fold_df[f"{prefix}_apyr"].std(ddof=0)),
        }
        for k in eval_cfg.top_ks:
            for t in eval_cfg.thresholds:
                col = f"{prefix}_top{k}_hit_ge_{int(t)}"
                out[f"top{k}_hit_ge_{int(t)}"] = float(fold_df[col].mean())
        return out

    report = {
        "input_csv": str(rec_cfg.dataset.input_csv),
        "rows": int(len(df)),
        "folds": int(len(fold_df)),
        "config": {
            "feature_profile": rec_cfg.feature.profile,
            "with_condition_props": bool(rec_cfg.feature.with_condition_props),
            "stage1_shortlist_k": int(rec_cfg.stage1.shortlist_k),
            "stage1_shrinkage_m": float(rec_cfg.stage1.shrinkage_m),
            "stage2_model_type": rec_cfg.stage2.model_type,
            "stage2_n_estimators": int(rec_cfg.stage2.n_estimators),
            "blend_weights": {
                "w_stage1": float(rec_cfg.blend.normalized()[0]),
                "w_stage2": float(rec_cfg.blend.normalized()[1]),
            },
            "max_folds": int(eval_cfg.max_folds),
        },
        "metrics": {
            "stage1": _summ("stage1"),
            "stage2": _summ("stage2"),
            "final": _summ("final"),
        },
        "comparisons": {
            "final_minus_stage2_apyr": float(_summ("final")["apyr_mean"] - _summ("stage2")["apyr_mean"]),
            "final_minus_stage2_spearman": float(_summ("final")["spearman_mean"] - _summ("stage2")["spearman_mean"]),
        },
    }

    out_json = Path(eval_cfg.output_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    fold_df.to_csv(out_json.with_suffix(".per_fold.csv"), index=False)
    if not pred_df.empty:
        pred_df.to_csv(out_json.with_suffix(".predictions.csv"), index=False)
    return report
