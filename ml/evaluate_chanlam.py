"""Evaluation harness for Chan-Lam condition-selection models."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import argparse
import json
import math
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.model_selection import LeaveOneGroupOut, ShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from ml.chanlam_pipeline import (
    _CONDITION_COLUMNS,
    _NUMERIC_COLUMNS,
    _TEXT_COLUMNS,
    build_feature_table,
    load_chanlam_dataset,
)


@dataclass
class EvalResult:
    spearman: float
    apyr_mean: float
    apyr_std: float
    n_groups: int


def _percentile_rank(values: np.ndarray, score: float) -> float:
    if values.size == 0:
        return 0.0
    return float(100.0 * np.mean(values <= score))


def compute_apyr(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    groups: Sequence[str],
    *,
    top_within: float = 5.0,
) -> tuple[float, float, int]:
    tmp = pd.DataFrame({"y": y_true, "pred": y_pred, "group": groups})
    apyr_values: List[float] = []
    for _, grp in tmp.groupby("group"):
        exp = grp["y"].to_numpy(dtype=float)
        pred = grp["pred"].to_numpy(dtype=float)
        if exp.size == 0:
            continue
        chosen = exp[pred >= np.max(pred) - top_within]
        if chosen.size == 0:
            continue
        ranks = [_percentile_rank(exp, val) for val in chosen]
        apyr_values.append(float(np.mean(ranks)))
    if not apyr_values:
        return 0.0, 0.0, 0
    return float(np.mean(apyr_values)), float(np.std(apyr_values)), len(apyr_values)


def evaluate_predictions(y_true: np.ndarray, y_pred: np.ndarray, groups: Sequence[str]) -> EvalResult:
    rho = spearmanr(y_true, y_pred).correlation
    spearman = 0.0 if rho is None or math.isnan(rho) else float(rho)
    apyr_mean, apyr_std, n_groups = compute_apyr(y_true, y_pred, groups)
    return EvalResult(spearman=spearman, apyr_mean=apyr_mean, apyr_std=apyr_std, n_groups=n_groups)


def _build_pipeline(
    *,
    categorical_cols: list[str],
    numeric_cols: list[str],
    text_cols: list[str],
    random_state: int,
    n_estimators: int,
) -> Pipeline:
    transformers: list[tuple[str, Any, Any]] = []
    if numeric_cols:
        transformers.append(("num", StandardScaler(), numeric_cols))
    if categorical_cols:
        transformers.append(("cat", OneHotEncoder(handle_unknown="ignore"), categorical_cols))
    for col in text_cols:
        transformers.append(
            (
                f"text_{col}",
                CountVectorizer(token_pattern=r"[^ ]+", binary=True),
                col,
            )
        )

    preprocess = ColumnTransformer(transformers=transformers, remainder="drop", sparse_threshold=0.3)
    model = RandomForestRegressor(
        n_estimators=n_estimators,
        min_samples_leaf=2,
        random_state=random_state,
        n_jobs=-1,
    )
    return Pipeline([("preprocess", preprocess), ("model", model)])


def _topk_metrics_for_fold(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    *,
    top_ks: Sequence[int],
    yield_thresholds: Sequence[float],
) -> Dict[str, float]:
    order = np.argsort(y_pred)[::-1]
    sorted_true = y_true[order]
    out: Dict[str, float] = {}
    for k in top_ks:
        chosen = sorted_true[:k]
        if chosen.size == 0:
            out[f"top{k}_max_yield"] = 0.0
            for t in yield_thresholds:
                out[f"top{k}_hit_ge_{int(t)}"] = 0.0
            continue
        out[f"top{k}_max_yield"] = float(np.max(chosen))
        for t in yield_thresholds:
            out[f"top{k}_hit_ge_{int(t)}"] = float(np.any(chosen >= t))
        out[f"top{k}_best_percentile"] = _percentile_rank(y_true, float(np.max(chosen)))
    return out


def run_logo_evaluation(
    X: pd.DataFrame,
    y: np.ndarray,
    groups: np.ndarray,
    *,
    pipeline_builder: Any,
    top_ks: Sequence[int],
    yield_thresholds: Sequence[float],
) -> Dict[str, Any]:
    logo = LeaveOneGroupOut()
    per_fold: List[Dict[str, Any]] = []
    y_true_all: List[float] = []
    y_pred_all: List[float] = []
    groups_all: List[str] = []

    for fold_idx, (train_idx, test_idx) in enumerate(logo.split(X, y, groups=groups), start=1):
        pipe = pipeline_builder()
        pipe.fit(X.iloc[train_idx], y[train_idx])
        pred = pipe.predict(X.iloc[test_idx])
        grp_name = str(groups[test_idx][0]) if len(test_idx) else "NA"
        result = evaluate_predictions(y[test_idx], pred, groups[test_idx])
        fold_row: Dict[str, Any] = {
            "fold": fold_idx,
            "group": grp_name,
            "rows": int(len(test_idx)),
            "spearman": result.spearman,
            "apyr_mean": result.apyr_mean,
            "apyr_std": result.apyr_std,
        }
        fold_row.update(
            _topk_metrics_for_fold(
                y[test_idx],
                pred,
                top_ks=top_ks,
                yield_thresholds=yield_thresholds,
            )
        )
        per_fold.append(fold_row)
        y_true_all.extend(y[test_idx].tolist())
        y_pred_all.extend(pred.tolist())
        groups_all.extend(groups[test_idx].tolist())

    global_eval = evaluate_predictions(np.array(y_true_all), np.array(y_pred_all), groups_all)
    fold_df = pd.DataFrame(per_fold)
    summary = {
        "folds": int(len(per_fold)),
        "global_spearman": global_eval.spearman,
        "global_apyr_mean": global_eval.apyr_mean,
        "global_apyr_std": global_eval.apyr_std,
        "global_apyr_n_groups": int(global_eval.n_groups),
        "per_fold": per_fold,
    }
    for k in top_ks:
        summary[f"top{k}_avg_best_percentile"] = float(fold_df[f"top{k}_best_percentile"].mean())
        summary[f"top{k}_avg_max_yield"] = float(fold_df[f"top{k}_max_yield"].mean())
        for t in yield_thresholds:
            summary[f"top{k}_hit_rate_ge_{int(t)}"] = float(fold_df[f"top{k}_hit_ge_{int(t)}"].mean())
    return summary


def run_random_split_evaluation(
    X: pd.DataFrame,
    y: np.ndarray,
    groups: np.ndarray,
    *,
    pipeline_builder: Any,
    n_splits: int,
    test_size: float,
    random_state: int,
) -> Dict[str, Any]:
    splitter = ShuffleSplit(n_splits=n_splits, test_size=test_size, random_state=random_state)
    rows: List[Dict[str, Any]] = []
    for fold_idx, (train_idx, test_idx) in enumerate(splitter.split(X), start=1):
        pipe = pipeline_builder()
        pipe.fit(X.iloc[train_idx], y[train_idx])
        pred = pipe.predict(X.iloc[test_idx])
        result = evaluate_predictions(y[test_idx], pred, groups[test_idx])
        rows.append(
            {
                "fold": fold_idx,
                "spearman": result.spearman,
                "apyr_mean": result.apyr_mean,
                "apyr_std": result.apyr_std,
            }
        )
    fold_df = pd.DataFrame(rows)
    return {
        "splits": n_splits,
        "test_size": test_size,
        "spearman_mean": float(fold_df["spearman"].mean()),
        "spearman_std": float(fold_df["spearman"].std(ddof=0)),
        "apyr_mean": float(fold_df["apyr_mean"].mean()),
        "apyr_std": float(fold_df["apyr_mean"].std(ddof=0)),
        "per_fold": rows,
    }


def evaluate_models(
    input_csv: Path,
    output_dir: Path,
    *,
    n_estimators: int,
    random_state: int,
    random_splits: int,
    random_test_size: float,
    top_ks: Sequence[int],
    yield_thresholds: Sequence[float],
) -> Dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_df = load_chanlam_dataset(input_csv)
    feat_df = build_feature_table(raw_df)

    y = feat_df["yield"].to_numpy(dtype=float)
    groups = feat_df["sulfonamide_smiles"].fillna("NA").astype(str).to_numpy()

    model_specs: Dict[str, Dict[str, list[str]]] = {
        "condition_only": {
            "categorical": list(_CONDITION_COLUMNS),
            "numeric": [],
            "text": [],
        },
        "condition_plus_substrate_ohe": {
            "categorical": list(_CONDITION_COLUMNS) + ["sulfonamide_smiles", "boronic_smiles"],
            "numeric": [],
            "text": [],
        },
        "descriptor_model": {
            "categorical": list(_CONDITION_COLUMNS),
            "numeric": list(_NUMERIC_COLUMNS),
            "text": list(_TEXT_COLUMNS),
        },
    }

    report: Dict[str, Any] = {
        "input_csv": str(input_csv),
        "rows": int(len(feat_df)),
        "unique_sulfonamides": int(feat_df["sulfonamide_smiles"].nunique()),
        "unique_boronic_partners": int(feat_df["boronic_smiles"].nunique()),
        "top_ks": list(top_ks),
        "yield_thresholds": list(yield_thresholds),
        "models": {},
    }

    for model_name, spec in model_specs.items():
        cols = spec["categorical"] + spec["numeric"] + spec["text"]
        X = feat_df[cols].copy()

        def make_pipe() -> Pipeline:
            return _build_pipeline(
                categorical_cols=spec["categorical"],
                numeric_cols=spec["numeric"],
                text_cols=spec["text"],
                random_state=random_state,
                n_estimators=n_estimators,
            )

        logo_result = run_logo_evaluation(
            X,
            y,
            groups,
            pipeline_builder=make_pipe,
            top_ks=top_ks,
            yield_thresholds=yield_thresholds,
        )
        random_result = run_random_split_evaluation(
            X,
            y,
            groups,
            pipeline_builder=make_pipe,
            n_splits=random_splits,
            test_size=random_test_size,
            random_state=random_state,
        )
        report["models"][model_name] = {
            "features": spec,
            "logo_cv": logo_result,
            "random_split_cv": random_result,
        }
        pd.DataFrame(logo_result["per_fold"]).to_csv(output_dir / f"{model_name}_logo_per_fold.csv", index=False)
        pd.DataFrame(random_result["per_fold"]).to_csv(output_dir / f"{model_name}_random_per_fold.csv", index=False)

    descriptor = report["models"]["descriptor_model"]["logo_cv"]
    cond_only = report["models"]["condition_only"]["logo_cv"]
    cond_plus_sub = report["models"]["condition_plus_substrate_ohe"]["logo_cv"]
    report["comparisons"] = {
        "descriptor_minus_condition_only_spearman": descriptor["global_spearman"] - cond_only["global_spearman"],
        "descriptor_minus_condition_only_apyr": descriptor["global_apyr_mean"] - cond_only["global_apyr_mean"],
        "descriptor_minus_condition_plus_substrate_spearman": descriptor["global_spearman"] - cond_plus_sub["global_spearman"],
        "descriptor_minus_condition_plus_substrate_apyr": descriptor["global_apyr_mean"] - cond_plus_sub["global_apyr_mean"],
    }

    with (output_dir / "evaluation_report.json").open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
    return report


def _parse_int_list(text: str) -> list[int]:
    return [int(v.strip()) for v in text.split(",") if v.strip()]


def _parse_float_list(text: str) -> list[float]:
    return [float(v.strip()) for v in text.split(",") if v.strip()]


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Evaluate Chan-Lam descriptor model vs baselines.")
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=Path("data/HTE_db/literature/ChanLam_dataset_converted (2)_canonical.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/ml/chanlam_eval"),
    )
    parser.add_argument("--n-estimators", type=int, default=250)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--random-splits", type=int, default=5)
    parser.add_argument("--random-test-size", type=float, default=0.2)
    parser.add_argument(
        "--top-ks",
        type=str,
        default="1,3,5",
        help="Comma-separated top-k values for hit-rate metrics.",
    )
    parser.add_argument(
        "--yield-thresholds",
        type=str,
        default="50,70",
        help="Comma-separated yield thresholds for top-k hit rates.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    top_ks = _parse_int_list(args.top_ks)
    thresholds = _parse_float_list(args.yield_thresholds)
    report = evaluate_models(
        input_csv=args.input_csv,
        output_dir=args.output_dir,
        n_estimators=args.n_estimators,
        random_state=args.random_state,
        random_splits=args.random_splits,
        random_test_size=args.random_test_size,
        top_ks=top_ks,
        yield_thresholds=thresholds,
    )
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

