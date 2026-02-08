"""Hybrid condition-recommendation PoC backtest on Chan-Lam CSV.

Hybrid score combines:
- ML prediction (descriptor + condition model)
- empirical condition prior from CSV (with shrinkage)
- nearest-reaction retrieval score from CSV precedents

Evaluation is offline Leave-One-Sulfonamide-Out (LOGO), so no wet-lab experiments
are required for this proof of concept.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from ml.chanlam_pipeline import build_feature_table, evaluate_predictions, load_chanlam_dataset
from ml.simple_chanlam_model import SIMPLE_NUMERIC_COLUMNS, SIMPLE_TEXT_COLUMNS
from ml.tune_and_holdout_chanlam import (
    _CONDITION_PROP_CAT_COLS,
    _CONDITION_PROP_NUM_COLS,
    add_condition_property_features,
)


_COND_COLS = ["catalyst", "base", "solvent"]
_REACTION_ID_COL = "reaction_smiles"
_GROUP_COL = "sulfonamide_smiles"


def _condition_key(df: pd.DataFrame) -> pd.Series:
    return df[_COND_COLS].fillna("NA").astype(str).agg("|".join, axis=1)


def _prepare_table(input_csv: Path) -> pd.DataFrame:
    raw = load_chanlam_dataset(input_csv)
    feat = build_feature_table(raw)
    feat = feat.copy()
    feat[_REACTION_ID_COL] = raw[_REACTION_ID_COL].fillna("NA").astype(str).to_numpy()
    return add_condition_property_features(feat)


def _ml_columns() -> Tuple[list[str], list[str], list[str]]:
    cat_cols = list(_COND_COLS) + list(_CONDITION_PROP_CAT_COLS)
    num_cols = list(SIMPLE_NUMERIC_COLUMNS) + list(_CONDITION_PROP_NUM_COLS)
    text_cols = list(SIMPLE_TEXT_COLUMNS)
    return cat_cols, num_cols, text_cols


def _build_ml_pipeline(*, n_estimators: int, random_state: int) -> Pipeline:
    cat_cols, num_cols, text_cols = _ml_columns()
    transformers: list[tuple[str, Any, Any]] = []
    if num_cols:
        transformers.append(("num", StandardScaler(), num_cols))
    if cat_cols:
        transformers.append(("cat", OneHotEncoder(handle_unknown="ignore"), cat_cols))
    for col in text_cols:
        transformers.append((f"text_{col}", CountVectorizer(token_pattern=r"[^ ]+", binary=True), col))
    preprocess = ColumnTransformer(transformers=transformers, remainder="drop", sparse_threshold=0.3)
    model = RandomForestRegressor(
        n_estimators=n_estimators,
        min_samples_leaf=2,
        random_state=random_state,
        n_jobs=-1,
    )
    return Pipeline([("preprocess", preprocess), ("model", model)])


def _reaction_vectorizer() -> tuple[ColumnTransformer, list[str], list[str]]:
    text_cols = list(SIMPLE_TEXT_COLUMNS)
    num_cols = list(SIMPLE_NUMERIC_COLUMNS)
    transformers: list[tuple[str, Any, Any]] = [("num", StandardScaler(), num_cols)]
    for col in text_cols:
        transformers.append((f"text_{col}", CountVectorizer(token_pattern=r"[^ ]+", binary=True), col))
    return ColumnTransformer(transformers=transformers, remainder="drop", sparse_threshold=0.3), text_cols, num_cols


def _build_prior(train_df: pd.DataFrame, *, shrinkage_m: float) -> Tuple[Dict[str, float], float]:
    global_mean = float(train_df["yield"].mean())
    stats = (
        train_df.assign(condition_key=_condition_key(train_df))
        .groupby("condition_key")["yield"]
        .agg(["mean", "count"])
        .reset_index()
    )
    prior: Dict[str, float] = {}
    for _, row in stats.iterrows():
        n = float(row["count"])
        mu = float(row["mean"])
        score = (n * mu + shrinkage_m * global_mean) / (n + shrinkage_m)
        prior[str(row["condition_key"])] = score
    return prior, global_mean


def _build_retrieval_assets(train_df: pd.DataFrame) -> Dict[str, Any]:
    # Unique reactions for similarity search.
    rxn_feature_cols = [_REACTION_ID_COL] + list(SIMPLE_TEXT_COLUMNS) + list(SIMPLE_NUMERIC_COLUMNS)
    unique_rxn = train_df[rxn_feature_cols].drop_duplicates(subset=[_REACTION_ID_COL]).reset_index(drop=True)

    vectorizer, _, _ = _reaction_vectorizer()
    rxn_matrix = vectorizer.fit_transform(unique_rxn)
    if not isinstance(rxn_matrix, csr_matrix):
        rxn_matrix = csr_matrix(rxn_matrix)

    yield_map = (
        train_df.assign(condition_key=_condition_key(train_df))
        .groupby([_REACTION_ID_COL, "condition_key"])["yield"]
        .mean()
    )
    return {
        "unique_rxn": unique_rxn,
        "vectorizer": vectorizer,
        "matrix": rxn_matrix,
        "yield_map": yield_map.to_dict(),
    }


def _retrieval_score_for_row(
    row: pd.Series,
    *,
    retrieval_assets: Dict[str, Any],
    condition_key: str,
    fallback: float,
    top_neighbors: int,
) -> float:
    unique_rxn = retrieval_assets["unique_rxn"]
    vectorizer: ColumnTransformer = retrieval_assets["vectorizer"]
    matrix: csr_matrix = retrieval_assets["matrix"]
    yield_map: Dict[tuple[str, str], float] = retrieval_assets["yield_map"]

    query_df = pd.DataFrame([row])[unique_rxn.columns]
    q = vectorizer.transform(query_df)
    sims = cosine_similarity(q, matrix).ravel()
    if sims.size == 0:
        return fallback

    idx = np.argsort(sims)[::-1][: max(1, int(top_neighbors))]
    weighted_sum = 0.0
    weight_total = 0.0
    for i in idx:
        sim = float(sims[i])
        if sim <= 0.0:
            continue
        rxn_id = str(unique_rxn.iloc[int(i)][_REACTION_ID_COL])
        y = yield_map.get((rxn_id, condition_key))
        if y is None:
            continue
        weighted_sum += sim * float(y)
        weight_total += sim
    if weight_total <= 0.0:
        return fallback
    return weighted_sum / weight_total


def _topk_hit_metrics(
    y_true: np.ndarray,
    scores: np.ndarray,
    reaction_ids: np.ndarray,
    *,
    top_ks: Sequence[int],
    thresholds: Sequence[float],
) -> Dict[str, float]:
    df = pd.DataFrame({"y": y_true, "score": scores, "reaction_id": reaction_ids})
    out: Dict[str, float] = {}
    grouped = list(df.groupby("reaction_id"))
    if not grouped:
        for k in top_ks:
            for t in thresholds:
                out[f"top{k}_hit_ge_{int(t)}"] = 0.0
        return out
    for k in top_ks:
        for t in thresholds:
            hits = []
            for _, grp in grouped:
                top = grp.sort_values("score", ascending=False).head(int(k))
                hits.append(float((top["y"] >= float(t)).any()))
            out[f"top{k}_hit_ge_{int(t)}"] = float(np.mean(hits))
    return out


def run_hybrid_poc(
    input_csv: Path,
    output_json: Path,
    *,
    n_estimators: int,
    random_state: int,
    shrinkage_m: float,
    w_ml: float,
    w_retrieval: float,
    w_prior: float,
    top_neighbors: int,
    top_ks: Sequence[int],
    thresholds: Sequence[float],
    max_folds: Optional[int],
) -> Dict[str, Any]:
    feat = _prepare_table(input_csv)
    feat = feat.copy()
    feat["condition_key"] = _condition_key(feat)

    cat_cols, num_cols, text_cols = _ml_columns()
    X_cols = cat_cols + num_cols + text_cols

    y = feat["yield"].to_numpy(dtype=float)
    groups = feat[_GROUP_COL].fillna("NA").astype(str).to_numpy()
    reaction_ids = feat[_REACTION_ID_COL].fillna("NA").astype(str).to_numpy()

    logo = LeaveOneGroupOut()
    fold_rows: List[Dict[str, Any]] = []
    pred_records: List[pd.DataFrame] = []

    for fold_idx, (train_idx, test_idx) in enumerate(logo.split(feat, y, groups=groups), start=1):
        if max_folds is not None and fold_idx > int(max_folds):
            break
        train_df = feat.iloc[train_idx].copy()
        test_df = feat.iloc[test_idx].copy()

        ml = _build_ml_pipeline(n_estimators=n_estimators, random_state=random_state)
        ml.fit(train_df[X_cols], train_df["yield"].to_numpy(dtype=float))
        ml_pred = np.asarray(ml.predict(test_df[X_cols]), dtype=float)

        prior_map, global_mean = _build_prior(train_df, shrinkage_m=shrinkage_m)
        prior_pred = np.asarray(
            [prior_map.get(str(c), global_mean) for c in test_df["condition_key"].astype(str)],
            dtype=float,
        )

        retrieval_assets = _build_retrieval_assets(train_df)
        retrieval_pred = np.asarray(
            [
                _retrieval_score_for_row(
                    row,
                    retrieval_assets=retrieval_assets,
                    condition_key=str(row["condition_key"]),
                    fallback=float(prior),
                    top_neighbors=top_neighbors,
                )
                for (_, row), prior in zip(test_df.iterrows(), prior_pred, strict=False)
            ],
            dtype=float,
        )

        hybrid_pred = w_ml * ml_pred + w_retrieval * retrieval_pred + w_prior * prior_pred

        fold_out = {
            "fold": int(fold_idx),
            "group": str(test_df[_GROUP_COL].iloc[0]) if not test_df.empty else "NA",
            "rows": int(len(test_df)),
        }

        for name, pred in (
            ("ml_only", ml_pred),
            ("prior_only", prior_pred),
            ("retrieval_only", retrieval_pred),
            ("hybrid", hybrid_pred),
        ):
            m = evaluate_predictions(test_df["yield"].to_numpy(dtype=float), pred, test_df[_GROUP_COL].astype(str).to_numpy())
            fold_out[f"{name}_spearman"] = float(m.spearman)
            fold_out[f"{name}_apyr"] = float(m.apyr_mean)
            topk = _topk_hit_metrics(
                test_df["yield"].to_numpy(dtype=float),
                pred,
                test_df[_REACTION_ID_COL].astype(str).to_numpy(),
                top_ks=top_ks,
                thresholds=thresholds,
            )
            for k, v in topk.items():
                fold_out[f"{name}_{k}"] = float(v)

        fold_rows.append(fold_out)
        pred_records.append(
            pd.DataFrame(
                {
                    "fold": fold_idx,
                    "group": test_df[_GROUP_COL].astype(str).to_numpy(),
                    "reaction_smiles": test_df[_REACTION_ID_COL].astype(str).to_numpy(),
                    "condition_key": test_df["condition_key"].astype(str).to_numpy(),
                    "y_true": test_df["yield"].to_numpy(dtype=float),
                    "ml_pred": ml_pred,
                    "prior_pred": prior_pred,
                    "retrieval_pred": retrieval_pred,
                    "hybrid_pred": hybrid_pred,
                }
            )
        )

    fold_df = pd.DataFrame(fold_rows)
    pred_df = pd.concat(pred_records, axis=0, ignore_index=True) if pred_records else pd.DataFrame()

    def _summary(prefix: str) -> Dict[str, Any]:
        out = {
            "spearman_mean": float(fold_df[f"{prefix}_spearman"].mean()),
            "spearman_std": float(fold_df[f"{prefix}_spearman"].std(ddof=0)),
            "apyr_mean": float(fold_df[f"{prefix}_apyr"].mean()),
            "apyr_std": float(fold_df[f"{prefix}_apyr"].std(ddof=0)),
        }
        for k in top_ks:
            for t in thresholds:
                col = f"{prefix}_top{k}_hit_ge_{int(t)}"
                out[f"top{k}_hit_ge_{int(t)}"] = float(fold_df[col].mean())
        return out

    report = {
        "input_csv": str(input_csv),
        "rows": int(len(feat)),
        "folds_run": int(len(fold_rows)),
        "settings": {
            "n_estimators": int(n_estimators),
            "random_state": int(random_state),
            "shrinkage_m": float(shrinkage_m),
            "weights": {"ml": float(w_ml), "retrieval": float(w_retrieval), "prior": float(w_prior)},
            "top_neighbors": int(top_neighbors),
            "top_ks": list(int(v) for v in top_ks),
            "thresholds": list(float(v) for v in thresholds),
            "max_folds": None if max_folds is None else int(max_folds),
        },
        "metrics": {
            "ml_only": _summary("ml_only"),
            "prior_only": _summary("prior_only"),
            "retrieval_only": _summary("retrieval_only"),
            "hybrid": _summary("hybrid"),
        },
        "comparisons": {
            "hybrid_minus_ml_apyr": float(_summary("hybrid")["apyr_mean"] - _summary("ml_only")["apyr_mean"]),
            "hybrid_minus_ml_spearman": float(_summary("hybrid")["spearman_mean"] - _summary("ml_only")["spearman_mean"]),
        },
    }

    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    if not pred_df.empty:
        pred_df.to_csv(output_json.with_suffix(".predictions.csv"), index=False)
    fold_df.to_csv(output_json.with_suffix(".per_fold.csv"), index=False)
    return report


def _parse_int_list(value: str) -> list[int]:
    return [int(v.strip()) for v in value.split(",") if v.strip()]


def _parse_float_list(value: str) -> list[float]:
    return [float(v.strip()) for v in value.split(",") if v.strip()]


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Hybrid CSV+ML recommendation PoC backtest.")
    p.add_argument(
        "--input-csv",
        type=Path,
        default=Path("data/HTE_db/literature/ChanLam_dataset_converted (2)_canonical.csv"),
    )
    p.add_argument(
        "--output-json",
        type=Path,
        default=Path("results/ml/chanlam_hybrid_poc/poc_report.json"),
    )
    p.add_argument("--n-estimators", type=int, default=180)
    p.add_argument("--random-state", type=int, default=42)
    p.add_argument("--shrinkage-m", type=float, default=20.0)
    p.add_argument("--w-ml", type=float, default=0.55)
    p.add_argument("--w-retrieval", type=float, default=0.30)
    p.add_argument("--w-prior", type=float, default=0.15)
    p.add_argument("--top-neighbors", type=int, default=8)
    p.add_argument("--top-ks", type=str, default="1,3")
    p.add_argument("--thresholds", type=str, default="50,70")
    p.add_argument("--max-folds", type=int, default=12, help="POC speed cap; use 0 for full LOGO.")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_parser().parse_args(argv)
    max_folds = None if int(args.max_folds) <= 0 else int(args.max_folds)
    report = run_hybrid_poc(
        input_csv=args.input_csv,
        output_json=args.output_json,
        n_estimators=args.n_estimators,
        random_state=args.random_state,
        shrinkage_m=args.shrinkage_m,
        w_ml=args.w_ml,
        w_retrieval=args.w_retrieval,
        w_prior=args.w_prior,
        top_neighbors=args.top_neighbors,
        top_ks=_parse_int_list(args.top_ks),
        thresholds=_parse_float_list(args.thresholds),
        max_folds=max_folds,
    )
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
