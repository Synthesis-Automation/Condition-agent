"""Simplified Chan-Lam model with high-impact ranking/condition-property upgrades.

Core feature scope:
- Reactant motifs (text tokens + motif counts)
- Spectator groups / formed motifs from reaction featurization
- Electronic/steric summaries (aryl/alkyl steric + aryl electronics)
- Core conditions: catalyst, base, solvent

Optional high-impact upgrades:
- chemistry-aware condition-property features
- grouped ranking model (LambdaMART via LightGBM ranker)
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
from lightgbm import LGBMRanker
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.model_selection import LeaveOneGroupOut, ShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from chemtools.featurizers.unified import featurize_reaction
from ml.chanlam_pipeline import (
    DescriptorBuilder,
    _CONDITION_COLUMNS,
    _extract_reaction_substrates,
    _tokenize_text_field,
    build_feature_table,
    evaluate_predictions,
    load_chanlam_dataset,
)
from ml.tune_and_holdout_chanlam import (
    _CONDITION_PROP_CAT_COLS,
    _CONDITION_PROP_NUM_COLS,
    add_condition_property_features,
)


SIMPLE_CONDITION_COLUMNS = ["catalyst", "base", "solvent"]
SIMPLE_TEXT_COLUMNS = [
    "formed_motifs_tokens",
    "spectator_groups_tokens",
    "sulf_motif_tokens",
    "bor_motif_tokens",
]
SIMPLE_NUMERIC_COLUMNS = [
    "sulf_motif_count",
    "sulf_aryl_steric_max",
    "sulf_alkyl_steric_max",
    "sulf_aryl_electronic_avg",
    "bor_motif_count",
    "bor_aryl_steric_max",
    "bor_alkyl_steric_max",
    "bor_aryl_electronic_avg",
]

BASE_ONLY_CONDITION_COLUMNS = ["base"]
BASE_ONLY_TEXT_COLUMNS = [
    "spectator_groups_tokens",
    "sulf_motif_tokens",
    "bor_motif_tokens",
]
BASE_ONLY_NUMERIC_COLUMNS = [
    "sulf_motif_count",
    "bor_motif_count",
]
BASE_ONLY_PROP_CAT_COLUMNS = [
    "base_family",
    "base_strength_band",
    "base_inorganic_band",
]
BASE_ONLY_PROP_NUM_COLUMNS = [
    "is_strong_base",
]


class BaseModelAdapter:
    def fit(self, X: pd.DataFrame, y: np.ndarray, groups: np.ndarray) -> "BaseModelAdapter":
        raise NotImplementedError

    def predict(self, X: pd.DataFrame) -> np.ndarray:
        raise NotImplementedError


class PipelineAdapter(BaseModelAdapter):
    def __init__(self, pipeline: Pipeline) -> None:
        self.pipeline = pipeline

    def fit(self, X: pd.DataFrame, y: np.ndarray, groups: np.ndarray) -> "PipelineAdapter":
        self.pipeline.fit(X, y)
        return self

    def predict(self, X: pd.DataFrame) -> np.ndarray:
        return np.asarray(self.pipeline.predict(X), dtype=float)


class LGBMRankerAdapter(BaseModelAdapter):
    def __init__(
        self,
        *,
        categorical_cols: list[str],
        numeric_cols: list[str],
        text_cols: list[str],
        n_estimators: int,
        random_state: int,
    ) -> None:
        transformers: list[tuple[str, Any, Any]] = []
        if numeric_cols:
            transformers.append(("num", StandardScaler(), numeric_cols))
        if categorical_cols:
            transformers.append(("cat", OneHotEncoder(handle_unknown="ignore"), categorical_cols))
        for col in text_cols:
            transformers.append((f"text_{col}", CountVectorizer(token_pattern=r"[^ ]+", binary=True), col))
        self.preprocess = ColumnTransformer(transformers=transformers, remainder="drop", sparse_threshold=0.3)
        self.model = LGBMRanker(
            objective="lambdarank",
            metric="ndcg",
            n_estimators=n_estimators,
            learning_rate=0.06,
            num_leaves=31,
            min_child_samples=20,
            random_state=random_state,
            verbosity=-1,
        )

    @staticmethod
    def _group_relevance(y_sorted: np.ndarray, groups_sorted: np.ndarray) -> np.ndarray:
        relevance = np.zeros_like(y_sorted, dtype=np.int32)
        _, counts = np.unique(groups_sorted, return_counts=True)
        start = 0
        for count in counts:
            stop = start + int(count)
            grp_y = y_sorted[start:stop]
            grp_rank = pd.Series(grp_y).rank(method="average", pct=True).to_numpy(dtype=float)
            relevance[start:stop] = np.clip(np.rint(grp_rank * 20.0), 1, 20).astype(np.int32)
            start = stop
        return relevance

    def fit(self, X: pd.DataFrame, y: np.ndarray, groups: np.ndarray) -> "LGBMRankerAdapter":
        X_mat = self.preprocess.fit_transform(X)
        g = np.asarray(groups, dtype=str)
        order = np.argsort(g, kind="mergesort")
        X_sorted = X_mat[order]
        y_sorted = np.asarray(y, dtype=float)[order]
        g_sorted = g[order]
        _, group_sizes = np.unique(g_sorted, return_counts=True)
        relevance = self._group_relevance(y_sorted, g_sorted)
        self.model.fit(X_sorted, relevance, group=group_sizes.tolist())
        return self

    def predict(self, X: pd.DataFrame) -> np.ndarray:
        X_mat = self.preprocess.transform(X)
        return np.asarray(self.model.predict(X_mat), dtype=float)


def _build_model_adapter(
    *,
    model_type: str,
    n_estimators: int,
    random_state: int,
    categorical_cols: list[str],
    numeric_cols: list[str],
    text_cols: list[str],
) -> BaseModelAdapter:
    if model_type == "lgbm_rank":
        return LGBMRankerAdapter(
            categorical_cols=categorical_cols,
            numeric_cols=numeric_cols,
            text_cols=text_cols,
            n_estimators=n_estimators,
            random_state=random_state,
        )
    if model_type != "rf_reg":
        raise ValueError(f"Unsupported model_type: {model_type}")
    transformers: list[tuple[str, Any, Any]] = []
    if numeric_cols:
        transformers.append(("num", StandardScaler(), numeric_cols))
    if categorical_cols:
        transformers.append(("cat", OneHotEncoder(handle_unknown="ignore"), categorical_cols))
    for col in text_cols:
        transformers.append((f"text_{col}", CountVectorizer(token_pattern=r"[^ ]+", binary=True), col))
    preprocess = ColumnTransformer(transformers=transformers, remainder="drop", sparse_threshold=0.3)
    model = RandomForestRegressor(
        n_estimators=n_estimators,
        min_samples_leaf=2,
        random_state=random_state,
        n_jobs=-1,
    )
    return PipelineAdapter(Pipeline([("preprocess", preprocess), ("model", model)]))


def _ensure_full_condition_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in _CONDITION_COLUMNS:
        if col not in out.columns:
            out[col] = "NA"
    return out


def _feature_schema(
    with_condition_props: bool,
    feature_profile: str,
) -> tuple[list[str], list[str], list[str]]:
    if feature_profile == "base_motif_spectator":
        cat_cols = list(BASE_ONLY_CONDITION_COLUMNS)
        num_cols = list(BASE_ONLY_NUMERIC_COLUMNS)
        text_cols = list(BASE_ONLY_TEXT_COLUMNS)
        if with_condition_props:
            cat_cols += list(BASE_ONLY_PROP_CAT_COLUMNS)
            num_cols += list(BASE_ONLY_PROP_NUM_COLUMNS)
        return cat_cols, num_cols, text_cols

    if feature_profile != "core_full":
        raise ValueError(f"Unsupported feature_profile: {feature_profile}")

    cat_cols = list(SIMPLE_CONDITION_COLUMNS)
    num_cols = list(SIMPLE_NUMERIC_COLUMNS)
    text_cols = list(SIMPLE_TEXT_COLUMNS)
    if with_condition_props:
        cat_cols += list(_CONDITION_PROP_CAT_COLS)
        num_cols += list(_CONDITION_PROP_NUM_COLS)
    return cat_cols, num_cols, text_cols


def _prepare_feature_table(
    input_csv: Path,
    *,
    with_condition_props: bool,
) -> pd.DataFrame:
    raw_df = load_chanlam_dataset(input_csv)
    feat_df = build_feature_table(raw_df)
    if with_condition_props:
        feat_df = add_condition_property_features(_ensure_full_condition_columns(feat_df))
    return feat_df


def _feature_frame(
    feat_df: pd.DataFrame,
    *,
    with_condition_props: bool,
    feature_profile: str,
) -> pd.DataFrame:
    cat_cols, num_cols, text_cols = _feature_schema(with_condition_props, feature_profile)
    return feat_df[cat_cols + num_cols + text_cols].copy()


def evaluate_simple_model(
    input_csv: Path,
    *,
    model_type: str,
    with_condition_props: bool,
    feature_profile: str,
    n_estimators: int,
    random_state: int,
    random_splits: int,
    random_test_size: float,
) -> Dict[str, Any]:
    feat_df = _prepare_feature_table(input_csv, with_condition_props=with_condition_props)
    X = _feature_frame(
        feat_df,
        with_condition_props=with_condition_props,
        feature_profile=feature_profile,
    )
    y = feat_df["yield"].to_numpy(dtype=float)
    groups = feat_df["sulfonamide_smiles"].fillna("NA").astype(str).to_numpy()
    cat_cols, num_cols, text_cols = _feature_schema(with_condition_props, feature_profile)

    splitter = ShuffleSplit(n_splits=random_splits, test_size=random_test_size, random_state=random_state)
    random_rows: List[Dict[str, Any]] = []
    for fold, (train_idx, test_idx) in enumerate(splitter.split(X), start=1):
        model = _build_model_adapter(
            model_type=model_type,
            n_estimators=n_estimators,
            random_state=random_state,
            categorical_cols=cat_cols,
            numeric_cols=num_cols,
            text_cols=text_cols,
        )
        model.fit(X.iloc[train_idx], y[train_idx], groups[train_idx])
        pred = model.predict(X.iloc[test_idx])
        metrics = evaluate_predictions(y[test_idx], pred, groups[test_idx])
        random_rows.append(
            {
                "fold": fold,
                "spearman": metrics.spearman,
                "apyr_mean": metrics.apyr_mean,
                "apyr_std": metrics.apyr_std,
            }
        )
    random_df = pd.DataFrame(random_rows)

    logo = LeaveOneGroupOut()
    logo_rows: List[Dict[str, Any]] = []
    y_all: List[float] = []
    pred_all: List[float] = []
    group_all: List[str] = []
    for fold, (train_idx, test_idx) in enumerate(logo.split(X, y, groups=groups), start=1):
        model = _build_model_adapter(
            model_type=model_type,
            n_estimators=n_estimators,
            random_state=random_state,
            categorical_cols=cat_cols,
            numeric_cols=num_cols,
            text_cols=text_cols,
        )
        model.fit(X.iloc[train_idx], y[train_idx], groups[train_idx])
        pred = model.predict(X.iloc[test_idx])
        metrics = evaluate_predictions(y[test_idx], pred, groups[test_idx])
        grp_name = str(groups[test_idx][0]) if len(test_idx) else "NA"
        logo_rows.append(
            {
                "fold": fold,
                "group": grp_name,
                "rows": int(len(test_idx)),
                "spearman": metrics.spearman,
                "apyr_mean": metrics.apyr_mean,
                "apyr_std": metrics.apyr_std,
            }
        )
        y_all.extend(y[test_idx].tolist())
        pred_all.extend(pred.tolist())
        group_all.extend(groups[test_idx].tolist())
    logo_global = evaluate_predictions(np.array(y_all), np.array(pred_all), group_all)

    return {
        "input_csv": str(input_csv),
        "rows": int(len(feat_df)),
        "model_type": model_type,
        "feature_profile": feature_profile,
        "with_condition_props": bool(with_condition_props),
        "features": {"categorical": cat_cols, "numeric": num_cols, "text": text_cols},
        "random_split_cv": {
            "splits": int(random_splits),
            "test_size": float(random_test_size),
            "spearman_mean": float(random_df["spearman"].mean()),
            "spearman_std": float(random_df["spearman"].std(ddof=0)),
            "apyr_mean": float(random_df["apyr_mean"].mean()),
            "apyr_std": float(random_df["apyr_mean"].std(ddof=0)),
            "per_fold": random_rows,
        },
        "logo_cv": {
            "folds": int(len(logo_rows)),
            "global_spearman": logo_global.spearman,
            "global_apyr_mean": logo_global.apyr_mean,
            "global_apyr_std": logo_global.apyr_std,
            "global_apyr_n_groups": int(logo_global.n_groups),
            "per_fold": logo_rows,
        },
    }


def _reaction_text_tokens(reaction_smiles: str) -> Dict[str, str]:
    out = {"formed_motifs_tokens": "", "spectator_groups_tokens": ""}
    if not reaction_smiles or ">>" not in reaction_smiles:
        return out
    payload = featurize_reaction(reaction_smiles, options={"detailed": True, "motif_site_filter": "substituent"})
    aggregates = (payload or {}).get("aggregates") or {}
    formed = aggregates.get("formed_motifs_center") or aggregates.get("formed_motifs") or []
    spectators = aggregates.get("spectator_groups_ranked") or aggregates.get("spectator_groups_combined") or []
    out["formed_motifs_tokens"] = _tokenize_text_field(" / ".join(str(v) for v in formed))
    out["spectator_groups_tokens"] = _tokenize_text_field(" / ".join(str(v) for v in spectators))
    return out


def _normalize_condition_space(df: pd.DataFrame, condition_cols: Sequence[str]) -> pd.DataFrame:
    out = df.copy()
    for col in condition_cols:
        out[col] = out[col].fillna("NA").astype(str)
    out = out[list(condition_cols)].drop_duplicates().reset_index(drop=True)
    return out


def recommend_conditions_for_reaction(
    reaction_smiles: str,
    *,
    input_csv: Path,
    model_type: str,
    with_condition_props: bool,
    feature_profile: str,
    n_estimators: int,
    random_state: int,
    top_k: int,
) -> pd.DataFrame:
    feat_df = _prepare_feature_table(input_csv, with_condition_props=with_condition_props)
    X_train = _feature_frame(
        feat_df,
        with_condition_props=with_condition_props,
        feature_profile=feature_profile,
    )
    y_train = feat_df["yield"].to_numpy(dtype=float)
    groups = feat_df["sulfonamide_smiles"].fillna("NA").astype(str).to_numpy()
    cat_cols, num_cols, text_cols = _feature_schema(with_condition_props, feature_profile)

    model = _build_model_adapter(
        model_type=model_type,
        n_estimators=n_estimators,
        random_state=random_state,
        categorical_cols=cat_cols,
        numeric_cols=num_cols,
        text_cols=text_cols,
    )
    model.fit(X_train, y_train, groups)

    sulfonamide_smiles, boronic_smiles = _extract_reaction_substrates(reaction_smiles)
    if not sulfonamide_smiles or not boronic_smiles:
        raise ValueError("Could not parse sulfonamide/boronic reactants from reaction_smiles.")

    descriptors = DescriptorBuilder().build_row_descriptors(sulfonamide_smiles, boronic_smiles)
    raw_df = load_chanlam_dataset(input_csv)
    condition_space = _normalize_condition_space(raw_df, [c for c in SIMPLE_CONDITION_COLUMNS if c in cat_cols])
    rxn_tokens = _reaction_text_tokens(reaction_smiles)

    rows: List[Dict[str, Any]] = []
    for _, cond in condition_space.iterrows():
        row: Dict[str, Any] = {
            "catalyst": str(cond.get("catalyst", "NA")),
            "base": str(cond.get("base", "NA")),
            "solvent": str(cond.get("solvent", "NA")),
            "formed_motifs_tokens": rxn_tokens["formed_motifs_tokens"],
            "spectator_groups_tokens": rxn_tokens["spectator_groups_tokens"],
            "sulf_motif_tokens": str(descriptors.get("sulf_motif_tokens", "")),
            "bor_motif_tokens": str(descriptors.get("bor_motif_tokens", "")),
            "sulf_motif_count": float(descriptors.get("sulf_motif_count", 0.0)),
            "sulf_aryl_steric_max": float(descriptors.get("sulf_aryl_steric_max", 0.0)),
            "sulf_alkyl_steric_max": float(descriptors.get("sulf_alkyl_steric_max", 0.0)),
            "sulf_aryl_electronic_avg": float(descriptors.get("sulf_aryl_electronic_avg", 5.0)),
            "bor_motif_count": float(descriptors.get("bor_motif_count", 0.0)),
            "bor_aryl_steric_max": float(descriptors.get("bor_aryl_steric_max", 0.0)),
            "bor_alkyl_steric_max": float(descriptors.get("bor_alkyl_steric_max", 0.0)),
            "bor_aryl_electronic_avg": float(descriptors.get("bor_aryl_electronic_avg", 5.0)),
        }
        rows.append(row)
    candidate_df = pd.DataFrame.from_records(rows)

    if with_condition_props:
        candidate_df = add_condition_property_features(_ensure_full_condition_columns(candidate_df))

    X_cand = candidate_df[cat_cols + num_cols + text_cols].copy()
    candidate_df["predicted_yield"] = model.predict(X_cand)
    candidate_df["reaction_smiles"] = reaction_smiles
    candidate_df["sulfonamide_smiles"] = sulfonamide_smiles
    candidate_df["boronic_smiles"] = boronic_smiles
    if feature_profile == "base_motif_spectator":
        candidate_df = candidate_df.drop(columns=["catalyst", "solvent"], errors="ignore")
    return candidate_df.sort_values("predicted_yield", ascending=False).reset_index(drop=True).head(max(1, int(top_k)))


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Simplified Chan-Lam model and condition recommendation.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    eval_p = sub.add_parser("evaluate", help="Evaluate simplified model on Chan-Lam dataset.")
    eval_p.add_argument(
        "--input-csv",
        type=Path,
        default=Path("data/HTE_db/literature/ChanLam_dataset_converted (2)_canonical.csv"),
    )
    eval_p.add_argument(
        "--output-json",
        type=Path,
        default=Path("results/ml/chanlam_simple/evaluation_report.json"),
    )
    eval_p.add_argument("--model-type", type=str, default="rf_reg", choices=["rf_reg", "lgbm_rank"])
    eval_p.add_argument(
        "--feature-profile",
        type=str,
        default="core_full",
        choices=["core_full", "base_motif_spectator"],
    )
    eval_p.add_argument(
        "--without-condition-props",
        action="store_true",
        help="Disable engineered condition-property features.",
    )
    eval_p.add_argument("--n-estimators", type=int, default=220)
    eval_p.add_argument("--random-state", type=int, default=42)
    eval_p.add_argument("--random-splits", type=int, default=5)
    eval_p.add_argument("--random-test-size", type=float, default=0.2)

    rec_p = sub.add_parser("recommend", help="Recommend top conditions for a new reaction.")
    rec_p.add_argument("--reaction-smiles", required=True, type=str)
    rec_p.add_argument(
        "--input-csv",
        type=Path,
        default=Path("data/HTE_db/literature/ChanLam_dataset_converted (2)_canonical.csv"),
    )
    rec_p.add_argument(
        "--output-csv",
        type=Path,
        default=Path("results/ml/chanlam_simple/new_reaction_top_conditions.csv"),
    )
    rec_p.add_argument("--model-type", type=str, default="rf_reg", choices=["rf_reg", "lgbm_rank"])
    rec_p.add_argument(
        "--feature-profile",
        type=str,
        default="core_full",
        choices=["core_full", "base_motif_spectator"],
    )
    rec_p.add_argument(
        "--without-condition-props",
        action="store_true",
        help="Disable engineered condition-property features.",
    )
    rec_p.add_argument("--n-estimators", type=int, default=220)
    rec_p.add_argument("--random-state", type=int, default=42)
    rec_p.add_argument("--top-k", type=int, default=10)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    with_condition_props = not bool(args.without_condition_props)

    if args.cmd == "evaluate":
        report = evaluate_simple_model(
            input_csv=args.input_csv,
            model_type=args.model_type,
            with_condition_props=with_condition_props,
            feature_profile=args.feature_profile,
            n_estimators=args.n_estimators,
            random_state=args.random_state,
            random_splits=args.random_splits,
            random_test_size=args.random_test_size,
        )
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(json.dumps(report, indent=2))
        return 0

    if args.cmd == "recommend":
        out = recommend_conditions_for_reaction(
            args.reaction_smiles,
            input_csv=args.input_csv,
            model_type=args.model_type,
            with_condition_props=with_condition_props,
            feature_profile=args.feature_profile,
            n_estimators=args.n_estimators,
            random_state=args.random_state,
            top_k=args.top_k,
        )
        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.output_csv, index=False)
        print(out.to_string(index=False))
        return 0

    parser.error("Unknown command")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
