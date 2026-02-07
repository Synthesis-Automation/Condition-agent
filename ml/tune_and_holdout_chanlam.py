"""Descriptor tuning with ranking objective + fixed external holdout evaluation."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import argparse
import json
import re
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
from lightgbm import LGBMRanker
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.model_selection import LeaveOneGroupOut, ShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from chemtools.reagent.reagent_v2 import ReagentTaxonomyV2
from ml.chanlam_pipeline import (
    _CONDITION_COLUMNS,
    _NUMERIC_COLUMNS,
    _TEXT_COLUMNS,
    build_feature_table,
    load_chanlam_dataset,
)
from ml.evaluate_chanlam import _topk_metrics_for_fold, evaluate_predictions


_ROLE_HINT_BY_COLUMN = {
    "catalyst": "metal_catalyst",
    "base": "base",
    "solvent": "solvent",
}

_CONDITION_PROP_CAT_COLS = [
    "catalyst_metal",
    "catalyst_family",
    "base_family",
    "base_strength_band",
    "base_inorganic_band",
    "solvent_family",
    "solvent_class",
    "solvent_proticity",
    "solvent_polarity_band",
]

_CONDITION_PROP_NUM_COLS = [
    "is_cu_catalyst",
    "is_strong_base",
    "is_polar_solvent",
]

_METAL_PATTERNS = {
    "Cu": re.compile(r"\bcu\b|copper", re.I),
    "Pd": re.compile(r"\bpd\b|palladium", re.I),
    "Ni": re.compile(r"\bni\b|nickel", re.I),
    "Co": re.compile(r"\bco\b|cobalt", re.I),
    "Fe": re.compile(r"\bfe\b|iron", re.I),
    "Rh": re.compile(r"\brh\b|rhodium", re.I),
    "Ir": re.compile(r"\bir\b|iridium", re.I),
    "Ru": re.compile(r"\bru\b|ruthenium", re.I),
}


@dataclass(frozen=True)
class CandidateSpec:
    name: str
    model_type: str
    include_substrate_ohe: bool
    include_condition_props: bool
    params: Dict[str, Any]


class BaseModelAdapter:
    """Fit/predict interface for evaluators."""

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
        params: Dict[str, Any],
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
            random_state=random_state,
            **params,
        )

    @staticmethod
    def _group_relevance(y_sorted: np.ndarray, groups_sorted: np.ndarray) -> np.ndarray:
        relevance = np.zeros_like(y_sorted, dtype=np.int32)
        _, counts = np.unique(groups_sorted, return_counts=True)
        start = 0
        for count in counts:
            stop = start + int(count)
            y_grp = y_sorted[start:stop]
            ranks = pd.Series(y_grp).rank(method="average", pct=True).to_numpy(dtype=float)
            relevance[start:stop] = np.clip(np.rint(ranks * 20.0), 1, 20).astype(np.int32)
            start = stop
        return relevance

    def fit(self, X: pd.DataFrame, y: np.ndarray, groups: np.ndarray) -> "LGBMRankerAdapter":
        X_mat = self.preprocess.fit_transform(X)
        groups_arr = np.asarray(groups, dtype=str)
        order = np.argsort(groups_arr, kind="mergesort")
        X_sorted = X_mat[order]
        y_sorted = np.asarray(y, dtype=float)[order]
        g_sorted = groups_arr[order]
        _, group_sizes = np.unique(g_sorted, return_counts=True)
        relevance = self._group_relevance(y_sorted, g_sorted)
        self.model.fit(X_sorted, relevance, group=group_sizes.tolist())
        return self

    def predict(self, X: pd.DataFrame) -> np.ndarray:
        X_mat = self.preprocess.transform(X)
        return np.asarray(self.model.predict(X_mat), dtype=float)


@lru_cache(maxsize=1)
def _taxonomy() -> ReagentTaxonomyV2:
    return ReagentTaxonomyV2.from_path()


def _normalize_name(value: Any) -> str:
    if value is None:
        return "NA"
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return "NA"
    return text


def _infer_catalyst_metal(catalyst: str) -> str:
    for metal, pattern in _METAL_PATTERNS.items():
        if pattern.search(catalyst):
            return metal
    return "other_or_unknown"


def _infer_base_strength(base: str) -> str:
    text = base.lower()
    if text == "na":
        return "unknown"
    strong_tokens = ("otbu", "ome", "oet", "oh", "hydride", "amide", "dbu", "dbn", "lithium")
    medium_tokens = ("co3", "po4", "f", "carbonate", "phosphate", "acetate")
    weak_tokens = ("et3n", "dipea", "pyridine", "dmap", "amine", "lutidine")
    if any(tok in text for tok in strong_tokens):
        return "strong"
    if any(tok in text for tok in medium_tokens):
        return "medium"
    if any(tok in text for tok in weak_tokens):
        return "weak"
    return "unknown"


def _infer_base_inorganic(base: str) -> str:
    text = base.lower()
    if text == "na":
        return "unknown"
    tokens = ("k", "na", "cs", "li", "carbonate", "phosphate", "hydroxide", "fluoride")
    return "inorganic" if any(tok in text for tok in tokens) else "organic_or_unknown"


def _infer_solvent_props(solvent: str) -> tuple[str, str, str]:
    text = solvent.lower()
    if text == "na":
        return "unknown", "unknown", "unknown"
    if text in {"dce", "dcm", "ch2cl2"}:
        return "chlorinated_aprotic", "aprotic", "medium"
    if text in {"mecn", "dmf", "dma", "dmso", "nmp"}:
        return "polar_aprotic", "aprotic", "high"
    if text in {"meoh", "etoh", "ipa", "h2o", "water"}:
        return "protic", "protic", "high"
    if text in {"etOAc".lower(), "thf", "dioxane", "mtbe"}:
        return "medium_aprotic", "aprotic", "medium"
    return "other", "unknown", "unknown"


def _classify_family(value: str, role_hint: str) -> str:
    if value == "NA":
        return "NA"
    try:
        result = _taxonomy().classify({"name": value, "cas": None, "smiles": None})
    except Exception:
        return "unknown"
    if result is None:
        return "unknown"
    if result.role_id != role_hint:
        return "role_mismatch"
    return result.family_id


def add_condition_property_features(feat_df: pd.DataFrame) -> pd.DataFrame:
    """Add engineered condition-property descriptors."""
    out = feat_df.copy()
    for col in _CONDITION_COLUMNS:
        out[col] = out[col].map(_normalize_name)

    catalyst_norm = out["catalyst"].astype(str)
    base_norm = out["base"].astype(str)
    solvent_norm = out["solvent"].astype(str)

    out["catalyst_metal"] = catalyst_norm.map(_infer_catalyst_metal)
    out["base_strength_band"] = base_norm.map(_infer_base_strength)
    out["base_inorganic_band"] = base_norm.map(_infer_base_inorganic)
    solvent_props = solvent_norm.map(_infer_solvent_props)
    out["solvent_class"] = [x[0] for x in solvent_props]
    out["solvent_proticity"] = [x[1] for x in solvent_props]
    out["solvent_polarity_band"] = [x[2] for x in solvent_props]

    out["catalyst_family"] = catalyst_norm.map(lambda x: _classify_family(x, _ROLE_HINT_BY_COLUMN["catalyst"]))
    out["base_family"] = base_norm.map(lambda x: _classify_family(x, _ROLE_HINT_BY_COLUMN["base"]))
    out["solvent_family"] = solvent_norm.map(lambda x: _classify_family(x, _ROLE_HINT_BY_COLUMN["solvent"]))

    out["is_cu_catalyst"] = (out["catalyst_metal"] == "Cu").astype(float)
    out["is_strong_base"] = (out["base_strength_band"] == "strong").astype(float)
    out["is_polar_solvent"] = out["solvent_class"].isin(["polar_aprotic", "protic"]).astype(float)

    for col in _CONDITION_PROP_NUM_COLS:
        out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0.0)
    for col in _CONDITION_PROP_CAT_COLS:
        out[col] = out[col].fillna("unknown").astype(str)
    return out


def _derive_external_holdout_groups(feat_df: pd.DataFrame) -> list[str]:
    cond_key = feat_df[list(_CONDITION_COLUMNS)].astype(str).agg("|".join, axis=1)
    cond_counts = (
        pd.DataFrame(
            {
                "sulfonamide_smiles": feat_df["sulfonamide_smiles"].astype(str),
                "cond_key": cond_key,
            }
        )
        .groupby("sulfonamide_smiles")["cond_key"]
        .nunique()
        .sort_values()
    )
    max_cov = int(cond_counts.max())
    return sorted(cond_counts[cond_counts < max_cov].index.tolist())


def _build_adapter(
    *,
    spec: CandidateSpec,
    categorical_cols: list[str],
    numeric_cols: list[str],
    text_cols: list[str],
    random_state: int,
) -> BaseModelAdapter:
    if spec.model_type == "lgbm_rank":
        return LGBMRankerAdapter(
            categorical_cols=categorical_cols,
            numeric_cols=numeric_cols,
            text_cols=text_cols,
            params=spec.params,
            random_state=random_state,
        )

    transformers: list[tuple[str, Any, Any]] = []
    if numeric_cols:
        transformers.append(("num", StandardScaler(), numeric_cols))
    if categorical_cols:
        transformers.append(("cat", OneHotEncoder(handle_unknown="ignore"), categorical_cols))
    for col in text_cols:
        transformers.append((f"text_{col}", CountVectorizer(token_pattern=r"[^ ]+", binary=True), col))

    preprocess = ColumnTransformer(transformers=transformers, remainder="drop", sparse_threshold=0.3)
    model_kwargs = dict(random_state=random_state, n_jobs=-1, **spec.params)
    if spec.model_type == "et":
        model = ExtraTreesRegressor(**model_kwargs)
    else:
        model = RandomForestRegressor(**model_kwargs)
    return PipelineAdapter(Pipeline([("preprocess", preprocess), ("model", model)]))


def _evaluate_random_cv(
    X: pd.DataFrame,
    y: np.ndarray,
    groups: np.ndarray,
    *,
    builder: Any,
    n_splits: int,
    test_size: float,
    random_state: int,
) -> Dict[str, Any]:
    splitter = ShuffleSplit(n_splits=n_splits, test_size=test_size, random_state=random_state)
    rows: List[Dict[str, Any]] = []
    for fold, (train_idx, test_idx) in enumerate(splitter.split(X), start=1):
        model = builder()
        model.fit(X.iloc[train_idx], y[train_idx], groups[train_idx])
        pred = model.predict(X.iloc[test_idx])
        metrics = evaluate_predictions(y[test_idx], pred, groups[test_idx])
        rows.append(
            {
                "fold": fold,
                "spearman": metrics.spearman,
                "apyr_mean": metrics.apyr_mean,
                "apyr_std": metrics.apyr_std,
            }
        )
    df = pd.DataFrame(rows)
    return {
        "splits": n_splits,
        "test_size": test_size,
        "spearman_mean": float(df["spearman"].mean()),
        "spearman_std": float(df["spearman"].std(ddof=0)),
        "apyr_mean": float(df["apyr_mean"].mean()),
        "apyr_std": float(df["apyr_mean"].std(ddof=0)),
        "per_fold": rows,
    }


def _evaluate_logo_cv(
    X: pd.DataFrame,
    y: np.ndarray,
    groups: np.ndarray,
    *,
    builder: Any,
    top_ks: Sequence[int],
    yield_thresholds: Sequence[float],
) -> Dict[str, Any]:
    logo = LeaveOneGroupOut()
    per_fold: List[Dict[str, Any]] = []
    all_true: List[float] = []
    all_pred: List[float] = []
    all_group: List[str] = []
    for fold, (train_idx, test_idx) in enumerate(logo.split(X, y, groups=groups), start=1):
        model = builder()
        model.fit(X.iloc[train_idx], y[train_idx], groups[train_idx])
        pred = model.predict(X.iloc[test_idx])
        grp = str(groups[test_idx][0]) if len(test_idx) else "NA"
        metrics = evaluate_predictions(y[test_idx], pred, groups[test_idx])
        row = {
            "fold": fold,
            "group": grp,
            "rows": int(len(test_idx)),
            "spearman": metrics.spearman,
            "apyr_mean": metrics.apyr_mean,
            "apyr_std": metrics.apyr_std,
        }
        row.update(
            _topk_metrics_for_fold(
                y[test_idx],
                pred,
                top_ks=top_ks,
                yield_thresholds=yield_thresholds,
            )
        )
        per_fold.append(row)
        all_true.extend(y[test_idx].tolist())
        all_pred.extend(pred.tolist())
        all_group.extend(groups[test_idx].tolist())
    global_metrics = evaluate_predictions(np.array(all_true), np.array(all_pred), all_group)
    df = pd.DataFrame(per_fold)
    out = {
        "folds": int(len(per_fold)),
        "global_spearman": global_metrics.spearman,
        "global_apyr_mean": global_metrics.apyr_mean,
        "global_apyr_std": global_metrics.apyr_std,
        "global_apyr_n_groups": global_metrics.n_groups,
        "per_fold": per_fold,
    }
    for k in top_ks:
        out[f"top{k}_avg_best_percentile"] = float(df[f"top{k}_best_percentile"].mean())
        out[f"top{k}_avg_max_yield"] = float(df[f"top{k}_max_yield"].mean())
        for t in yield_thresholds:
            out[f"top{k}_hit_rate_ge_{int(t)}"] = float(df[f"top{k}_hit_ge_{int(t)}"].mean())
    return out


def _evaluate_external_holdout(
    X: pd.DataFrame,
    y: np.ndarray,
    groups: np.ndarray,
    *,
    holdout_groups: Sequence[str],
    builder: Any,
    top_ks: Sequence[int],
    yield_thresholds: Sequence[float],
) -> Dict[str, Any]:
    holdout_set = set(holdout_groups)
    test_mask = np.array([g in holdout_set for g in groups], dtype=bool)
    train_mask = ~test_mask
    model = builder()
    model.fit(X.iloc[train_mask], y[train_mask], groups[train_mask])
    pred = model.predict(X.iloc[test_mask])

    global_metrics = evaluate_predictions(y[test_mask], pred, groups[test_mask])
    rows: List[Dict[str, Any]] = []
    test_df = pd.DataFrame({"group": groups[test_mask], "y_true": y[test_mask], "y_pred": pred})
    for group, grp_df in test_df.groupby("group"):
        row = {"group": str(group), "rows": int(len(grp_df))}
        local = evaluate_predictions(grp_df["y_true"].to_numpy(), grp_df["y_pred"].to_numpy(), [group] * len(grp_df))
        row["spearman"] = local.spearman
        row["apyr_mean"] = local.apyr_mean
        row["apyr_std"] = local.apyr_std
        row.update(
            _topk_metrics_for_fold(
                grp_df["y_true"].to_numpy(),
                grp_df["y_pred"].to_numpy(),
                top_ks=top_ks,
                yield_thresholds=yield_thresholds,
            )
        )
        rows.append(row)
    df = pd.DataFrame(rows)
    out = {
        "train_rows": int(np.sum(train_mask)),
        "test_rows": int(np.sum(test_mask)),
        "test_groups": int(len(holdout_groups)),
        "global_spearman": global_metrics.spearman,
        "global_apyr_mean": global_metrics.apyr_mean,
        "global_apyr_std": global_metrics.apyr_std,
        "global_apyr_n_groups": global_metrics.n_groups,
        "per_group": rows,
    }
    for k in top_ks:
        out[f"top{k}_hit_rate_ge_50"] = float(df[f"top{k}_hit_ge_50"].mean())
        out[f"top{k}_hit_rate_ge_70"] = float(df[f"top{k}_hit_ge_70"].mean())
        out[f"top{k}_avg_best_percentile"] = float(df[f"top{k}_best_percentile"].mean())
    return out


def _candidate_specs(*, fast_mode: bool) -> list[CandidateSpec]:
    if fast_mode:
        return [
            CandidateSpec(
                name="et_desc_props_plus_sub_fast",
                model_type="et",
                include_substrate_ohe=True,
                include_condition_props=True,
                params={"n_estimators": 220, "min_samples_leaf": 1, "max_features": "sqrt", "max_depth": None},
            ),
            CandidateSpec(
                name="lgbm_rank_desc_props_plus_sub_fast",
                model_type="lgbm_rank",
                include_substrate_ohe=True,
                include_condition_props=True,
                params={"n_estimators": 180, "learning_rate": 0.06, "num_leaves": 31, "min_child_samples": 20},
            ),
            CandidateSpec(
                name="rf_desc_props_plus_sub_fast",
                model_type="rf",
                include_substrate_ohe=True,
                include_condition_props=True,
                params={"n_estimators": 240, "min_samples_leaf": 2, "max_features": "sqrt", "max_depth": None},
            ),
        ]
    return [
        CandidateSpec(
            name="rf_desc_props",
            model_type="rf",
            include_substrate_ohe=False,
            include_condition_props=True,
            params={"n_estimators": 400, "min_samples_leaf": 2, "max_features": 0.7, "max_depth": None},
        ),
        CandidateSpec(
            name="et_desc_props",
            model_type="et",
            include_substrate_ohe=False,
            include_condition_props=True,
            params={"n_estimators": 400, "min_samples_leaf": 1, "max_features": "sqrt", "max_depth": None},
        ),
        CandidateSpec(
            name="rf_desc_props_plus_sub",
            model_type="rf",
            include_substrate_ohe=True,
            include_condition_props=True,
            params={"n_estimators": 450, "min_samples_leaf": 1, "max_features": "sqrt", "max_depth": None},
        ),
        CandidateSpec(
            name="et_desc_props_plus_sub",
            model_type="et",
            include_substrate_ohe=True,
            include_condition_props=True,
            params={"n_estimators": 500, "min_samples_leaf": 1, "max_features": "sqrt", "max_depth": None},
        ),
        CandidateSpec(
            name="lgbm_rank_desc_props",
            model_type="lgbm_rank",
            include_substrate_ohe=False,
            include_condition_props=True,
            params={"n_estimators": 300, "learning_rate": 0.05, "num_leaves": 63, "min_child_samples": 20},
        ),
        CandidateSpec(
            name="lgbm_rank_desc_props_plus_sub",
            model_type="lgbm_rank",
            include_substrate_ohe=True,
            include_condition_props=True,
            params={"n_estimators": 350, "learning_rate": 0.05, "num_leaves": 63, "min_child_samples": 20},
        ),
    ]


def run_tuning_and_holdout(
    input_csv: Path,
    output_dir: Path,
    *,
    random_state: int,
    random_splits: int,
    random_test_size: float,
    top_ks: Sequence[int],
    yield_thresholds: Sequence[float],
    fast_mode: bool,
    max_finalists: int,
) -> Dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_df = load_chanlam_dataset(input_csv)
    feat_df = add_condition_property_features(build_feature_table(raw_df))
    y = feat_df["yield"].to_numpy(dtype=float)
    groups = feat_df["sulfonamide_smiles"].fillna("NA").astype(str).to_numpy()

    baseline_cat_cols = list(_CONDITION_COLUMNS) + ["sulfonamide_smiles", "boronic_smiles"]
    baseline_X = feat_df[baseline_cat_cols].copy()

    baseline_spec = CandidateSpec(
        name="baseline_condition_plus_substrate_ohe",
        model_type="rf",
        include_substrate_ohe=True,
        include_condition_props=False,
        params={"n_estimators": 350, "min_samples_leaf": 2, "max_features": "sqrt", "max_depth": None},
    )

    def baseline_builder() -> BaseModelAdapter:
        return _build_adapter(
            spec=baseline_spec,
            categorical_cols=baseline_cat_cols,
            numeric_cols=[],
            text_cols=[],
            random_state=random_state,
        )

    baseline_logo = _evaluate_logo_cv(
        baseline_X,
        y,
        groups,
        builder=baseline_builder,
        top_ks=top_ks,
        yield_thresholds=yield_thresholds,
    )

    descriptor_base_cols = list(_CONDITION_COLUMNS) + list(_NUMERIC_COLUMNS) + list(_TEXT_COLUMNS)
    descriptor_plus_sub_cols = descriptor_base_cols + ["sulfonamide_smiles", "boronic_smiles"]
    descriptor_prop_cols = list(_CONDITION_PROP_CAT_COLS) + list(_CONDITION_PROP_NUM_COLS)

    candidate_reports: List[Dict[str, Any]] = []
    for spec in _candidate_specs(fast_mode=fast_mode):
        use_cols = descriptor_plus_sub_cols if spec.include_substrate_ohe else descriptor_base_cols
        if spec.include_condition_props:
            use_cols = use_cols + descriptor_prop_cols
        X = feat_df[use_cols].copy()

        cat_cols = list(_CONDITION_COLUMNS)
        num_cols = list(_NUMERIC_COLUMNS)
        text_cols = list(_TEXT_COLUMNS)
        if spec.include_substrate_ohe:
            cat_cols = cat_cols + ["sulfonamide_smiles", "boronic_smiles"]
        if spec.include_condition_props:
            cat_cols = cat_cols + list(_CONDITION_PROP_CAT_COLS)
            num_cols = num_cols + list(_CONDITION_PROP_NUM_COLS)

        def make_builder(
            _spec: CandidateSpec = spec,
            _cat_cols: list[str] = cat_cols,
            _num_cols: list[str] = num_cols,
            _text_cols: list[str] = text_cols,
        ) -> BaseModelAdapter:
            return _build_adapter(
                spec=_spec,
                categorical_cols=_cat_cols,
                numeric_cols=_num_cols,
                text_cols=_text_cols,
                random_state=random_state,
            )

        random_eval = _evaluate_random_cv(
            X,
            y,
            groups,
            builder=make_builder,
            n_splits=random_splits,
            test_size=random_test_size,
            random_state=random_state,
        )
        candidate_reports.append(
            {
                "name": spec.name,
                "spec": {
                    "model_type": spec.model_type,
                    "include_substrate_ohe": spec.include_substrate_ohe,
                    "include_condition_props": spec.include_condition_props,
                    "params": dict(spec.params),
                },
                "random_split_cv": random_eval,
            }
        )

    candidate_reports = sorted(
        candidate_reports,
        key=lambda r: (r["random_split_cv"]["apyr_mean"], r["random_split_cv"]["spearman_mean"]),
        reverse=True,
    )
    finalists = candidate_reports[:max_finalists]

    final_with_logo: List[Dict[str, Any]] = []
    for entry in finalists:
        spec_dict = entry["spec"]
        spec = CandidateSpec(
            name=entry["name"],
            model_type=spec_dict["model_type"],
            include_substrate_ohe=spec_dict["include_substrate_ohe"],
            include_condition_props=spec_dict["include_condition_props"],
            params=spec_dict["params"],
        )
        use_cols = descriptor_plus_sub_cols if spec.include_substrate_ohe else descriptor_base_cols
        if spec.include_condition_props:
            use_cols = use_cols + descriptor_prop_cols
        X = feat_df[use_cols].copy()

        cat_cols = list(_CONDITION_COLUMNS)
        num_cols = list(_NUMERIC_COLUMNS)
        text_cols = list(_TEXT_COLUMNS)
        if spec.include_substrate_ohe:
            cat_cols = cat_cols + ["sulfonamide_smiles", "boronic_smiles"]
        if spec.include_condition_props:
            cat_cols = cat_cols + list(_CONDITION_PROP_CAT_COLS)
            num_cols = num_cols + list(_CONDITION_PROP_NUM_COLS)

        def make_builder(
            _spec: CandidateSpec = spec,
            _cat_cols: list[str] = cat_cols,
            _num_cols: list[str] = num_cols,
            _text_cols: list[str] = text_cols,
        ) -> BaseModelAdapter:
            return _build_adapter(
                spec=_spec,
                categorical_cols=_cat_cols,
                numeric_cols=_num_cols,
                text_cols=_text_cols,
                random_state=random_state,
            )

        logo_eval = _evaluate_logo_cv(
            X,
            y,
            groups,
            builder=make_builder,
            top_ks=top_ks,
            yield_thresholds=yield_thresholds,
        )
        row = dict(entry)
        row["logo_cv"] = logo_eval
        final_with_logo.append(row)

    best = sorted(
        final_with_logo,
        key=lambda r: (r["logo_cv"]["global_apyr_mean"], r["logo_cv"]["global_spearman"]),
        reverse=True,
    )[0]

    holdout_groups = _derive_external_holdout_groups(feat_df)
    (output_dir / "external_holdout_sulfonamides.txt").write_text("\n".join(holdout_groups), encoding="utf-8")

    baseline_holdout = _evaluate_external_holdout(
        baseline_X,
        y,
        groups,
        holdout_groups=holdout_groups,
        builder=baseline_builder,
        top_ks=top_ks,
        yield_thresholds=yield_thresholds,
    )

    best_spec = CandidateSpec(
        name=best["name"],
        model_type=best["spec"]["model_type"],
        include_substrate_ohe=best["spec"]["include_substrate_ohe"],
        include_condition_props=best["spec"]["include_condition_props"],
        params=best["spec"]["params"],
    )
    best_cols = descriptor_plus_sub_cols if best_spec.include_substrate_ohe else descriptor_base_cols
    if best_spec.include_condition_props:
        best_cols = best_cols + descriptor_prop_cols
    best_X = feat_df[best_cols].copy()
    best_cat_cols = list(_CONDITION_COLUMNS)
    best_num_cols = list(_NUMERIC_COLUMNS)
    best_text_cols = list(_TEXT_COLUMNS)
    if best_spec.include_substrate_ohe:
        best_cat_cols = best_cat_cols + ["sulfonamide_smiles", "boronic_smiles"]
    if best_spec.include_condition_props:
        best_cat_cols = best_cat_cols + list(_CONDITION_PROP_CAT_COLS)
        best_num_cols = best_num_cols + list(_CONDITION_PROP_NUM_COLS)

    def best_builder() -> BaseModelAdapter:
        return _build_adapter(
            spec=best_spec,
            categorical_cols=best_cat_cols,
            numeric_cols=best_num_cols,
            text_cols=best_text_cols,
            random_state=random_state,
        )

    best_holdout = _evaluate_external_holdout(
        best_X,
        y,
        groups,
        holdout_groups=holdout_groups,
        builder=best_builder,
        top_ks=top_ks,
        yield_thresholds=yield_thresholds,
    )

    report = {
        "input_csv": str(input_csv),
        "rows": int(len(feat_df)),
        "unique_sulfonamides": int(feat_df["sulfonamide_smiles"].nunique()),
        "external_holdout_groups_count": int(len(holdout_groups)),
        "external_holdout_list": holdout_groups,
        "feature_blocks": {
            "condition_property_categorical": list(_CONDITION_PROP_CAT_COLS),
            "condition_property_numeric": list(_CONDITION_PROP_NUM_COLS),
        },
        "fast_mode": bool(fast_mode),
        "max_finalists": int(max_finalists),
        "baseline_condition_plus_substrate_logo": {
            "global_spearman": baseline_logo["global_spearman"],
            "global_apyr_mean": baseline_logo["global_apyr_mean"],
            "global_apyr_std": baseline_logo["global_apyr_std"],
        },
        "tuning_candidates_random": candidate_reports,
        "tuning_finalists_logo": final_with_logo,
        "best_descriptor_model": best,
        "comparisons": {
            "best_descriptor_minus_baseline_logo_spearman": best["logo_cv"]["global_spearman"] - baseline_logo["global_spearman"],
            "best_descriptor_minus_baseline_logo_apyr": best["logo_cv"]["global_apyr_mean"] - baseline_logo["global_apyr_mean"],
        },
        "external_holdout_evaluation": {
            "baseline_condition_plus_substrate": baseline_holdout,
            "best_descriptor_model": best_holdout,
            "best_minus_baseline_spearman": best_holdout["global_spearman"] - baseline_holdout["global_spearman"],
            "best_minus_baseline_apyr": best_holdout["global_apyr_mean"] - baseline_holdout["global_apyr_mean"],
        },
    }

    with (output_dir / "tuning_and_holdout_report.json").open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
    pd.DataFrame(baseline_holdout["per_group"]).to_csv(output_dir / "baseline_external_holdout_per_group.csv", index=False)
    pd.DataFrame(best_holdout["per_group"]).to_csv(output_dir / "best_descriptor_external_holdout_per_group.csv", index=False)
    return report


def _parse_int_list(text: str) -> list[int]:
    return [int(v.strip()) for v in text.split(",") if v.strip()]


def _parse_float_list(text: str) -> list[float]:
    return [float(v.strip()) for v in text.split(",") if v.strip()]


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Tune descriptor model (with ranking) and run fixed external holdout.",
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=Path("data/HTE_db/literature/ChanLam_dataset_converted (2)_canonical.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/ml/chanlam_tune_holdout"),
    )
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--random-splits", type=int, default=3)
    parser.add_argument("--random-test-size", type=float, default=0.2)
    parser.add_argument("--top-ks", type=str, default="1,3,5")
    parser.add_argument("--yield-thresholds", type=str, default="50,70")
    parser.add_argument("--fast", action="store_true", help="Use a reduced candidate set for faster turnaround.")
    parser.add_argument("--max-finalists", type=int, default=4, help="Top-N candidates (from random CV) to run LOGO on.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    report = run_tuning_and_holdout(
        input_csv=args.input_csv,
        output_dir=args.output_dir,
        random_state=args.random_state,
        random_splits=args.random_splits,
        random_test_size=args.random_test_size,
        top_ks=_parse_int_list(args.top_ks),
        yield_thresholds=_parse_float_list(args.yield_thresholds),
        fast_mode=bool(args.fast),
        max_finalists=int(args.max_finalists),
    )
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
