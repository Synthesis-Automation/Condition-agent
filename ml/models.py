"""Model building blocks for two-stage condition recommendation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict

import numpy as np
import pandas as pd
from lightgbm import LGBMRanker
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from ml.features import FeatureSpec, feature_frame


def _condition_key_from_columns(df: pd.DataFrame, condition_cols: list[str]) -> pd.Series:
    out = df.copy()
    for col in condition_cols:
        out[col] = out[col].fillna("NA").astype(str)
    return out[condition_cols].astype(str).agg("|".join, axis=1)


@dataclass
class ConditionPriorRanker:
    """Fast stage-1 ranker based on shrinkage-smoothed condition priors."""

    condition_cols: list[str]
    shrinkage_m: float = 20.0
    prior_map: Dict[str, float] | None = None
    global_mean: float = 0.0

    def fit(self, df: pd.DataFrame, y_col: str = "yield") -> "ConditionPriorRanker":
        work = df.copy()
        work["condition_key"] = _condition_key_from_columns(work, self.condition_cols)
        self.global_mean = float(pd.to_numeric(work[y_col], errors="coerce").fillna(0.0).mean())
        stats = work.groupby("condition_key")[y_col].agg(["mean", "count"]).reset_index()
        out: Dict[str, float] = {}
        for _, row in stats.iterrows():
            n = float(row["count"])
            mu = float(row["mean"])
            out[str(row["condition_key"])] = float((n * mu + self.shrinkage_m * self.global_mean) / (n + self.shrinkage_m))
        self.prior_map = out
        return self

    def score(self, df: pd.DataFrame) -> np.ndarray:
        if self.prior_map is None:
            raise RuntimeError("ConditionPriorRanker is not fit.")
        keys = _condition_key_from_columns(df, self.condition_cols).astype(str).tolist()
        return np.asarray([self.prior_map.get(k, self.global_mean) for k in keys], dtype=float)


class Stage2Model:
    """Stage-2 descriptor model (RF regression or LGBM ranking)."""

    def __init__(
        self,
        *,
        model_type: str,
        feature_spec: FeatureSpec,
        n_estimators: int,
        random_state: int,
    ) -> None:
        self.model_type = model_type
        self.feature_spec = feature_spec
        self.n_estimators = int(n_estimators)
        self.random_state = int(random_state)
        self.preprocess = self._build_preprocess(feature_spec)
        self.pipeline: Pipeline | None = None
        self.ranker: LGBMRanker | None = None

        if model_type == "rf_reg":
            rf = RandomForestRegressor(
                n_estimators=self.n_estimators,
                min_samples_leaf=2,
                random_state=self.random_state,
                n_jobs=-1,
            )
            self.pipeline = Pipeline([("preprocess", self.preprocess), ("model", rf)])
        elif model_type == "lgbm_rank":
            self.ranker = LGBMRanker(
                objective="lambdarank",
                metric="ndcg",
                n_estimators=self.n_estimators,
                learning_rate=0.06,
                num_leaves=31,
                min_child_samples=20,
                random_state=self.random_state,
                verbosity=-1,
            )
        else:
            raise ValueError(f"Unsupported model type: {model_type}")

    @staticmethod
    def _build_preprocess(spec: FeatureSpec) -> ColumnTransformer:
        transformers: list[tuple[str, Any, Any]] = []
        if spec.numeric_cols:
            transformers.append(("num", StandardScaler(), spec.numeric_cols))
        if spec.categorical_cols:
            transformers.append(("cat", OneHotEncoder(handle_unknown="ignore"), spec.categorical_cols))
        for col in spec.text_cols:
            transformers.append((f"text_{col}", CountVectorizer(token_pattern=r"[^ ]+", binary=True), col))
        return ColumnTransformer(transformers=transformers, remainder="drop", sparse_threshold=0.3)

    @staticmethod
    def _group_relevance(y_sorted: np.ndarray, groups_sorted: np.ndarray) -> np.ndarray:
        relevance = np.zeros_like(y_sorted, dtype=np.int32)
        _, counts = np.unique(groups_sorted, return_counts=True)
        start = 0
        for count in counts:
            stop = start + int(count)
            grp_y = y_sorted[start:stop]
            ranks = pd.Series(grp_y).rank(method="average", pct=True).to_numpy(dtype=float)
            relevance[start:stop] = np.clip(np.rint(ranks * 20.0), 1, 20).astype(np.int32)
            start = stop
        return relevance

    def fit(self, df: pd.DataFrame, y: np.ndarray, groups: np.ndarray) -> "Stage2Model":
        X = feature_frame(df, self.feature_spec)
        if self.model_type == "rf_reg":
            assert self.pipeline is not None
            self.pipeline.fit(X, y)
            return self

        assert self.ranker is not None
        X_mat = self.preprocess.fit_transform(X)
        groups_arr = np.asarray(groups, dtype=str)
        order = np.argsort(groups_arr, kind="mergesort")
        X_sorted = X_mat[order]
        y_sorted = np.asarray(y, dtype=float)[order]
        g_sorted = groups_arr[order]
        _, group_sizes = np.unique(g_sorted, return_counts=True)
        rel = self._group_relevance(y_sorted, g_sorted)
        self.ranker.fit(X_sorted, rel, group=group_sizes.tolist())
        return self

    def predict(self, df: pd.DataFrame) -> np.ndarray:
        X = feature_frame(df, self.feature_spec)
        if self.model_type == "rf_reg":
            assert self.pipeline is not None
            return np.asarray(self.pipeline.predict(X), dtype=float)

        assert self.ranker is not None
        X_mat = self.preprocess.transform(X)
        return np.asarray(self.ranker.predict(X_mat), dtype=float)
