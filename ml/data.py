"""Dataset IO and split utilities for rebuilt ML system."""

from __future__ import annotations

from pathlib import Path
from typing import Iterator, Tuple

import pandas as pd
from sklearn.model_selection import LeaveOneGroupOut

from ml.chemistry import build_feature_table, load_chanlam_dataset
from ml.condition_features import add_condition_property_features
from ml.contracts import DatasetConfig


_COND_COLS = ["catalyst", "base", "solvent"]


def condition_key(df: pd.DataFrame) -> pd.Series:
    """Stable key for a condition tuple."""
    out = df.copy()
    for col in _COND_COLS:
        out[col] = out[col].fillna("NA").astype(str)
    return out[_COND_COLS].astype(str).agg("|".join, axis=1)


def load_chanlam_feature_table(
    cfg: DatasetConfig,
    *,
    with_condition_props: bool,
) -> pd.DataFrame:
    """Load canonical Chan-Lam CSV and build model-ready feature table."""
    raw_df = load_chanlam_dataset(Path(cfg.input_csv))
    feat_df = build_feature_table(raw_df)
    feat_df = feat_df.copy()
    feat_df[cfg.reaction_col] = raw_df[cfg.reaction_col].fillna("NA").astype(str).to_numpy()
    if with_condition_props:
        feat_df = add_condition_property_features(feat_df)
    return feat_df


def logo_splits(df: pd.DataFrame, *, group_col: str) -> Iterator[Tuple[int, pd.Index, pd.Index]]:
    """Generate LOGO split indices with fold numbers."""
    groups = df[group_col].fillna("NA").astype(str).to_numpy()
    logo = LeaveOneGroupOut()
    for fold_idx, (train_idx, test_idx) in enumerate(logo.split(df, groups=groups), start=1):
        yield fold_idx, pd.Index(train_idx), pd.Index(test_idx)


def assert_no_reaction_leakage(
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    *,
    reaction_col: str,
) -> None:
    """Raise if same reaction appears in both train and test partitions."""
    train_rxn = set(train_df[reaction_col].fillna("NA").astype(str).tolist())
    test_rxn = set(test_df[reaction_col].fillna("NA").astype(str).tolist())
    overlap = train_rxn.intersection(test_rxn)
    if overlap:
        example = sorted(overlap)[0]
        raise ValueError(f"Reaction leakage detected for '{example}' ({len(overlap)} overlapping reactions).")
