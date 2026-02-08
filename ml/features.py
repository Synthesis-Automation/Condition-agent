"""Feature profile definitions for the rebuilt ML pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

import pandas as pd

from ml.condition_features import CONDITION_PROP_CAT_COLS, CONDITION_PROP_NUM_COLS


_CORE_COND = ["catalyst", "base", "solvent"]
_CORE_TEXT = [
    "formed_motifs_tokens",
    "spectator_groups_tokens",
    "sulf_motif_tokens",
    "bor_motif_tokens",
]
_CORE_NUM = [
    "sulf_motif_count",
    "sulf_aryl_steric_max",
    "sulf_alkyl_steric_max",
    "sulf_aryl_electronic_avg",
    "bor_motif_count",
    "bor_aryl_steric_max",
    "bor_alkyl_steric_max",
    "bor_aryl_electronic_avg",
]

_BASE_COND = ["base"]
_BASE_TEXT = ["spectator_groups_tokens", "sulf_motif_tokens", "bor_motif_tokens"]
_BASE_NUM = ["sulf_motif_count", "bor_motif_count"]
_BASE_PROP_CAT = ["base_family", "base_strength_band", "base_inorganic_band"]
_BASE_PROP_NUM = ["is_strong_base"]


@dataclass(frozen=True)
class FeatureSpec:
    categorical_cols: List[str]
    numeric_cols: List[str]
    text_cols: List[str]

    @property
    def all_cols(self) -> List[str]:
        return list(self.categorical_cols) + list(self.numeric_cols) + list(self.text_cols)


def resolve_feature_spec(
    *,
    profile: str,
    with_condition_props: bool,
) -> FeatureSpec:
    """Resolve named feature profile into explicit columns."""
    if profile == "base_motif_spectator":
        cat = list(_BASE_COND)
        num = list(_BASE_NUM)
        text = list(_BASE_TEXT)
        if with_condition_props:
            cat += list(_BASE_PROP_CAT)
            num += list(_BASE_PROP_NUM)
        return FeatureSpec(cat, num, text)

    if profile == "core_full":
        cat = list(_CORE_COND)
        num = list(_CORE_NUM)
        text = list(_CORE_TEXT)
        if with_condition_props:
            cat += list(CONDITION_PROP_CAT_COLS)
            num += list(CONDITION_PROP_NUM_COLS)
        return FeatureSpec(cat, num, text)

    raise ValueError(f"Unsupported feature profile: {profile}")


def ensure_feature_columns(df: pd.DataFrame, spec: FeatureSpec) -> pd.DataFrame:
    """Ensure all columns in spec exist with sensible defaults."""
    out = df.copy()
    for col in spec.categorical_cols:
        if col not in out.columns:
            out[col] = "NA"
        out[col] = out[col].fillna("NA").astype(str)
    for col in spec.numeric_cols:
        if col not in out.columns:
            out[col] = 0.0
        out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0.0)
    for col in spec.text_cols:
        if col not in out.columns:
            out[col] = ""
        out[col] = out[col].fillna("").astype(str)
    return out


def feature_frame(df: pd.DataFrame, spec: FeatureSpec) -> pd.DataFrame:
    out = ensure_feature_columns(df, spec)
    return out[spec.all_cols].copy()


def build_condition_library(df: pd.DataFrame, *, condition_cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in condition_cols:
        out[col] = out[col].fillna("NA").astype(str)
    return out[condition_cols].drop_duplicates().reset_index(drop=True)


def profile_condition_columns(profile: str) -> list[str]:
    if profile == "base_motif_spectator":
        return ["base"]
    if profile == "core_full":
        return ["catalyst", "base", "solvent"]
    raise ValueError(f"Unsupported feature profile: {profile}")


def profile_defaults() -> Dict[str, str]:
    return {
        "profile": "core_full",
        "with_condition_props": "true",
    }
