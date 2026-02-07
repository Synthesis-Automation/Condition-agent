"""
Training pipeline for Chan-Lam condition selection with chemistry-aware descriptors.

This module builds features from:
- Taxonomy/motif outputs (reactant motifs, formed motifs, spectator groups)
- Functional-group spectators
- Electronic/steric summaries from chemtools featurizers
- Condition variables (catalyst/base/solvent/reagent slots)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import argparse
import json
import math
import re
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from joblib import dump
from scipy.stats import spearmanr
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.model_selection import LeaveOneGroupOut, ShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from chemtools.featurizers.formatters.molecule import build_molecule_bundle
from chemtools.util.functional_groups import get_functional_groups


_BORON_TOKEN_RE = re.compile(r"B\(|\[B")
_TOKEN_SPLIT_RE = re.compile(r"[|/,;]+")

_CONDITION_COLUMNS = [
    "catalyst",
    "ligand",
    "base",
    "acid",
    "oxidant",
    "reductant",
    "additive",
    "condensation_agent",
    "other_reagent",
    "solvent",
]

_TEXT_COLUMNS = [
    "formed_motifs_tokens",
    "spectator_groups_tokens",
    "sulf_motif_tokens",
    "sulf_fg_tokens",
    "bor_motif_tokens",
    "bor_fg_tokens",
]

_NUMERIC_COLUMNS = [
    "sulf_motif_count",
    "sulf_aryl_steric_max",
    "sulf_alkyl_steric_max",
    "sulf_aryl_electronic_avg",
    "sulf_aromatic_ring_count",
    "sulf_heavy_atom_count",
    "sulf_molecular_weight",
    "sulf_logP",
    "sulf_TPSA",
    "sulf_HBA",
    "sulf_HBD",
    "sulf_rotatable_bonds",
    "sulf_fraction_csp3",
    "bor_motif_count",
    "bor_aryl_steric_max",
    "bor_alkyl_steric_max",
    "bor_aryl_electronic_avg",
    "bor_aromatic_ring_count",
    "bor_heavy_atom_count",
    "bor_molecular_weight",
    "bor_logP",
    "bor_TPSA",
    "bor_HBA",
    "bor_HBD",
    "bor_rotatable_bonds",
    "bor_fraction_csp3",
]


@dataclass
class EvaluationResult:
    """Container for split-wise evaluation metrics."""

    spearman: float
    apyr_mean: float
    apyr_std: float
    n_groups: int


def _tokenize_text_field(value: Any) -> str:
    """Normalize motif/spectator text fields into CountVectorizer-ready tokens."""
    if value is None:
        return ""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    tokens = [tok.strip().replace(" ", "_") for tok in _TOKEN_SPLIT_RE.split(text) if tok.strip()]
    if not tokens:
        return ""
    return " ".join(sorted(set(tokens)))


def _extract_reaction_substrates(reaction_smiles: str) -> Tuple[Optional[str], Optional[str]]:
    """Extract sulfonamide and boron partner SMILES from reaction SMILES."""
    if not isinstance(reaction_smiles, str) or ">>" not in reaction_smiles:
        return None, None
    left = reaction_smiles.split(">>", 1)[0]
    reactants = [r.strip() for r in left.split(".") if r.strip()]
    if not reactants:
        return None, None

    boronic = None
    non_boron: List[str] = []
    for token in reactants:
        if _BORON_TOKEN_RE.search(token):
            boronic = token if boronic is None else boronic
        else:
            non_boron.append(token)

    sulfonamide = None
    for token in non_boron:
        if "S(=O)(=O)" in token and "N" in token:
            sulfonamide = token
            break
    if sulfonamide is None and non_boron:
        sulfonamide = non_boron[0]
    return sulfonamide, boronic


def _collect_motif_ids(bundle: Dict[str, Any]) -> List[str]:
    motif_ids: List[str] = []
    for key in ("motifs", "context_motifs"):
        for entry in bundle.get(key, []) or []:
            if not isinstance(entry, dict):
                continue
            motif_id = str(entry.get("compound_id") or entry.get("id") or "").strip()
            if motif_id:
                motif_ids.append(motif_id)
    return sorted(set(motif_ids))


def _extract_score_values(payload: Any) -> List[float]:
    """Extract nested score values from steric/electronic payload entries."""
    scores: List[float] = []
    if payload is None:
        return scores
    if isinstance(payload, dict):
        value = payload.get("score_0_10")
        if isinstance(value, (int, float)):
            scores.append(float(value))
        return scores
    if isinstance(payload, list):
        for item in payload:
            scores.extend(_extract_score_values(item))
    return scores


def _safe_mean(values: Sequence[float], default: float = 0.0) -> float:
    if not values:
        return default
    return float(np.mean(values))


def _safe_max(values: Sequence[float], default: float = 0.0) -> float:
    if not values:
        return default
    return float(np.max(values))


class DescriptorBuilder:
    """Cache-aware molecular descriptor builder for sulfonamide/boronic partners."""

    def __init__(self) -> None:
        self._cache: Dict[str, Dict[str, Any]] = {}

    def _describe_smiles(self, smiles: Optional[str], prefix: str) -> Dict[str, Any]:
        if not smiles:
            return self._empty_descriptor(prefix)
        if smiles in self._cache:
            cached = self._cache[smiles]
            return {f"{prefix}_{k}": v for k, v in cached.items()}

        bundle = build_molecule_bundle(smiles, options={"include_rdkit": True})
        motif_ids = _collect_motif_ids(bundle)
        fg_tokens = sorted(get_functional_groups(smiles))

        aryl_steric_scores: List[float] = []
        alkyl_steric_scores: List[float] = []
        for entry in (bundle.get("steric", {}) or {}).get("aryl", []) or []:
            aryl_steric_scores.extend(_extract_score_values(entry.get("result")))
        for entry in (bundle.get("steric", {}) or {}).get("alkyl", []) or []:
            alkyl_steric_scores.extend(_extract_score_values(entry.get("result")))

        electronic_scores: List[float] = []
        for entry in (bundle.get("electronics", {}) or {}).get("aryl", []) or []:
            electronic_scores.extend(_extract_score_values(entry.get("result")))

        rdkit_props = bundle.get("rdkit_properties", {}) or {}
        aryl_analysis = bundle.get("aryl_analysis", {}) or {}

        descriptor = {
            "motif_tokens": " ".join(motif_ids),
            "fg_tokens": " ".join(fg_tokens),
            "motif_count": float(len(motif_ids)),
            "aryl_steric_max": _safe_max(aryl_steric_scores, default=0.0),
            "alkyl_steric_max": _safe_max(alkyl_steric_scores, default=0.0),
            "aryl_electronic_avg": _safe_mean(electronic_scores, default=5.0),
            "aromatic_ring_count": float(aryl_analysis.get("aromatic_ring_count", 0.0) or 0.0),
            "heavy_atom_count": float(rdkit_props.get("heavy_atom_count", 0.0) or 0.0),
            "molecular_weight": float(rdkit_props.get("molecular_weight", 0.0) or 0.0),
            "logP": float(rdkit_props.get("logP", 0.0) or 0.0),
            "TPSA": float(rdkit_props.get("TPSA", 0.0) or 0.0),
            "HBA": float(rdkit_props.get("HBA", 0.0) or 0.0),
            "HBD": float(rdkit_props.get("HBD", 0.0) or 0.0),
            "rotatable_bonds": float(rdkit_props.get("rotatable_bonds", 0.0) or 0.0),
            "fraction_csp3": float(rdkit_props.get("fraction_csp3", 0.0) or 0.0),
        }
        self._cache[smiles] = descriptor
        return {f"{prefix}_{k}": v for k, v in descriptor.items()}

    @staticmethod
    def _empty_descriptor(prefix: str) -> Dict[str, Any]:
        out = {}
        for name in (
            "motif_tokens",
            "fg_tokens",
            "motif_count",
            "aryl_steric_max",
            "alkyl_steric_max",
            "aryl_electronic_avg",
            "aromatic_ring_count",
            "heavy_atom_count",
            "molecular_weight",
            "logP",
            "TPSA",
            "HBA",
            "HBD",
            "rotatable_bonds",
            "fraction_csp3",
        ):
            value: Any = "" if name.endswith("_tokens") else 0.0
            out[f"{prefix}_{name}"] = value
        out[f"{prefix}_aryl_electronic_avg"] = 5.0
        return out

    def build_row_descriptors(
        self,
        sulfonamide_smiles: Optional[str],
        boronic_smiles: Optional[str],
    ) -> Dict[str, Any]:
        result = {}
        result.update(self._describe_smiles(sulfonamide_smiles, "sulf"))
        result.update(self._describe_smiles(boronic_smiles, "bor"))
        return result


def load_chanlam_dataset(csv_path: Path) -> pd.DataFrame:
    """Load and normalize Chan-Lam dataset table."""
    df = pd.read_csv(csv_path)
    if "yield" not in df.columns:
        raise ValueError("Input CSV must contain a 'yield' column.")
    if "reaction_smiles" not in df.columns:
        raise ValueError("Input CSV must contain a 'reaction_smiles' column.")

    parsed = df["reaction_smiles"].astype(str).map(_extract_reaction_substrates)
    df["sulfonamide_smiles"] = [s for s, _ in parsed]
    df["boronic_smiles"] = [b for _, b in parsed]
    return df


def build_feature_table(df: pd.DataFrame) -> pd.DataFrame:
    """Build model-ready feature table with chemistry descriptors."""
    builder = DescriptorBuilder()
    records: List[Dict[str, Any]] = []
    for _, row in df.iterrows():
        record: Dict[str, Any] = {
            "yield": float(row["yield"]),
            "sulfonamide_smiles": row.get("sulfonamide_smiles"),
            "boronic_smiles": row.get("boronic_smiles"),
        }
        for col in _CONDITION_COLUMNS:
            value = row.get(col)
            record[col] = "NA" if pd.isna(value) else str(value)

        record["formed_motifs_tokens"] = _tokenize_text_field(row.get("formed_motifs"))
        record["spectator_groups_tokens"] = _tokenize_text_field(row.get("spectator_groups"))

        descriptor = builder.build_row_descriptors(
            record["sulfonamide_smiles"],
            record["boronic_smiles"],
        )
        record.update(descriptor)
        records.append(record)
    feat_df = pd.DataFrame.from_records(records)
    for col in _NUMERIC_COLUMNS:
        feat_df[col] = pd.to_numeric(feat_df[col], errors="coerce").fillna(0.0)
    for col in _TEXT_COLUMNS:
        feat_df[col] = feat_df[col].fillna("").astype(str)
    return feat_df


def _percentile_rank(values: np.ndarray, score: float) -> float:
    """Percentile rank (weak) in [0, 100]."""
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
    """
    Compute Average Percentile Yield Ranking (APYR).

    For each group, select rows within `top_within` predicted yield units
    from the group's max prediction and compute percentile rank of the
    corresponding experimental yields.
    """
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
        ranks = [_percentile_rank(exp, val) for val in chosen_exp]
        apyr_values.append(float(np.mean(ranks)))

    if not apyr_values:
        return 0.0, 0.0, 0
    return float(np.mean(apyr_values)), float(np.std(apyr_values)), len(apyr_values)


def evaluate_predictions(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    groups: Sequence[str],
) -> EvaluationResult:
    """Compute Spearman and APYR metrics for a prediction set."""
    rho = spearmanr(y_true, y_pred).correlation
    spearman = 0.0 if rho is None or math.isnan(rho) else float(rho)
    apyr_mean, apyr_std, n_groups = compute_apyr(y_true, y_pred, groups)
    return EvaluationResult(
        spearman=spearman,
        apyr_mean=apyr_mean,
        apyr_std=apyr_std,
        n_groups=n_groups,
    )


def build_model_pipeline() -> Pipeline:
    """Build sklearn pipeline with mixed numeric/categorical/text feature handling."""
    preprocess = ColumnTransformer(
        transformers=[
            ("num", StandardScaler(), _NUMERIC_COLUMNS),
            ("cat", OneHotEncoder(handle_unknown="ignore"), _CONDITION_COLUMNS),
            (
                "formed_vec",
                CountVectorizer(token_pattern=r"[^ ]+", binary=True),
                "formed_motifs_tokens",
            ),
            (
                "spectator_vec",
                CountVectorizer(token_pattern=r"[^ ]+", binary=True),
                "spectator_groups_tokens",
            ),
            (
                "sulf_motif_vec",
                CountVectorizer(token_pattern=r"[^ ]+", binary=True),
                "sulf_motif_tokens",
            ),
            (
                "sulf_fg_vec",
                CountVectorizer(token_pattern=r"[^ ]+", binary=True),
                "sulf_fg_tokens",
            ),
            (
                "bor_motif_vec",
                CountVectorizer(token_pattern=r"[^ ]+", binary=True),
                "bor_motif_tokens",
            ),
            (
                "bor_fg_vec",
                CountVectorizer(token_pattern=r"[^ ]+", binary=True),
                "bor_fg_tokens",
            ),
        ],
        remainder="drop",
        sparse_threshold=0.3,
    )

    model = RandomForestRegressor(
        n_estimators=400,
        min_samples_leaf=2,
        random_state=42,
        n_jobs=-1,
    )
    return Pipeline([("preprocess", preprocess), ("model", model)])


def run_random_split_cv(
    X: pd.DataFrame,
    y: np.ndarray,
    groups: Sequence[str],
    *,
    n_splits: int,
    test_size: float,
    random_state: int,
) -> Dict[str, Any]:
    """Evaluate repeated random splits."""
    splitter = ShuffleSplit(
        n_splits=n_splits,
        test_size=test_size,
        random_state=random_state,
    )
    fold_metrics: List[Dict[str, Any]] = []
    for fold_id, (train_idx, test_idx) in enumerate(splitter.split(X), start=1):
        pipeline = build_model_pipeline()
        pipeline.fit(X.iloc[train_idx], y[train_idx])
        pred = pipeline.predict(X.iloc[test_idx])
        result = evaluate_predictions(y[test_idx], pred, np.array(groups)[test_idx])
        fold_metrics.append(
            {
                "fold": fold_id,
                "spearman": result.spearman,
                "apyr_mean": result.apyr_mean,
                "apyr_std": result.apyr_std,
                "n_groups": result.n_groups,
            }
        )

    spearman_values = [m["spearman"] for m in fold_metrics]
    apyr_values = [m["apyr_mean"] for m in fold_metrics]
    return {
        "splits": n_splits,
        "test_size": test_size,
        "fold_metrics": fold_metrics,
        "spearman_mean": float(np.mean(spearman_values)),
        "spearman_std": float(np.std(spearman_values)),
        "apyr_mean": float(np.mean(apyr_values)),
        "apyr_std": float(np.std(apyr_values)),
    }


def run_logo_cv(
    X: pd.DataFrame,
    y: np.ndarray,
    groups: Sequence[str],
) -> Dict[str, Any]:
    """Leave-one-sulfonamide-out cross validation."""
    logo = LeaveOneGroupOut()
    fold_metrics: List[Dict[str, Any]] = []
    all_true: List[float] = []
    all_pred: List[float] = []
    all_groups: List[str] = []

    for fold_id, (train_idx, test_idx) in enumerate(logo.split(X, y, groups=groups), start=1):
        pipeline = build_model_pipeline()
        pipeline.fit(X.iloc[train_idx], y[train_idx])
        pred = pipeline.predict(X.iloc[test_idx])
        test_groups = np.array(groups)[test_idx]

        result = evaluate_predictions(y[test_idx], pred, test_groups)
        fold_metrics.append(
            {
                "fold": fold_id,
                "group": str(test_groups[0]) if len(test_groups) else "NA",
                "rows": int(len(test_idx)),
                "spearman": result.spearman,
                "apyr_mean": result.apyr_mean,
                "apyr_std": result.apyr_std,
            }
        )

        all_true.extend(y[test_idx].tolist())
        all_pred.extend(pred.tolist())
        all_groups.extend(test_groups.tolist())

    global_result = evaluate_predictions(
        np.array(all_true, dtype=float),
        np.array(all_pred, dtype=float),
        all_groups,
    )
    return {
        "folds": len(fold_metrics),
        "fold_metrics": fold_metrics,
        "global_spearman": global_result.spearman,
        "global_apyr_mean": global_result.apyr_mean,
        "global_apyr_std": global_result.apyr_std,
        "global_apyr_n_groups": global_result.n_groups,
    }


def fit_final_model(X: pd.DataFrame, y: np.ndarray) -> Pipeline:
    """Train final model on all rows."""
    pipeline = build_model_pipeline()
    pipeline.fit(X, y)
    return pipeline


def export_feature_importance(pipeline: Pipeline, output_csv: Path, top_n: int = 200) -> None:
    """Save top feature importances from random forest model."""
    model = pipeline.named_steps["model"]
    preprocess = pipeline.named_steps["preprocess"]
    feature_names = preprocess.get_feature_names_out()
    importances = model.feature_importances_
    imp_df = pd.DataFrame({"feature": feature_names, "importance": importances})
    imp_df = imp_df.sort_values("importance", ascending=False).head(top_n)
    imp_df.to_csv(output_csv, index=False)


def train_and_evaluate(
    input_csv: Path,
    output_dir: Path,
    *,
    n_random_splits: int,
    random_test_size: float,
    random_state: int,
) -> Dict[str, Any]:
    """Run the full descriptor + training + evaluation + export workflow."""
    output_dir.mkdir(parents=True, exist_ok=True)

    raw_df = load_chanlam_dataset(input_csv)
    feat_df = build_feature_table(raw_df)

    X = feat_df[_CONDITION_COLUMNS + _TEXT_COLUMNS + _NUMERIC_COLUMNS].copy()
    y = feat_df["yield"].to_numpy(dtype=float)
    groups = feat_df["sulfonamide_smiles"].fillna("NA").astype(str).to_numpy()

    random_cv = run_random_split_cv(
        X,
        y,
        groups,
        n_splits=n_random_splits,
        test_size=random_test_size,
        random_state=random_state,
    )
    logo_cv = run_logo_cv(X, y, groups)

    final_model = fit_final_model(X, y)
    model_path = output_dir / "chanlam_rf_model.joblib"
    dump(final_model, model_path)

    predictions = final_model.predict(X)
    pred_df = feat_df[
        ["yield", "sulfonamide_smiles", "boronic_smiles", "catalyst", "base", "solvent"]
    ].copy()
    pred_df["predicted_yield"] = predictions
    pred_df.to_csv(output_dir / "training_set_predictions.csv", index=False)

    export_feature_importance(final_model, output_dir / "feature_importance_top200.csv")

    metrics = {
        "input_csv": str(input_csv),
        "rows": int(len(feat_df)),
        "unique_sulfonamides": int(feat_df["sulfonamide_smiles"].nunique()),
        "unique_boronic_partners": int(feat_df["boronic_smiles"].nunique()),
        "unique_conditions": int(
            feat_df[["catalyst", "base", "solvent"]].astype(str).agg("|".join, axis=1).nunique()
        ),
        "random_split_cv": random_cv,
        "leave_one_sulfonamide_out_cv": logo_cv,
        "artifacts": {
            "model": str(model_path),
            "predictions": str(output_dir / "training_set_predictions.csv"),
            "feature_importance": str(output_dir / "feature_importance_top200.csv"),
        },
    }
    with (output_dir / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2)
    return metrics


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Train Chan-Lam condition-selection model with chemistry descriptors."
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=Path("data/HTE_db/literature/ChanLam_dataset_converted (2)_canonical.csv"),
        help="Path to converted Chan-Lam dataset CSV.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/ml/chanlam"),
        help="Directory for model and metrics artifacts.",
    )
    parser.add_argument(
        "--n-random-splits",
        type=int,
        default=10,
        help="Number of repeated random train/test splits.",
    )
    parser.add_argument(
        "--random-test-size",
        type=float,
        default=0.2,
        help="Random split test fraction.",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random seed for reproducibility.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    metrics = train_and_evaluate(
        input_csv=args.input_csv,
        output_dir=args.output_dir,
        n_random_splits=args.n_random_splits,
        random_test_size=args.random_test_size,
        random_state=args.random_state,
    )
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
