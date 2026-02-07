"""Descriptor-model tuning and fixed external-holdout evaluation for Chan-Lam."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import argparse
import json
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from ml.chanlam_pipeline import (
    _CONDITION_COLUMNS,
    _NUMERIC_COLUMNS,
    _TEXT_COLUMNS,
    build_feature_table,
    load_chanlam_dataset,
)
from ml.evaluate_chanlam import run_logo_evaluation, run_random_split_evaluation


@dataclass(frozen=True)
class CandidateSpec:
    name: str
    model_type: str
    n_estimators: int
    min_samples_leaf: int
    max_features: float | str
    max_depth: Optional[int]
    include_substrate_ohe: bool


def _build_pipeline(
    *,
    categorical_cols: list[str],
    numeric_cols: list[str],
    text_cols: list[str],
    random_state: int,
    spec: CandidateSpec,
) -> Pipeline:
    transformers: list[tuple[str, Any, Any]] = []
    if numeric_cols:
        transformers.append(("num", StandardScaler(), numeric_cols))
    if categorical_cols:
        transformers.append(("cat", OneHotEncoder(handle_unknown="ignore"), categorical_cols))
    for col in text_cols:
        transformers.append((f"text_{col}", CountVectorizer(token_pattern=r"[^ ]+", binary=True), col))

    preprocess = ColumnTransformer(transformers=transformers, remainder="drop", sparse_threshold=0.3)
    model_kwargs = dict(
        n_estimators=spec.n_estimators,
        min_samples_leaf=spec.min_samples_leaf,
        max_features=spec.max_features,
        max_depth=spec.max_depth,
        random_state=random_state,
        n_jobs=-1,
    )
    if spec.model_type == "et":
        model = ExtraTreesRegressor(**model_kwargs)
    else:
        model = RandomForestRegressor(**model_kwargs)
    return Pipeline([("preprocess", preprocess), ("model", model)])


def _derive_external_holdout_groups(feat_df: pd.DataFrame) -> list[str]:
    """Fixed holdout list: substrates with partial condition coverage."""
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
    holdout = sorted(cond_counts[cond_counts < max_cov].index.tolist())
    return holdout


def _evaluate_external_holdout(
    X: pd.DataFrame,
    y: np.ndarray,
    groups: np.ndarray,
    holdout_groups: Sequence[str],
    *,
    pipeline_builder: Any,
    top_ks: Sequence[int],
    yield_thresholds: Sequence[float],
) -> Dict[str, Any]:
    from ml.evaluate_chanlam import evaluate_predictions, _topk_metrics_for_fold

    holdout_set = set(holdout_groups)
    test_mask = np.array([g in holdout_set for g in groups], dtype=bool)
    train_mask = ~test_mask
    if not np.any(test_mask):
        raise ValueError("No external holdout rows selected.")

    pipe = pipeline_builder()
    pipe.fit(X.iloc[train_mask], y[train_mask])
    pred = pipe.predict(X.iloc[test_mask])

    eval_global = evaluate_predictions(y[test_mask], pred, groups[test_mask])
    per_group_rows: List[Dict[str, Any]] = []
    test_df = pd.DataFrame(
        {
            "group": groups[test_mask],
            "y_true": y[test_mask],
            "y_pred": pred,
        }
    )
    for group, grp in test_df.groupby("group"):
        row = {
            "group": str(group),
            "rows": int(len(grp)),
        }
        local = evaluate_predictions(grp["y_true"].to_numpy(), grp["y_pred"].to_numpy(), [group] * len(grp))
        row["spearman"] = local.spearman
        row["apyr_mean"] = local.apyr_mean
        row["apyr_std"] = local.apyr_std
        row.update(
            _topk_metrics_for_fold(
                grp["y_true"].to_numpy(),
                grp["y_pred"].to_numpy(),
                top_ks=top_ks,
                yield_thresholds=yield_thresholds,
            )
        )
        per_group_rows.append(row)

    per_group_df = pd.DataFrame(per_group_rows)
    out = {
        "train_rows": int(np.sum(train_mask)),
        "test_rows": int(np.sum(test_mask)),
        "test_groups": int(len(holdout_groups)),
        "global_spearman": eval_global.spearman,
        "global_apyr_mean": eval_global.apyr_mean,
        "global_apyr_std": eval_global.apyr_std,
        "global_apyr_n_groups": eval_global.n_groups,
        "per_group": per_group_rows,
    }
    for k in top_ks:
        out[f"top{k}_hit_rate_ge_50"] = float(per_group_df[f"top{k}_hit_ge_50"].mean())
        out[f"top{k}_hit_rate_ge_70"] = float(per_group_df[f"top{k}_hit_ge_70"].mean())
        out[f"top{k}_avg_best_percentile"] = float(per_group_df[f"top{k}_best_percentile"].mean())
    return out


def _candidate_specs() -> list[CandidateSpec]:
    return [
        CandidateSpec("rf_desc_a", "rf", 250, 2, "sqrt", None, False),
        CandidateSpec("rf_desc_b", "rf", 400, 1, "sqrt", None, False),
        CandidateSpec("rf_desc_c", "rf", 400, 2, 0.7, None, False),
        CandidateSpec("et_desc_a", "et", 350, 1, "sqrt", None, False),
        CandidateSpec("rf_desc_plus_sub_a", "rf", 350, 2, "sqrt", None, True),
        CandidateSpec("rf_desc_plus_sub_b", "rf", 500, 1, "sqrt", None, True),
        CandidateSpec("et_desc_plus_sub_a", "et", 450, 1, "sqrt", None, True),
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
) -> Dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_df = load_chanlam_dataset(input_csv)
    feat_df = build_feature_table(raw_df)
    y = feat_df["yield"].to_numpy(dtype=float)
    groups = feat_df["sulfonamide_smiles"].fillna("NA").astype(str).to_numpy()

    baseline_categorical = list(_CONDITION_COLUMNS) + ["sulfonamide_smiles", "boronic_smiles"]
    baseline_X = feat_df[baseline_categorical].copy()

    def baseline_builder() -> Pipeline:
        spec = CandidateSpec(
            name="baseline_ohe",
            model_type="rf",
            n_estimators=350,
            min_samples_leaf=2,
            max_features="sqrt",
            max_depth=None,
            include_substrate_ohe=True,
        )
        return _build_pipeline(
            categorical_cols=baseline_categorical,
            numeric_cols=[],
            text_cols=[],
            random_state=random_state,
            spec=spec,
        )

    baseline_logo = run_logo_evaluation(
        baseline_X,
        y,
        groups,
        pipeline_builder=baseline_builder,
        top_ks=top_ks,
        yield_thresholds=yield_thresholds,
    )

    descriptor_base_cols = list(_CONDITION_COLUMNS) + list(_NUMERIC_COLUMNS) + list(_TEXT_COLUMNS)
    descriptor_plus_sub_cols = descriptor_base_cols + ["sulfonamide_smiles", "boronic_smiles"]

    candidate_reports: List[Dict[str, Any]] = []
    for spec in _candidate_specs():
        use_cols = descriptor_plus_sub_cols if spec.include_substrate_ohe else descriptor_base_cols
        X = feat_df[use_cols].copy()
        cat_cols = list(_CONDITION_COLUMNS)
        if spec.include_substrate_ohe:
            cat_cols = cat_cols + ["sulfonamide_smiles", "boronic_smiles"]

        def make_builder(_spec: CandidateSpec = spec, _cat_cols: list[str] = cat_cols) -> Pipeline:
            return _build_pipeline(
                categorical_cols=_cat_cols,
                numeric_cols=list(_NUMERIC_COLUMNS),
                text_cols=list(_TEXT_COLUMNS),
                random_state=random_state,
                spec=_spec,
            )

        random_eval = run_random_split_evaluation(
            X,
            y,
            groups,
            pipeline_builder=make_builder,
            n_splits=random_splits,
            test_size=random_test_size,
            random_state=random_state,
        )
        candidate_reports.append(
            {
                "name": spec.name,
                "spec": {
                    "model_type": spec.model_type,
                    "n_estimators": spec.n_estimators,
                    "min_samples_leaf": spec.min_samples_leaf,
                    "max_features": spec.max_features,
                    "max_depth": spec.max_depth,
                    "include_substrate_ohe": spec.include_substrate_ohe,
                },
                "random_split_cv": random_eval,
            }
        )

    # Pre-screen top candidates by random APYR, then run full LOGO on top 3.
    candidate_reports = sorted(
        candidate_reports,
        key=lambda r: (
            r["random_split_cv"]["apyr_mean"],
            r["random_split_cv"]["spearman_mean"],
        ),
        reverse=True,
    )
    finalists = candidate_reports[:3]

    final_with_logo: List[Dict[str, Any]] = []
    for entry in finalists:
        spec_dict = entry["spec"]
        spec = CandidateSpec(
            name=entry["name"],
            model_type=spec_dict["model_type"],
            n_estimators=spec_dict["n_estimators"],
            min_samples_leaf=spec_dict["min_samples_leaf"],
            max_features=spec_dict["max_features"],
            max_depth=spec_dict["max_depth"],
            include_substrate_ohe=spec_dict["include_substrate_ohe"],
        )
        use_cols = descriptor_plus_sub_cols if spec.include_substrate_ohe else descriptor_base_cols
        X = feat_df[use_cols].copy()
        cat_cols = list(_CONDITION_COLUMNS)
        if spec.include_substrate_ohe:
            cat_cols = cat_cols + ["sulfonamide_smiles", "boronic_smiles"]

        def make_builder(_spec: CandidateSpec = spec, _cat_cols: list[str] = cat_cols) -> Pipeline:
            return _build_pipeline(
                categorical_cols=_cat_cols,
                numeric_cols=list(_NUMERIC_COLUMNS),
                text_cols=list(_TEXT_COLUMNS),
                random_state=random_state,
                spec=_spec,
            )

        logo_eval = run_logo_evaluation(
            X,
            y,
            groups,
            pipeline_builder=make_builder,
            top_ks=top_ks,
            yield_thresholds=yield_thresholds,
        )
        row = dict(entry)
        row["logo_cv"] = logo_eval
        final_with_logo.append(row)

    best = sorted(
        final_with_logo,
        key=lambda r: (
            r["logo_cv"]["global_apyr_mean"],
            r["logo_cv"]["global_spearman"],
        ),
        reverse=True,
    )[0]

    # External holdout evaluation with fixed unseen substrate list.
    holdout_groups = _derive_external_holdout_groups(feat_df)
    holdout_path = output_dir / "external_holdout_sulfonamides.txt"
    holdout_path.write_text("\n".join(holdout_groups), encoding="utf-8")

    # baseline holdout
    baseline_holdout = _evaluate_external_holdout(
        baseline_X,
        y,
        groups,
        holdout_groups=holdout_groups,
        pipeline_builder=baseline_builder,
        top_ks=top_ks,
        yield_thresholds=yield_thresholds,
    )

    # best descriptor holdout
    best_spec_dict = best["spec"]
    best_spec = CandidateSpec(
        name=best["name"],
        model_type=best_spec_dict["model_type"],
        n_estimators=best_spec_dict["n_estimators"],
        min_samples_leaf=best_spec_dict["min_samples_leaf"],
        max_features=best_spec_dict["max_features"],
        max_depth=best_spec_dict["max_depth"],
        include_substrate_ohe=best_spec_dict["include_substrate_ohe"],
    )
    best_cols = descriptor_plus_sub_cols if best_spec.include_substrate_ohe else descriptor_base_cols
    best_X = feat_df[best_cols].copy()
    best_cat_cols = list(_CONDITION_COLUMNS)
    if best_spec.include_substrate_ohe:
        best_cat_cols = best_cat_cols + ["sulfonamide_smiles", "boronic_smiles"]

    def best_builder() -> Pipeline:
        return _build_pipeline(
            categorical_cols=best_cat_cols,
            numeric_cols=list(_NUMERIC_COLUMNS),
            text_cols=list(_TEXT_COLUMNS),
            random_state=random_state,
            spec=best_spec,
        )

    best_holdout = _evaluate_external_holdout(
        best_X,
        y,
        groups,
        holdout_groups=holdout_groups,
        pipeline_builder=best_builder,
        top_ks=top_ks,
        yield_thresholds=yield_thresholds,
    )

    report = {
        "input_csv": str(input_csv),
        "rows": int(len(feat_df)),
        "unique_sulfonamides": int(feat_df["sulfonamide_smiles"].nunique()),
        "external_holdout_groups_count": int(len(holdout_groups)),
        "external_holdout_list": holdout_groups,
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
        description="Tune descriptor model and run fixed external holdout evaluation.",
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
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    top_ks = _parse_int_list(args.top_ks)
    yield_thresholds = _parse_float_list(args.yield_thresholds)
    report = run_tuning_and_holdout(
        input_csv=args.input_csv,
        output_dir=args.output_dir,
        random_state=args.random_state,
        random_splits=args.random_splits,
        random_test_size=args.random_test_size,
        top_ks=top_ks,
        yield_thresholds=yield_thresholds,
    )
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

