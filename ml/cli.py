"""Unified CLI for rebuilt Chan-Lam ML system."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional, Sequence

from joblib import dump, load

from ml.contracts import BlendConfig, DatasetConfig, EvalConfig, FeatureConfig, RecommenderConfig, Stage1Config, Stage2Config
from ml.data import load_chanlam_feature_table
from ml.eval import run_logo_backtest
from ml.recommender import TwoStageRecommender
from ml.tuning import parse_float_list, parse_int_list, tune_from_prediction_csv


def _bool_from_flag(value: bool) -> bool:
    return not bool(value)


def _build_recommender_config(args: argparse.Namespace) -> RecommenderConfig:
    return RecommenderConfig(
        dataset=DatasetConfig(input_csv=Path(args.input_csv)),
        feature=FeatureConfig(
            profile=args.feature_profile,
            with_condition_props=_bool_from_flag(args.without_condition_props),
        ),
        stage1=Stage1Config(
            shortlist_k=int(args.shortlist_k),
            shrinkage_m=float(args.shrinkage_m),
        ),
        stage2=Stage2Config(
            model_type=args.model_type,
            n_estimators=int(args.n_estimators),
            random_state=int(args.random_state),
        ),
        blend=BlendConfig(
            w_stage1=float(args.w_stage1),
            w_stage2=float(args.w_stage2),
        ),
    )


def _add_shared_args(p: argparse.ArgumentParser) -> None:
    ds = DatasetConfig()
    feat = FeatureConfig()
    st2 = Stage2Config()
    st1 = Stage1Config()
    blend = BlendConfig()
    p.add_argument(
        "--input-csv",
        type=str,
        default=str(ds.input_csv),
    )
    p.add_argument("--feature-profile", type=str, default=feat.profile, choices=["core_full", "base_motif_spectator"])
    p.add_argument("--without-condition-props", action="store_true")
    p.add_argument("--model-type", type=str, default=st2.model_type, choices=["rf_reg", "lgbm_rank"])
    p.add_argument("--n-estimators", type=int, default=st2.n_estimators)
    p.add_argument("--random-state", type=int, default=st2.random_state)
    p.add_argument("--shortlist-k", type=int, default=st1.shortlist_k)
    p.add_argument("--shrinkage-m", type=float, default=st1.shrinkage_m)
    p.add_argument("--w-stage1", type=float, default=blend.w_stage1)
    p.add_argument("--w-stage2", type=float, default=blend.w_stage2)


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Rebuilt Chan-Lam ML system CLI.")
    sub = p.add_subparsers(dest="cmd", required=True)

    train_p = sub.add_parser("train", help="Train two-stage recommender and save artifact.")
    _add_shared_args(train_p)
    train_p.add_argument("--artifact", type=Path, default=Path("results/ml/chanlam_rebuild/recommender.joblib"))

    eval_p = sub.add_parser("evaluate", help="Run LOGO backtest with two-stage recommender.")
    _add_shared_args(eval_p)
    eval_p.add_argument("--output-json", type=Path, default=Path("results/ml/chanlam_rebuild/evaluation_report.json"))
    eval_p.add_argument("--top-ks", type=str, default="1,3,5")
    eval_p.add_argument("--thresholds", type=str, default="50,70")
    eval_p.add_argument("--max-folds", type=int, default=0)

    tune_p = sub.add_parser(
        "tune-defaults",
        help="Auto-tune shortlist and blend defaults from LOGO predictions.",
    )
    _add_shared_args(tune_p)
    tune_p.add_argument(
        "--eval-output-json",
        type=Path,
        default=Path("results/ml/chanlam_rebuild/evaluation_for_default_tuning.json"),
    )
    tune_p.add_argument(
        "--output-json",
        type=Path,
        default=Path("results/ml/chanlam_rebuild/default_tuning_report.json"),
    )
    tune_p.add_argument("--top-ks", type=str, default="1,3,5")
    tune_p.add_argument("--thresholds", type=str, default="50,70")
    tune_p.add_argument("--max-folds", type=int, default=0)
    tune_p.add_argument(
        "--shortlist-grid",
        type=str,
        default="5,10,15,20,25,30",
    )
    tune_p.add_argument(
        "--w-stage1-grid",
        type=str,
        default="0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0",
    )
    tune_p.add_argument(
        "--objective-mode",
        type=str,
        default="composite",
        choices=["composite", "apyr_lexicographic"],
    )
    tune_p.add_argument("--objective-w-apyr", type=float, default=0.50)
    tune_p.add_argument("--objective-w-top1", type=float, default=0.25)
    tune_p.add_argument("--objective-w-spearman", type=float, default=0.15)
    tune_p.add_argument("--objective-w-top3", type=float, default=0.05)
    tune_p.add_argument("--objective-w-stage2-pref", type=float, default=0.05)
    tune_p.add_argument(
        "--min-stage2-weight",
        type=float,
        default=0.10,
        help="Composite mode: enforce minimum stage2 blend weight to avoid stage1-only collapse.",
    )

    rec_p = sub.add_parser("recommend", help="Recommend conditions for a new reaction using saved artifact.")
    rec_p.add_argument("--artifact", type=Path, default=Path("results/ml/chanlam_rebuild/recommender.joblib"))
    rec_p.add_argument("--reaction-smiles", type=str, required=True)
    rec_p.add_argument("--top-k", type=int, default=10)
    rec_p.add_argument(
        "--output-csv",
        type=Path,
        default=Path("results/ml/chanlam_rebuild/new_reaction_top_conditions.csv"),
    )
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_parser().parse_args(argv)

    if args.cmd == "train":
        cfg = _build_recommender_config(args)
        df = load_chanlam_feature_table(cfg.dataset, with_condition_props=cfg.feature.with_condition_props)
        rec = TwoStageRecommender.create(cfg).fit(df)
        args.artifact.parent.mkdir(parents=True, exist_ok=True)
        dump(rec, args.artifact)
        summary = {
            "artifact": str(args.artifact),
            "rows": int(len(df)),
            "feature_profile": cfg.feature.profile,
            "with_condition_props": bool(cfg.feature.with_condition_props),
            "stage2_model_type": cfg.stage2.model_type,
            "condition_library_size": int(len(rec.condition_library)),
        }
        (args.artifact.with_suffix(".json")).write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print(json.dumps(summary, indent=2))
        return 0

    if args.cmd == "evaluate":
        cfg = _build_recommender_config(args)
        eval_cfg = EvalConfig(
            output_json=args.output_json,
            top_ks=tuple(parse_int_list(args.top_ks)),
            thresholds=tuple(parse_float_list(args.thresholds)),
            max_folds=int(args.max_folds),
        )
        report = run_logo_backtest(cfg, eval_cfg)
        print(json.dumps(report, indent=2))
        return 0

    if args.cmd == "tune-defaults":
        cfg = _build_recommender_config(args)
        eval_cfg = EvalConfig(
            output_json=args.eval_output_json,
            top_ks=tuple(parse_int_list(args.top_ks)),
            thresholds=tuple(parse_float_list(args.thresholds)),
            max_folds=int(args.max_folds),
        )
        eval_report = run_logo_backtest(cfg, eval_cfg)
        pred_csv = eval_cfg.output_json.with_suffix(".predictions.csv")
        tune_report = tune_from_prediction_csv(
            pred_csv,
            output_json=args.output_json,
            shortlist_grid=parse_int_list(args.shortlist_grid),
            w_stage1_grid=parse_float_list(args.w_stage1_grid),
            objective_mode=args.objective_mode,
            objective_weights={
                "apyr": float(args.objective_w_apyr),
                "top1": float(args.objective_w_top1),
                "spearman": float(args.objective_w_spearman),
                "top3": float(args.objective_w_top3),
                "stage2_pref": float(args.objective_w_stage2_pref),
            },
            min_stage2_weight=float(args.min_stage2_weight),
        )
        merged = {
            "evaluation_report": eval_report,
            "tuning_report": tune_report,
            "recommended_defaults": {
                "shortlist_k": tune_report["best"]["shortlist_k"],
                "w_stage1": tune_report["best"]["w_stage1"],
                "w_stage2": tune_report["best"]["w_stage2"],
            },
        }
        args.output_json.write_text(json.dumps(merged, indent=2), encoding="utf-8")
        print(json.dumps(merged, indent=2))
        return 0

    if args.cmd == "recommend":
        rec: TwoStageRecommender = load(args.artifact)
        out = rec.recommend(args.reaction_smiles, top_k=int(args.top_k))
        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.output_csv, index=False)
        print(out.to_string(index=False))
        return 0

    raise ValueError(f"Unknown command: {args.cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
