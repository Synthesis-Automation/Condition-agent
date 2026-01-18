"""
Evaluate the impact of spectator group weighting on HTE recommendations.

This script uses a train/test split of a single HTE dataset (default: C-S coupling)
and compares retrieval metrics with and without spectator group weighting.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
import tempfile
from pathlib import Path
from typing import Iterable, Optional, Tuple, Dict, Any, List

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from chemtools.HTE import HTERecommender


def _parse_reaction_smiles(reaction_smiles: str) -> Tuple[str, Optional[str], Optional[str]]:
    text = (reaction_smiles or "").strip()
    if not text:
        return "", None, None
    if ">>" in text:
        reactants_part, product = text.split(">>", 1)
        reactants = [r for r in reactants_part.split(".") if r]
        reactant_a = reactants[0] if reactants else ""
        reactant_b = ".".join(reactants[1:]) if len(reactants) > 1 else None
        product_smiles = product.strip() or None
        return reactant_a, reactant_b, product_smiles
    reactants = [r for r in text.split(".") if r]
    reactant_a = reactants[0] if reactants else ""
    reactant_b = ".".join(reactants[1:]) if len(reactants) > 1 else None
    return reactant_a, reactant_b, None


def _normalize_text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and pd.isna(value):
        return ""
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none"}:
        return ""
    return text


def _condition_key(
    catalyst: Any,
    ligand: Any,
    base: Any,
    solvent: Any,
) -> Tuple[str, str, str, str]:
    return (
        _normalize_text(catalyst),
        _normalize_text(ligand),
        _normalize_text(base),
        _normalize_text(solvent),
    )


def _row_condition_key(row: pd.Series) -> Tuple[str, str, str, str]:
    return _condition_key(
        row.get("catalyst"),
        row.get("ligand"),
        row.get("base"),
        row.get("solvent"),
    )


def _rec_condition_key(rec: Any) -> Tuple[str, str, str, str]:
    return _condition_key(rec.catalyst, rec.ligand, rec.base, rec.solvent)


def _rank_for_key(recs: Iterable[Any], target_key: Tuple[str, str, str, str]) -> Optional[int]:
    for idx, rec in enumerate(recs, start=1):
        if _rec_condition_key(rec) == target_key:
            return idx
    return None


def _compute_metrics(ranks: List[Optional[int]], ks: Iterable[int]) -> Dict[str, Any]:
    total = len(ranks)
    if total == 0:
        return {
            "count": 0,
            "coverage": 0.0,
            "mrr": 0.0,
            "avg_rank": None,
        }
    hits = [r for r in ranks if r is not None]
    coverage = len(hits) / total if total else 0.0
    mrr = sum(1.0 / r for r in hits) / len(hits) if hits else 0.0
    avg_rank = sum(hits) / len(hits) if hits else None
    metrics: Dict[str, Any] = {
        "count": total,
        "coverage": coverage,
        "mrr": mrr,
        "avg_rank": avg_rank,
    }
    for k in ks:
        metrics[f"hit@{k}"] = sum(1 for r in hits if r <= k) / total if total else 0.0
    return metrics


def _summarize_rank_deltas(
    ranks_with: List[Optional[int]],
    ranks_without: List[Optional[int]],
) -> Dict[str, Any]:
    improved = 0
    worsened = 0
    unchanged = 0
    gained = 0
    lost = 0
    deltas: List[int] = []
    for with_rank, without_rank in zip(ranks_with, ranks_without):
        if with_rank is None and without_rank is None:
            unchanged += 1
            continue
        if with_rank is None:
            lost += 1
            continue
        if without_rank is None:
            gained += 1
            continue
        delta = without_rank - with_rank
        deltas.append(delta)
        if delta > 0:
            improved += 1
        elif delta < 0:
            worsened += 1
        else:
            unchanged += 1
    avg_delta = sum(deltas) / len(deltas) if deltas else 0.0
    return {
        "improved": improved,
        "worsened": worsened,
        "unchanged": unchanged,
        "gained": gained,
        "lost": lost,
        "avg_rank_delta": avg_delta,
    }


def _save_train_csv(df: pd.DataFrame, output_path: Optional[Path]) -> Path:
    if output_path is None:
        handle = tempfile.NamedTemporaryFile(delete=False, suffix=".csv")
        output_path = Path(handle.name)
        handle.close()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    return output_path


def _prepare_split(
    df: pd.DataFrame,
    test_fraction: float,
    test_size: Optional[int],
    seed: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    indices = list(df.index)
    rng = random.Random(seed)
    rng.shuffle(indices)
    if test_size is not None:
        test_count = min(test_size, len(indices))
    else:
        test_count = max(1, int(len(indices) * test_fraction))
    test_set = set(indices[:test_count])
    test_df = df.loc[list(test_set)].copy()
    train_df = df.loc[[idx for idx in indices if idx not in test_set]].copy()
    return train_df, test_df


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Evaluate spectator group weighting for HTE recommendations."
    )
    parser.add_argument(
        "--input",
        default="data/HTE_db/literature/C-S Coupling_canonical.csv",
        help="Path to HTE CSV file with spectator_groups (default: C-S coupling).",
    )
    parser.add_argument(
        "--test-fraction",
        type=float,
        default=0.2,
        help="Fraction of data to reserve for testing (default: 0.2).",
    )
    parser.add_argument(
        "--test-size",
        type=int,
        default=None,
        help="Exact number of test rows (overrides --test-fraction).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=13,
        help="Random seed for split reproducibility.",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=10,
        help="Top-k recommendations to consider (default: 10).",
    )
    parser.add_argument(
        "--min-experiments",
        type=int,
        default=1,
        help="Minimum experiments per condition (default: 1).",
    )
    parser.add_argument(
        "--include-missing-spectator",
        action="store_true",
        help="Include rows with missing spectator_groups in the evaluation.",
    )
    parser.add_argument(
        "--train-output",
        default=None,
        help="Optional path to save the train CSV used for evaluation.",
    )
    parser.add_argument(
        "--output",
        default="results/spectator_groups_eval.json",
        help="Path to write JSON metrics (default: results/spectator_groups_eval.json).",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    df = pd.read_csv(input_path)
    df = df[df["reaction_smiles"].notna()]
    if not args.include_missing_spectator:
        df = df[df["spectator_groups"].apply(lambda v: bool(_normalize_text(v)))]

    train_df, test_df = _prepare_split(df, args.test_fraction, args.test_size, args.seed)
    train_path = _save_train_csv(train_df, Path(args.train_output) if args.train_output else None)

    train_condition_keys = {
        _row_condition_key(row) for _, row in train_df.iterrows()
    }

    recommender = HTERecommender(str(train_path))

    ranks_with: List[Optional[int]] = []
    ranks_without: List[Optional[int]] = []
    skipped_unseen = 0
    skipped_invalid = 0
    skipped_no_recs = 0

    for _, row in test_df.iterrows():
        key = _row_condition_key(row)
        if key not in train_condition_keys:
            skipped_unseen += 1
            continue
        reactant_a, reactant_b, product = _parse_reaction_smiles(
            _normalize_text(row.get("reaction_smiles"))
        )
        if not reactant_a:
            skipped_invalid += 1
            continue

        result_with = recommender.recommend(
            reactant_a_smiles=reactant_a,
            reactant_b_smiles=reactant_b,
            product_smiles=product,
            top_k=args.top_k,
            min_experiments=args.min_experiments,
            use_spectator_groups=True,
        )
        result_without = recommender.recommend(
            reactant_a_smiles=reactant_a,
            reactant_b_smiles=reactant_b,
            product_smiles=product,
            top_k=args.top_k,
            min_experiments=args.min_experiments,
            use_spectator_groups=False,
        )

        if not result_with.recommendations and not result_without.recommendations:
            skipped_no_recs += 1
            continue

        ranks_with.append(_rank_for_key(result_with.recommendations, key))
        ranks_without.append(_rank_for_key(result_without.recommendations, key))

    metrics_with = _compute_metrics(ranks_with, ks=[1, 3, 5, 10])
    metrics_without = _compute_metrics(ranks_without, ks=[1, 3, 5, 10])
    deltas = _summarize_rank_deltas(ranks_with, ranks_without)

    report = {
        "input": str(input_path),
        "train_rows": int(len(train_df)),
        "test_rows": int(len(test_df)),
        "evaluated_rows": int(len(ranks_with)),
        "skipped_unseen_conditions": int(skipped_unseen),
        "skipped_invalid_smiles": int(skipped_invalid),
        "skipped_no_recommendations": int(skipped_no_recs),
        "settings": {
            "top_k": args.top_k,
            "min_experiments": args.min_experiments,
            "test_fraction": args.test_fraction,
            "test_size": args.test_size,
            "seed": args.seed,
            "include_missing_spectator": args.include_missing_spectator,
        },
        "with_spectator_groups": metrics_with,
        "without_spectator_groups": metrics_without,
        "rank_deltas": deltas,
    }

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)

    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
