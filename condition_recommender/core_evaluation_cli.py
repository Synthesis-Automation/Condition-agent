"""CLI for leakage-safe reaction-core retrieval calibration."""

from __future__ import annotations

import argparse
import json

from .core_evaluation import evaluate_reaction_core_index


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Evaluate exact/local/context reaction-core retrieval"
    )
    parser.add_argument("records_path")
    parser.add_argument("output_dir")
    parser.add_argument("--test-fraction", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--top-k", type=int, default=5)
    parser.add_argument("--minimum-pool-size", type=int)
    parser.add_argument(
        "--split-mode",
        choices=(
            "grouped_random",
            "scaffold_disjoint",
            "source_disjoint",
            "forward_time",
        ),
        default="grouped_random",
    )
    args = parser.parse_args()
    report = evaluate_reaction_core_index(
        args.records_path,
        args.output_dir,
        test_fraction=args.test_fraction,
        seed=args.seed,
        top_k=args.top_k,
        minimum_pool_size=args.minimum_pool_size,
        split_mode=args.split_mode,
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
