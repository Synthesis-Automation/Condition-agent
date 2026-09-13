"""CLI for canonical-reaction-group held-out recommendation evaluation."""

from __future__ import annotations

import argparse
import json

from .evaluation import evaluate_generic_index


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Evaluate generic recommendation without reaction leakage"
    )
    parser.add_argument(
        "records_path", help="Canonical records.jsonl or persisted generic index"
    )
    parser.add_argument("output_dir", help="Destination for evaluation artifacts")
    parser.add_argument("--test-fraction", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--top-k", type=int, default=5)
    parser.add_argument("--minimum-pool-size", type=int)
    parser.add_argument("--experimental-shared-core", action="store_true")
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
    parser.add_argument(
        "--retrieval-strategy",
        choices=(
            "hybrid",
            "family_only",
            "generic_only",
            "transformation_prior",
            "legacy_pilot",
        ),
        default="hybrid",
    )
    args = parser.parse_args()
    report = evaluate_generic_index(
        args.records_path,
        args.output_dir,
        test_fraction=args.test_fraction,
        seed=args.seed,
        top_k=args.top_k,
        minimum_pool_size=args.minimum_pool_size,
        split_mode=args.split_mode,
        retrieval_strategy=args.retrieval_strategy,
        experimental_shared_core=args.experimental_shared_core,
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
