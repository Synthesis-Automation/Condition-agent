"""CLI for Phase 1 molecular-feature machine and chemist evaluation."""

from __future__ import annotations

import argparse
import json

from .molecular_feature_evaluation import evaluate_molecular_features


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Evaluate molecular reactivity features and build review artifacts"
    )
    parser.add_argument("output_dir")
    parser.add_argument("--benchmark")
    args = parser.parse_args()
    options = {"benchmark_path": args.benchmark} if args.benchmark else {}
    report = evaluate_molecular_features(args.output_dir, **options)
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
