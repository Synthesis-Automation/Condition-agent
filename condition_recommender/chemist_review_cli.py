"""CLI for blind chemist-review packet generation."""

from __future__ import annotations

import argparse
import json

from .chemist_review import generate_chemist_review_packet


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a blind generic-condition chemist review packet"
    )
    parser.add_argument("index_path")
    parser.add_argument("output_dir")
    parser.add_argument("--max-cases", type=int, default=300)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--top-k", type=int, default=3)
    parser.add_argument("--minimum-pool-size", type=int, default=1)
    args = parser.parse_args()
    report = generate_chemist_review_packet(
        args.index_path,
        args.output_dir,
        max_cases=args.max_cases,
        seed=args.seed,
        top_k=args.top_k,
        minimum_pool_size=args.minimum_pool_size,
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
