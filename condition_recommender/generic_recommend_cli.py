"""CLI for type-agnostic recommendation from converted generic records."""

from __future__ import annotations

import argparse
import json

from .generic_api import recommend_generic_conditions


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recommend canonical condition recipes by reaction signature"
    )
    parser.add_argument("reaction_smiles")
    parser.add_argument(
        "--records",
        default="results/generic_conversion/records.jsonl",
        help="Canonical records.jsonl or persisted generic index JSON",
    )
    parser.add_argument("--top-k", type=int, default=5)
    parser.add_argument("--minimum-pool-size", type=int)
    args = parser.parse_args()
    result = recommend_generic_conditions(
        args.reaction_smiles,
        records_path=args.records,
        top_k=args.top_k,
        minimum_pool_size=args.minimum_pool_size,
    )
    print(json.dumps(result.to_dict(), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
