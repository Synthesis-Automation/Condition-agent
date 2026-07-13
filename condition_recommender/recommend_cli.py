"""CLI for reaction-label condition recommendation."""

from __future__ import annotations

import argparse
import json

from .label_api import recommend_conditions_from_labels


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recommend condition recipes from reactive-taxonomy labels"
    )
    parser.add_argument("reaction_smiles")
    parser.add_argument(
        "--records",
        default="datasets/reaction_label/v2.1_cleaned.csv",
        help="Cleaned reaction-label CSV",
    )
    parser.add_argument("--top-k", type=int, default=5)
    args = parser.parse_args()
    result = recommend_conditions_from_labels(
        args.reaction_smiles,
        records_path=args.records,
        top_k=args.top_k,
    )
    print(json.dumps(result.to_dict(), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
