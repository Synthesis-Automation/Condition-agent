"""CLI for the standalone Suzuki recommendation baseline."""

from __future__ import annotations

import argparse
import json

from .api import recommend_conditions


def main() -> None:
    parser = argparse.ArgumentParser(description="Recommend Suzuki condition recipes")
    parser.add_argument("reaction_smiles")
    parser.add_argument("--records", default="results/suzuki_recommendation_pilot/verified.csv")
    parser.add_argument("--top-k", type=int, default=5)
    args = parser.parse_args()
    result = recommend_conditions(args.reaction_smiles, records_path=args.records, top_k=args.top_k)
    print(json.dumps(result.to_dict(), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
