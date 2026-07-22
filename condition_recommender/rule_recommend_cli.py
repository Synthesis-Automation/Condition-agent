"""CLI for structural expert-rule condition recommendations."""

from __future__ import annotations

import argparse
import json

from .rules import recommend_rule_conditions


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recommend expert condition templates from taxonomy facts"
    )
    parser.add_argument("reaction_smiles")
    parser.add_argument(
        "--include-draft",
        action="store_true",
        help="Include explicitly marked draft rules and templates for review",
    )
    args = parser.parse_args()
    result = recommend_rule_conditions(
        args.reaction_smiles,
        include_draft=args.include_draft,
    )
    print(json.dumps(result.to_dict(), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
