"""CLI for type-agnostic recommendation from converted generic records."""

from __future__ import annotations

import argparse
import json

from reactive_taxonomy import RxnMapperProvider

from .generic_api import recommend_generic_conditions
from .models import ChemistRankingPreferences
from .ranking_preferences import (
    RANKING_COMPONENTS,
    available_ranking_profiles,
    resolve_ranking_preferences,
)


def _ranking_preferences(args: argparse.Namespace) -> ChemistRankingPreferences:
    base = resolve_ranking_preferences(
        ChemistRankingPreferences(profile_id=args.ranking_profile)
    )
    if not args.ranking_weight:
        return base
    weights = dict(base.weights)
    for assignment in args.ranking_weight:
        try:
            name, raw_value = assignment.split("=", 1)
            value = float(raw_value)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid ranking weight {assignment!r}; use COMPONENT=VALUE"
            ) from exc
        if name not in RANKING_COMPONENTS:
            raise ValueError(
                f"Unknown ranking component {name!r}; choose from "
                + ", ".join(RANKING_COMPONENTS)
            )
        weights[name] = value
    return ChemistRankingPreferences(
        profile_id=f"{args.ranking_profile}:custom",
        weights=weights,
        customized=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recommend canonical condition recipes by reaction signature"
    )
    parser.add_argument("reaction_smiles")
    parser.add_argument(
        "--records",
        default="results/generic_conversion/records.jsonl",
        help="Canonical records/manifest or a persisted SQLite runtime index",
    )
    parser.add_argument("--top-k", type=int, default=5)
    parser.add_argument("--minimum-pool-size", type=int)
    parser.add_argument(
        "--use-rxnmapper",
        action="store_true",
        help=(
            "Use optional RXNMapper evidence for an unresolved/ambiguous query; "
            "results remain expert-review required"
        ),
    )
    parser.add_argument(
        "--unrestricted",
        action="store_true",
        help=(
            "Use the paired review-core index and relax fallback gates; expert "
            "review is mandatory"
        ),
    )
    profile_ids = tuple(
        value["profile_id"] for value in available_ranking_profiles()
    )
    parser.add_argument(
        "--ranking-profile",
        choices=profile_ids,
        default="default",
        help="Chemist-facing ranking priority preset",
    )
    parser.add_argument(
        "--ranking-weight",
        action="append",
        metavar="COMPONENT=VALUE",
        help=(
            "Override one preset weight; repeat for multiple components. "
            "Values are normalized and cannot change chemistry gates."
        ),
    )
    args = parser.parse_args()
    if args.use_rxnmapper and not RxnMapperProvider.is_available():
        parser.error(
            "RXNMapper is not installed; run "
            "'python -m pip install -r requirements-mapping.txt'"
        )
    try:
        preferences = _ranking_preferences(args)
    except ValueError as exc:
        parser.error(str(exc))
    result = recommend_generic_conditions(
        args.reaction_smiles,
        records_path=args.records,
        top_k=args.top_k,
        minimum_pool_size=args.minimum_pool_size,
        use_rxnmapper=args.use_rxnmapper,
        unrestricted_fallback=args.unrestricted,
        ranking_preferences=preferences,
    )
    print(json.dumps(result.to_dict(), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
