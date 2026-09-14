"""Freeze, evaluate, and render a single-step retrosynthesis comparison."""

from __future__ import annotations

import argparse
import json

from core_retrosynthesis.strategy_evaluation import (
    evaluate_strategy_panel,
    freeze_strategy_panel,
)
from core_retrosynthesis.strategy_review import render_strategy_review


def main() -> None:
    """Run the selected evaluation stage without replacing frozen inputs."""

    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    freeze = commands.add_parser("freeze")
    freeze.add_argument("source")
    freeze.add_argument("output")
    freeze.add_argument("--size", type=int, default=1500)
    freeze.add_argument("--queries", type=int, default=60)
    freeze.add_argument("--seed", type=int, default=2026091401)
    freeze.add_argument("--exclude", action="append", default=[])
    evaluate = commands.add_parser("evaluate")
    evaluate.add_argument("panel")
    review = commands.add_parser("review")
    review.add_argument("report")
    args = parser.parse_args()
    if args.command == "freeze":
        print(freeze_strategy_panel(
            args.source, args.output, size=args.size, query_limit=args.queries,
            seed=args.seed, exclusion_paths=args.exclude,
        ))
    elif args.command == "evaluate":
        report = evaluate_strategy_panel(args.panel)
        print(json.dumps(report["summary"], indent=2))
    else:
        print(render_strategy_review(args.report))


if __name__ == "__main__":
    main()
