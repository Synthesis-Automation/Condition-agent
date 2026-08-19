"""Minimal command line entry point for condition recommendation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from chem_coworker.contracts import CompletionChoice, ConditionRequest
from chem_coworker.service import ConditionCoworker

from .interactive import InteractiveSettings, run_interactive


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="chem-coworker",
        description="Structure-first reaction-condition recommendation",
    )
    parser.add_argument(
        "reaction_smiles",
        nargs="?",
        help="Reaction SMILES; omit it to start the interactive app",
    )
    parser.add_argument(
        "--index",
        type=Path,
        default=None,
        help=(
            "Recommendation index; by default the largest compatible full/compact "
            "artifact is selected"
        ),
    )
    parser.add_argument("--top-k", type=int, default=5)
    parser.add_argument("--minimum-pool-size", type=int)
    parser.add_argument("--ranking-profile", default="default")
    parser.add_argument("--unrestricted-fallback", action="store_true")
    parser.add_argument("--use-rxnmapper", action="store_true")
    parser.add_argument("--json", action="store_true", dest="as_json")
    parser.add_argument(
        "--interactive",
        action="store_true",
        help="Continue interactively after an optional initial reaction",
    )
    parser.add_argument(
        "--plain",
        action="store_true",
        help="Disable the enhanced terminal UI",
    )
    parser.add_argument(
        "--no-history",
        action="store_true",
        help="Do not persist reaction input history",
    )
    parser.add_argument(
        "--completion",
        action="append",
        default=[],
        metavar="REQUIREMENT=OPTION",
        help="Confirm an option; prefix a custom identifier with custom:",
    )
    return parser


def _completion_choices(values: Sequence[str]) -> tuple[CompletionChoice, ...]:
    choices = []
    for value in values:
        requirement, separator, selection = value.partition("=")
        if not separator or not requirement or not selection:
            raise ValueError("completion must use REQUIREMENT=OPTION")
        if selection.startswith("custom:"):
            choices.append(
                CompletionChoice(
                    requirement_id=requirement,
                    custom_identifier=selection.removeprefix("custom:"),
                )
            )
        else:
            choices.append(
                CompletionChoice(requirement_id=requirement, option_id=selection)
            )
    return tuple(choices)


def main(argv: Sequence[str] | None = None) -> int:
    """Run one recommendation or an interactive session."""

    parser = _parser()
    args = parser.parse_args(argv)
    try:
        if args.index is None:
            coworker = ConditionCoworker.from_default(
                use_rxnmapper=args.use_rxnmapper,
                include_review=args.unrestricted_fallback,
            )
        else:
            coworker = ConditionCoworker.from_path(
                args.index,
                use_rxnmapper=args.use_rxnmapper,
                include_review=args.unrestricted_fallback,
            )
    except (OSError, RuntimeError, ValueError) as exc:
        parser.error(str(exc))
    for warning in coworker.startup_warnings:
        print(f"Warning: {warning}")
    if args.interactive or args.reaction_smiles is None:
        if args.completion:
            parser.error("--completion is only supported in one-shot mode")
        return run_interactive(
            coworker,
            settings=InteractiveSettings(
                top_k=args.top_k,
                minimum_pool_size=args.minimum_pool_size,
                ranking_profile=args.ranking_profile,
                unrestricted_fallback=args.unrestricted_fallback,
                as_json=args.as_json,
            ),
            initial_reaction=args.reaction_smiles,
            enhanced=False if args.plain else None,
            persistent_history=not args.no_history,
        )
    response = coworker.recommend(
        ConditionRequest(
            reaction_smiles=args.reaction_smiles,
            top_k=args.top_k,
            minimum_pool_size=args.minimum_pool_size,
            unrestricted_fallback=args.unrestricted_fallback,
            ranking_profile=args.ranking_profile,
            completion_choices=_completion_choices(args.completion),
        )
    )
    print(
        json.dumps(response.to_dict(), ensure_ascii=False, indent=2)
        if args.as_json
        else response.answer
    )
    return 0 if response.valid else 2
