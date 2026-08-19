"""Minimal command line entry point for condition recommendation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from chem_coworker.contracts import (
    CompletionChoice,
    ConditionRequest,
    ConditionReviewSettings,
)
from chem_coworker.review import LLMConditionReviewer
from chem_coworker.service import ConditionCoworker

from .config import load_config
from .interactive import InteractiveSettings, run_interactive
from .models import infer_provider, provider_model_set, selectable_models


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
    parser.add_argument(
        "--review",
        choices=("off", "auto", "always"),
        default=None,
        help="LLM review mode (saved default: auto)",
    )
    parser.add_argument("--model", help="LLM review model")
    parser.add_argument("--provider", choices=("openai", "aliyun"))
    parser.add_argument(
        "--reasoning-effort",
        choices=("none", "low", "medium", "high", "xhigh", "max"),
    )
    parser.add_argument("--review-candidates", type=int)
    parser.add_argument("--review-max-tokens", type=int)
    parser.add_argument(
        "--review-order",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Apply the advisory LLM verdict order to presentation",
    )
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


def _review_settings(args: argparse.Namespace) -> ConditionReviewSettings:
    saved = load_config()
    model = args.model or str(saved["name"])
    known_models = {item["name"] for item in selectable_models()}
    if args.model and model not in known_models:
        raise ValueError(f"Unknown review model: {model}")
    if args.provider:
        provider = args.provider
    elif args.model:
        provider = infer_provider(model)
    else:
        provider = str(saved["provider"])
    if args.model and model not in provider_model_set(provider):
        raise ValueError(f"Review model {model} does not belong to provider {provider}")
    if model not in provider_model_set(provider):
        model = next(
            item["name"] for item in selectable_models() if item["provider"] == provider
        )
    return ConditionReviewSettings(
        mode=args.review or str(saved["review_mode"]),
        provider=provider,
        model=model,
        reasoning_effort=(args.reasoning_effort or str(saved["reasoning_effort"])),
        max_candidates=(
            args.review_candidates
            if args.review_candidates is not None
            else int(saved["review_candidates"])
        ),
        max_output_tokens=(
            args.review_max_tokens
            if args.review_max_tokens is not None
            else int(saved["review_max_tokens"])
        ),
        apply_order=(
            args.review_order
            if args.review_order is not None
            else bool(saved["apply_review_order"])
        ),
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Run one recommendation or an interactive session."""

    parser = _parser()
    args = parser.parse_args(argv)
    try:
        review_settings = _review_settings(args)
    except ValueError as exc:
        parser.error(str(exc))
    reviewer = LLMConditionReviewer()
    try:
        if args.index is None:
            coworker = ConditionCoworker.from_default(
                use_rxnmapper=args.use_rxnmapper,
                include_review=args.unrestricted_fallback,
                reviewer=reviewer,
            )
        else:
            coworker = ConditionCoworker.from_path(
                args.index,
                use_rxnmapper=args.use_rxnmapper,
                include_review=args.unrestricted_fallback,
                reviewer=reviewer,
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
                review_mode=review_settings.mode,
                provider=review_settings.provider,
                model=review_settings.model,
                reasoning_effort=review_settings.reasoning_effort,
                review_candidates=review_settings.max_candidates,
                review_max_tokens=review_settings.max_output_tokens,
                apply_review_order=review_settings.apply_order,
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
            review=review_settings,
        )
    )
    print(
        json.dumps(response.to_dict(), ensure_ascii=False, indent=2)
        if args.as_json
        else response.answer
    )
    return 0 if response.valid else 2
