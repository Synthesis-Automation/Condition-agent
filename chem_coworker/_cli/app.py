"""Minimal command line entry point for condition recommendation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from condition_registry import ConditionConstraintSet, normalize_condition_constraint
from core_retrosynthesis import RouteSearchPolicy

from chem_coworker.assistance import (
    AssistanceController,
    AssistanceRequest,
    ConditionCapabilities,
    ConfirmedConstraint,
    MultistepCapabilities,
    OpenAICompatibleAssistanceTransport,
    RetrosynthesisCapabilities,
    render_assistance,
)

from chem_coworker.contracts import (
    CompletionChoice,
    ConditionRequest,
    ConditionReviewSettings,
    MultistepRetrosynthesisRequest,
)
from chem_coworker.multistep import MultistepRetrosynthesisCoworker
from chem_coworker.multistep_review import LLMMultistepReviewer
from chem_coworker.review import LLMConditionReviewer
from chem_coworker.retrosynthesis import RetrosynthesisCoworker
from chem_coworker.retrosynthesis_review import LLMRetrosynthesisReviewer
from chem_coworker.service import ConditionCoworker
from chem_coworker.contracts import RetrosynthesisRequest

from .config import load_config
from .interactive import InteractiveSettings, run_interactive
from .models import infer_provider, provider_model_set, selectable_models


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="chem-coworker",
        description="Condition recommendation and validated retrosynthesis",
    )
    parser.add_argument(
        "reaction_smiles",
        nargs="?",
        help="Reaction SMILES, or target SMILES in a retrosynthesis mode",
    )
    parser.add_argument(
        "--mode",
        choices=("conditions", "retro", "multistep"),
        default="conditions",
        help="Workflow for positional input and initial interactive mode",
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
    parser.add_argument(
        "--retrosynthesis-library",
        type=Path,
        help="Generic operator library for one-step retrosynthesis",
    )
    parser.add_argument("--max-realizations", type=int, default=3)
    parser.add_argument("--max-templates-to-apply", type=int, default=500)
    parser.add_argument("--max-candidates-to-validate", type=int, default=100)
    parser.add_argument(
        "--stock-index",
        type=Path,
        help="Supplier stock portfolio or literature molecule index",
    )
    parser.add_argument("--max-depth", type=int, choices=(2, 3), default=3)
    parser.add_argument("--per-step-top-k", type=int, default=5)
    parser.add_argument("--beam-width", type=int, default=20)
    parser.add_argument("--max-expansions", type=int, default=12)
    parser.add_argument(
        "--route-max-templates",
        type=int,
        default=40,
        help="Template-application budget for each multistep expansion",
    )
    parser.add_argument(
        "--route-validate-limit",
        type=int,
        default=10,
        help="Candidate-validation budget for each multistep expansion",
    )
    parser.add_argument(
        "--guidance",
        default="",
        help="Advisory route-ranking preferences for multistep LLM review",
    )
    parser.add_argument(
        "--retro-conditions",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Attach deterministic condition evidence to retro strategies",
    )
    parser.add_argument("--retro-condition-top-k", type=int, default=3)
    parser.add_argument(
        "--include-l0",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Allow the broadest validated operator fallback in retro mode",
    )
    parser.add_argument(
        "--retro-context",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use precedent context in deterministic retro ranking",
    )
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
    parser.add_argument(
        "--assistance",
        choices=("off", "shadow", "advisory"),
        default="off",
        help="Experimental common assistance controller (default: off)",
    )
    parser.add_argument(
        "--save-assistance-trace",
        type=Path,
        help="Opt in to saving the experimental assistance trace as JSON",
    )
    parser.add_argument(
        "--exclude-condition",
        action="append",
        default=[],
        help="Explicitly exclude a registry-resolved condition substance",
    )
    parser.add_argument("--maximum-temperature-c", type=float)
    parser.add_argument("--required-atmosphere")
    parser.add_argument("--required-solvent")
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


def _condition_constraints(args: argparse.Namespace) -> ConditionConstraintSet:
    raw = [
        ("excluded_substance", value) for value in args.exclude_condition
    ]
    raw.extend(
        (kind, value)
        for kind, value in (
            ("maximum_temperature_c", args.maximum_temperature_c),
            ("required_atmosphere", args.required_atmosphere),
            ("required_solvent", args.required_solvent),
        )
        if value is not None
    )
    constraints = []
    for kind, value in raw:
        resolution = normalize_condition_constraint(
            kind,  # type: ignore[arg-type]
            value,
            provenance="explicit_user",
        )
        if resolution.status != "resolved" or resolution.constraint is None:
            detail = ", ".join(resolution.warnings) or resolution.status
            raise ValueError(f"could not resolve {kind}={value!r}: {detail}")
        constraints.append(resolution.constraint)
    return ConditionConstraintSet(tuple(constraints))


def _run_assistance(
    args: argparse.Namespace,
    review_settings: ConditionReviewSettings,
    coworker: ConditionCoworker,
    retrosynthesis_coworker: RetrosynthesisCoworker | None,
    multistep_coworker: MultistepRetrosynthesisCoworker | None,
    condition_constraints: ConditionConstraintSet,
):
    controller = AssistanceController(
        transport=OpenAICompatibleAssistanceTransport(),
        condition_capabilities=ConditionCapabilities(coworker),
        retrosynthesis_capabilities=(
            RetrosynthesisCapabilities(retrosynthesis_coworker)
            if retrosynthesis_coworker is not None
            else None
        ),
        multistep_capabilities=(
            MultistepCapabilities(multistep_coworker)
            if multistep_coworker is not None
            else None
        ),
    )
    constraints = (
        ConfirmedConstraint(
            constraint_id="cli.top_k",
            owner="condition_recommender",
            kind="top_k",
            value=args.top_k,
            provenance="explicit_user",
        ),
        ConfirmedConstraint(
            constraint_id="cli.ranking_profile",
            owner="condition_recommender",
            kind="ranking_profile",
            value=args.ranking_profile,
            provenance="explicit_user",
        ),
    )
    route_policy = RouteSearchPolicy(
        initial_max_depth=args.max_depth,
        initial_beam_width=args.beam_width,
        initial_max_expansions=args.max_expansions,
        maximum_max_depth=args.max_depth,
        maximum_beam_width=args.beam_width,
        maximum_max_expansions=args.max_expansions,
    )
    request = AssistanceRequest(
        objective={
            "conditions": "Explain and compare deterministic condition recommendations",
            "retro": "Explain and compare forward-validated disconnection strategies",
            "multistep": "Explain and compare deterministic route-search results",
        }[args.mode],
        mode=args.mode,
        structure_input=args.reaction_smiles,
        constraints=constraints if args.mode == "conditions" else (),
        condition_constraints=condition_constraints,
        route_search_policy=route_policy,
        provider_settings={
            "provider": review_settings.provider,
            "model": review_settings.model,
            "reasoning_effort": review_settings.reasoning_effort,
            "max_output_tokens": review_settings.max_output_tokens,
        },
    )
    return controller, controller.run(request)


def main(argv: Sequence[str] | None = None) -> int:
    """Run one recommendation or an interactive session."""

    parser = _parser()
    args = parser.parse_args(argv)
    if args.mode == "multistep" and not 1 <= args.top_k <= 10:
        parser.error("--top-k must be between 1 and 10 in multistep mode")
    try:
        review_settings = _review_settings(args)
        condition_constraints = _condition_constraints(args)
    except ValueError as exc:
        parser.error(str(exc))
    reviewer = LLMConditionReviewer()
    retrosynthesis_reviewer = LLMRetrosynthesisReviewer()
    multistep_reviewer = LLMMultistepReviewer()
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
    retrosynthesis_coworker = None
    needs_retrosynthesis = bool(
        args.mode in {"retro", "multistep"}
        or args.reaction_smiles is None
        or args.interactive
        or args.retrosynthesis_library is not None
    )
    if needs_retrosynthesis:
        try:
            if args.retrosynthesis_library is None:
                retrosynthesis_coworker = RetrosynthesisCoworker.from_default(
                    condition_recommender=coworker.recommender,
                    reviewer=retrosynthesis_reviewer,
                )
            else:
                retrosynthesis_coworker = RetrosynthesisCoworker.from_path(
                    args.retrosynthesis_library,
                    condition_recommender=coworker.recommender,
                    reviewer=retrosynthesis_reviewer,
                )
        except (OSError, RuntimeError, ValueError) as exc:
            if args.mode in {"retro", "multistep"} or args.retrosynthesis_library is not None:
                parser.error(str(exc))
            print(f"Warning: retrosynthesis unavailable: {exc}")
    if retrosynthesis_coworker is not None:
        for warning in retrosynthesis_coworker.startup_warnings:
            print(f"Warning: {warning}")
    multistep_coworker = None
    if retrosynthesis_coworker is not None:
        try:
            multistep_coworker = (
                MultistepRetrosynthesisCoworker.from_retrosynthesis_coworker(
                    retrosynthesis_coworker,
                    stock_path=args.stock_index,
                    reviewer=multistep_reviewer,
                )
            )
        except (OSError, RuntimeError, ValueError) as exc:
            if args.mode == "multistep":
                parser.error(str(exc))
            print(f"Warning: multistep retrosynthesis unavailable: {exc}")
    if args.interactive or args.reaction_smiles is None:
        if args.assistance != "off":
            parser.error("experimental assistance currently requires one-shot input")
        if args.completion:
            parser.error("--completion is only supported in one-shot mode")
        return run_interactive(
            coworker,
            retrosynthesis_coworker=retrosynthesis_coworker,
            multistep_coworker=multistep_coworker,
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
                mode=args.mode,
                max_realizations_per_strategy=args.max_realizations,
                max_templates_to_apply=args.max_templates_to_apply,
                max_candidates_to_validate=args.max_candidates_to_validate,
                include_l0=args.include_l0,
                use_retro_context=args.retro_context,
                include_retro_conditions=args.retro_conditions,
                retro_condition_top_k=args.retro_condition_top_k,
                multistep_max_depth=args.max_depth,
                multistep_per_step_top_k=args.per_step_top_k,
                multistep_beam_width=args.beam_width,
                multistep_max_expansions=args.max_expansions,
                multistep_max_templates_to_apply=args.route_max_templates,
                multistep_max_candidates_to_validate=args.route_validate_limit,
                strategic_guidance=args.guidance,
            ),
            initial_reaction=args.reaction_smiles,
            enhanced=False if args.plain else None,
            persistent_history=not args.no_history,
        )
    if args.assistance != "off":
        if args.completion:
            parser.error("--completion is not yet supported with --assistance")
        assistance_controller, run = _run_assistance(
            args,
            review_settings,
            coworker,
            retrosynthesis_coworker,
            multistep_coworker,
            condition_constraints,
        )
        if (
            args.assistance == "advisory"
            and run.state.status == "needs_user_input"
            and run.state.unresolved_questions
        ):
            print(run.state.unresolved_questions[0].prompt)
            raw_value = input("> ").strip()
            if not raw_value:
                parser.error("a non-empty clarification answer is required")
            try:
                run = assistance_controller.resume_with_condition_constraint(
                    run.state,
                    raw_value,
                )
            except ValueError as exc:
                parser.error(str(exc))
        if args.save_assistance_trace is not None:
            args.save_assistance_trace.write_text(
                run.state.to_json() + "\n",
                encoding="utf-8",
            )
        if args.assistance == "advisory":
            print(
                json.dumps(run.to_dict(), ensure_ascii=False, indent=2)
                if args.as_json
                else render_assistance(run.state)
            )
            valid = getattr(run.authoritative_result, "valid", True)
            return 0 if valid and run.state.status == "completed" else 2
    if args.mode == "multistep":
        if args.completion:
            parser.error("--completion is not supported in multistep mode")
        if multistep_coworker is None:
            parser.error("multistep retrosynthesis is unavailable")
        response = multistep_coworker.plan(
            MultistepRetrosynthesisRequest(
                target_smiles=args.reaction_smiles,
                top_k=args.top_k,
                max_depth=args.max_depth,
                per_step_top_k=args.per_step_top_k,
                beam_width=args.beam_width,
                max_expansions=args.max_expansions,
                max_templates_to_apply=args.route_max_templates,
                max_candidates_to_validate=args.route_validate_limit,
                use_context=args.retro_context,
                include_l0=args.include_l0,
                include_conditions=args.retro_conditions,
                condition_top_k=args.retro_condition_top_k,
                strategic_guidance=args.guidance,
                review=review_settings,
            )
        )
        print(
            json.dumps(response.to_dict(), ensure_ascii=False, indent=2)
            if args.as_json
            else response.answer
        )
        return 0 if response.valid else 2
    if args.mode == "retro":
        if args.completion:
            parser.error("--completion is not supported in retro mode")
        if retrosynthesis_coworker is None:
            parser.error("retrosynthesis is unavailable")
        response = retrosynthesis_coworker.disconnect(
            RetrosynthesisRequest(
                target_smiles=args.reaction_smiles,
                top_k=args.top_k,
                max_realizations_per_strategy=args.max_realizations,
                max_templates_to_apply=args.max_templates_to_apply,
                max_candidates_to_validate=args.max_candidates_to_validate,
                use_context=args.retro_context,
                include_l0=args.include_l0,
                include_conditions=args.retro_conditions,
                condition_top_k=args.retro_condition_top_k,
                condition_minimum_pool_size=args.minimum_pool_size,
                unrestricted_condition_fallback=args.unrestricted_fallback,
                review=review_settings,
            )
        )
        print(
            json.dumps(response.to_dict(), ensure_ascii=False, indent=2)
            if args.as_json
            else response.answer
        )
        return 0 if response.valid else 2
    response = coworker.recommend(
        ConditionRequest(
            reaction_smiles=args.reaction_smiles,
            top_k=args.top_k,
            minimum_pool_size=args.minimum_pool_size,
            unrestricted_fallback=args.unrestricted_fallback,
            ranking_profile=args.ranking_profile,
            completion_choices=_completion_choices(args.completion),
            condition_constraints=condition_constraints,
            review=review_settings,
        )
    )
    print(
        json.dumps(response.to_dict(), ensure_ascii=False, indent=2)
        if args.as_json
        else response.answer
    )
    return 0 if response.valid else 2
