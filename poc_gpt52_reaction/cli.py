"""CLI entrypoint for standalone GPT-5.2-style reaction analysis PoC."""

from __future__ import annotations

import argparse
import json
from typing import Any

from .analyzer import analyze_reaction


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run deterministic reaction analysis with optional gpt-5.2 refinement.",
    )
    parser.add_argument("reaction", help="Reaction SMILES (reactant>>product)")
    parser.add_argument("--use-llm", action="store_true", help="Enable optional LLM refinement.")
    parser.add_argument("--provider", default="openai", help="LLM provider (default: openai)")
    parser.add_argument("--model", default="gpt-5.2", help="LLM model (default: gpt-5.2)")
    parser.add_argument("--json", action="store_true", help="Print full JSON output.")
    return parser


def _print_human(result: dict[str, Any]) -> None:
    print(f"Reaction: {result['reaction_smiles']}")
    print(f"Top call: {result['best_hypothesis']['reaction_class']} (confidence={result['confidence']:.2f})")
    print(f"Summary: {result['summary']}")
    print("Reasoning steps:")
    for idx, step in enumerate(result["reasoning_steps"], start=1):
        print(f"  {idx}. {step['name']}: {step['conclusion']}")


def main() -> int:
    args = _build_parser().parse_args()
    result = analyze_reaction(
        args.reaction,
        use_llm=args.use_llm,
        provider=args.provider,
        model=args.model,
    ).to_dict()

    if args.json:
        print(json.dumps(result, ensure_ascii=True, indent=2, sort_keys=True))
    else:
        _print_human(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
