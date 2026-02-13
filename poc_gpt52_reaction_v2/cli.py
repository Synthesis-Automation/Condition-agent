"""CLI entrypoint for general-purpose reaction analysis PoC v2."""

from __future__ import annotations

import argparse
import json
from typing import Any, Dict

from .analyzer import analyze_reaction_general


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run taxonomy-first reaction analysis with optional gpt-5.2 reranking.",
    )
    parser.add_argument("reaction", help="Reaction SMILES (supports multi-component sides)")
    parser.add_argument("--use-llm", action="store_true", help="Enable optional LLM reranking.")
    parser.add_argument("--provider", default="openai", help="LLM provider (default: openai)")
    parser.add_argument("--model", default="gpt-5.2", help="LLM model (default: gpt-5.2)")
    parser.add_argument("--max-candidates", type=int, default=8, help="Max deterministic candidates.")
    parser.add_argument("--min-confidence", type=float, default=0.5, help="Minimum confidence gate.")
    parser.add_argument("--json", action="store_true", help="Print full JSON result.")
    return parser


def _print_human(payload: Dict[str, Any]) -> None:
    decision = payload["decision"]
    validation = payload["validation"]
    candidates = payload["taxonomy_candidates"]
    print(f"Reaction: {payload['reaction_smiles']}")
    print(
        f"Decision: {decision['reaction_type']} "
        f"(confidence={decision['confidence']:.2f}, source={decision['source']})"
    )
    print(f"Rationale: {decision['rationale']}")
    print(f"Validator: {'pass' if validation['passed'] else 'fail'}")
    if validation["issues"]:
        print("Issues:")
        for issue in validation["issues"]:
            print(f"  - {issue}")
    print(f"Candidates ({len(candidates)}):")
    for candidate in candidates[:5]:
        print(f"  - {candidate['reaction_type']} ({candidate['deterministic_score']:.2f})")


def main() -> int:
    args = _build_parser().parse_args()
    result = analyze_reaction_general(
        args.reaction,
        use_llm=args.use_llm,
        provider=args.provider,
        model=args.model,
        max_candidates=args.max_candidates,
        min_confidence=args.min_confidence,
    )
    payload = result.to_dict()

    if args.json:
        print(json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True))
    else:
        _print_human(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
