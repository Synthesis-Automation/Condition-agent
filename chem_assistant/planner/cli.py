"""
Command-line runner for the auto-conditions pipeline.

Examples:
    python -m chem_assistant.planner.cli --reaction "Brc1ccccc1.N1CCOCC1>>Brc1ccccc1N1CCOCC1" --top-k 5
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict

from .api import ReactionInput, auto_conditions


def _serialize(obj: Any) -> Any:
    if hasattr(obj, "model_dump"):
        return obj.model_dump()
    if isinstance(obj, (list, tuple)):
        return [_serialize(x) for x in obj]
    if isinstance(obj, dict):
        return {k: _serialize(v) for k, v in obj.items()}
    return obj


def run_cli() -> int:
    parser = argparse.ArgumentParser(
        description="Run deterministic auto-conditions recommendations."
    )
    parser.add_argument(
        "--reaction",
        required=True,
        help="Reaction SMILES (reactants>>products).",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=5,
        help="Number of DRFP protocol candidates to retrieve.",
    )
    parser.add_argument(
        "--max-protocols",
        type=int,
        default=3,
        help="Maximum number of protocol outputs to format.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print full JSON payload instead of a summary.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Optional path to write JSON output.",
    )

    args = parser.parse_args()

    reaction = ReactionInput(reaction_smiles=args.reaction)
    result = auto_conditions(
        reaction,
        top_k_protocols=args.top_k,
        max_protocols=args.max_protocols,
        build_protocols=True,
    )

    payload: Dict[str, Any] = _serialize(result)

    if args.out:
        args.out.write_text(json.dumps(payload, indent=2))

    if args.json:
        print(json.dumps(payload, indent=2))
    else:
        print(f"Family: {result.family.family} (confidence={result.family.confidence})")
        print(f"Rule candidates: {len(result.rule_candidates)}")
        print(f"Protocol candidates: {len(result.protocol_candidates)}")
        if result.hte_summary:
            print(f"HTE type: {result.hte_summary.reaction_type}")
        print("\nTop protocols:")
        for proto in result.protocols:
            print(f"- {proto.candidate_id} ({proto.source})")
            for step in proto.additions[:6]:
                print(f"    • {step}")
            if len(proto.additions) > 6:
                print("    …")

    return 0


if __name__ == "__main__":
    raise SystemExit(run_cli())
