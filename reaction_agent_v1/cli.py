"""CLI for reaction agent foundation."""

from __future__ import annotations

import argparse
import json
from typing import Any, Dict

from .gateway import ReactionAgentGateway


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run reaction agent foundation loop.")
    parser.add_argument("reaction", help="Reaction SMILES")
    parser.add_argument("--session-id", default=None, help="Optional session identifier")
    parser.add_argument("--max-steps", type=int, default=12, help="Maximum planner steps")
    parser.add_argument("--min-confidence", type=float, default=0.5, help="Validator minimum confidence")
    parser.add_argument("--json", action="store_true", help="Print full JSON output")
    return parser


def _print_human(payload: Dict[str, Any]) -> None:
    final_decision = payload["final_decision"]
    evidence = payload.get("evidence", {})
    diff_ready = bool((evidence.get("diff") or {}).get("principal_pair"))
    print(f"Session: {payload['session_id']}")
    print(f"Reaction: {payload['reaction_smiles']}")
    print(
        f"Final decision: {final_decision.get('reaction_type', 'unknown')} "
        f"(confidence={float(final_decision.get('confidence', 0.0)):.2f})"
    )
    print(f"Status: {payload['status']}")
    print(f"Evidence diff ready: {diff_ready}")
    print(f"Coverage suggestions: {len(payload['coverage_suggestions'])}")
    print(f"Tool artifacts: {sorted(payload.get('tool_artifacts', {}).keys())}")
    print("Trace:")
    for row in payload["trace"]:
        print(
            f"  - step={row['step_index']} action={row['action']} "
            f"status={row['status']} tool={row.get('tool_name')}"
        )


def main() -> int:
    args = _build_parser().parse_args()
    gateway = ReactionAgentGateway(min_confidence=args.min_confidence)
    result = gateway.run(
        args.reaction,
        session_id=args.session_id,
        max_steps=args.max_steps,
    )
    payload = result.to_dict()
    if args.json:
        print(json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True))
    else:
        _print_human(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
