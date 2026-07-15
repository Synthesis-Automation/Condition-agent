"""CLI for Phase 2 reaction-edit machine and chemist evaluation."""

from __future__ import annotations

import argparse
import json

from .reaction_edit_evaluation import evaluate_reaction_edits


def main() -> int:
    """Run the Phase 2 benchmark and print its machine report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", help="Directory for evaluation artifacts")
    parser.add_argument("--benchmark", help="Optional benchmark manifest path")
    args = parser.parse_args()
    options = {"benchmark_path": args.benchmark} if args.benchmark else {}
    report = evaluate_reaction_edits(args.output_dir, **options)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["machine_gate_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
