from __future__ import annotations

"""Command line interface for the SchemeConditionDB matcher."""

import argparse
import json
import sys
from pathlib import Path
from typing import Sequence

from .loader import SchemeConditionDBError, load_db
from .matcher import match
from .types import MatchResult, SchemeConditionDB, SelectorRuleDB


def _print_text(result: MatchResult) -> None:
    """Emit a minimal human-readable summary."""

    print(f"Reaction type: {result.reaction_type}")
    print(f"Match type: {result.match_type}")
    entry_label = result.entry_name or "(unnamed)"
    print(f"Entry: {result.entry_id} {entry_label}")
    print("Conditions:")
    for key, value in result.conditions.items():
        if isinstance(value, dict):
            serialized = json.dumps(value, indent=2, sort_keys=True)
            print(f"  {key}: {serialized}")
        else:
            print(f"  {key}: {value}")


def _run_match(args: argparse.Namespace) -> int:
    try:
        db = load_db(args.db)
        if isinstance(db, SchemeConditionDB):
            if args.rxn is None:
                raise ValueError("--rxn is required when matching SchemeConditionDB payloads")
            result = match(db, args.rxn)
        elif isinstance(db, SelectorRuleDB):
            if args.features is None:
                raise ValueError("--features is required when matching selector rule payloads")
            features_path: Path = args.features
            if not features_path.exists():
                raise ValueError(f"Features file not found: {features_path}")
            features_text = features_path.read_text(encoding="utf-8")
            try:
                features_payload = json.loads(features_text)
            except json.JSONDecodeError as exc:  # pragma: no cover - defensive
                raise ValueError("Failed to parse features JSON") from exc
            if not isinstance(features_payload, dict):
                raise ValueError("Features JSON must be an object at the top level")
            result = match(db, features=features_payload)
        else:  # pragma: no cover - defensive guard
            raise TypeError(f"Unsupported database type: {type(db)!r}")
    except (SchemeConditionDBError, ValueError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    payload = result.to_json_dict()
    if args.no_trace:
        payload.pop("trace", None)

    if args.json:
        indent = args.indent if args.indent is not None else 2
        print(json.dumps(payload, indent=indent, sort_keys=True))
    else:
        _print_text(result)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="scdb", description="Rule database matcher")
    subparsers = parser.add_subparsers(dest="command", required=True)

    match_parser = subparsers.add_parser("match", help="Match reaction conditions for a rule database")
    match_parser.add_argument("--db", required=True, type=Path, help="Path to a rule database JSON file")
    match_parser.add_argument("--rxn", help="Reaction SMILES (required for scheme databases)")
    match_parser.add_argument(
        "--features",
        type=Path,
        help="Path to JSON feature payload (required for selector databases)",
    )
    match_parser.add_argument("--json", action="store_true", help="Emit JSON instead of text")
    match_parser.add_argument("--no-trace", action="store_true", help="Exclude trace information from JSON output")
    match_parser.add_argument("--indent", type=int, default=None, help="Indent level for JSON output")
    match_parser.set_defaults(func=_run_match)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
