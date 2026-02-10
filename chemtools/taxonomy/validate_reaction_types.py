#!/usr/bin/env python3
"""
Validate reaction_types.v4.0.json reaction constraints schema.

This utility enforces a stable constraints schema across all reaction families.

Usage:
    python -m chemtools.taxonomy.validate_reaction_types
    python -m chemtools.taxonomy.validate_reaction_types --fix
    python -m chemtools.taxonomy.validate_reaction_types --check-only
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List

from chemtools.taxonomy.reaction_catalog import (
    REACTION_CONSTRAINT_KEYS,
    REACTION_TYPES_FILE,
    normalize_reaction_constraints,
)


def _load_payload(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _validate_constraints(payload: Dict[str, Any]) -> Dict[str, Any]:
    reactions = payload.get("reaction_types") or []
    issues: List[str] = []
    updates: List[str] = []
    allowed = set(REACTION_CONSTRAINT_KEYS)

    for idx, entry in enumerate(reactions):
        if not isinstance(entry, dict):
            issues.append(f"reaction[{idx}] is not an object")
            continue
        rid = str(entry.get("id") or f"<index:{idx}>")
        raw_constraints = entry.get("constraints")
        if raw_constraints is None:
            issues.append(f"{rid}: missing constraints")
        elif not isinstance(raw_constraints, dict):
            issues.append(f"{rid}: constraints must be an object")
            raw_constraints = {}

        raw_constraints = dict(raw_constraints or {})
        unknown = sorted(k for k in raw_constraints if k not in allowed)
        if unknown:
            issues.append(f"{rid}: unknown constraint keys: {', '.join(unknown)}")

        normalized = normalize_reaction_constraints(raw_constraints)
        # Keep only schema keys in normalized output.
        if raw_constraints != normalized:
            entry["constraints"] = normalized
            updates.append(rid)
        else:
            # Ensure explicit constraints field exists for all families.
            entry["constraints"] = normalized

    return {"issues": issues, "updates": updates, "count": len(reactions)}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate reaction taxonomy constraints schema.",
    )
    parser.add_argument(
        "--path",
        type=Path,
        default=REACTION_TYPES_FILE,
        help=f"Path to reaction types JSON (default: {REACTION_TYPES_FILE})",
    )
    parser.add_argument(
        "--fix",
        action="store_true",
        help="Write normalized constraints back to file.",
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Return non-zero if any schema issues are found.",
    )
    args = parser.parse_args()

    payload = _load_payload(args.path)
    result = _validate_constraints(payload)

    print(f"Checked {result['count']} reaction families.")
    print(f"Issues: {len(result['issues'])}")
    if result["issues"]:
        for issue in result["issues"][:30]:
            print(f"  - {issue}")
        if len(result["issues"]) > 30:
            print(f"  ... and {len(result['issues']) - 30} more")

    print(f"Normalized updates needed: {len(result['updates'])}")
    if args.fix and result["updates"]:
        with args.path.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
        print(f"Wrote normalized constraints to {args.path}")

    if args.check_only and result["issues"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
