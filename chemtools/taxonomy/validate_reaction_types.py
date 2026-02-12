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

from chemtools.taxonomy import loader as taxonomy_loader
from chemtools.taxonomy.reaction_catalog import (
    REACTION_CONSTRAINT_KEYS,
    REACTION_SYNTHON_SLOT_KEYS,
    REACTION_TYPES_FILE,
    normalize_reaction_constraints,
    normalize_reaction_synthons,
)


def _load_payload(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _validate_constraints(payload: Dict[str, Any]) -> Dict[str, Any]:
    reactions = payload.get("reaction_types") or []
    issues: List[str] = []
    updates: List[str] = []
    allowed = set(REACTION_CONSTRAINT_KEYS)
    allowed_synthon_slot_keys = set(REACTION_SYNTHON_SLOT_KEYS)
    known_synthon_ids = {
        str(entry.get("id")).strip()
        for entry in (taxonomy_loader.load_synthons().get("synthons") or [])
        if isinstance(entry, dict) and str(entry.get("id") or "").strip()
    }

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

        raw_synthons = entry.get("synthons")
        if raw_synthons is not None and not isinstance(raw_synthons, dict):
            issues.append(f"{rid}: synthons must be an object")
            raw_synthons = {}
        raw_synthons = dict(raw_synthons or {})
        for slot_name, slot_spec in raw_synthons.items():
            if not isinstance(slot_spec, dict):
                continue
            unknown_slot_keys = sorted(
                key for key in slot_spec.keys()
                if key not in allowed_synthon_slot_keys
            )
            if unknown_slot_keys:
                issues.append(
                    f"{rid}: synthons.{slot_name} unknown keys: {', '.join(unknown_slot_keys)}"
                )
        normalized_synthons = normalize_reaction_synthons(raw_synthons)
        unknown_synthon_ids = sorted(
            synthon_id
            for slot_req in normalized_synthons.values()
            for synthon_id in slot_req.allowed
            if known_synthon_ids and synthon_id not in known_synthon_ids
        )
        if unknown_synthon_ids:
            issues.append(
                f"{rid}: unknown synthon ids: {', '.join(_dedupe(unknown_synthon_ids))}"
            )
        if raw_synthons != {
            slot: {
                "include": value.allowed,
                "min_hits": value.min_hits,
                "min_reactants": value.min_reactants,
            }
            for slot, value in normalized_synthons.items()
        }:
            if raw_synthons:
                updates.append(f"{rid}:synthons")
        if raw_synthons:
            # Keep source-like structure but normalized values.
            entry["synthons"] = {
                slot: {
                    "include": value.allowed,
                    "min_hits": value.min_hits,
                    "min_reactants": value.min_reactants,
                }
                for slot, value in normalized_synthons.items()
            }

    return {"issues": issues, "updates": updates, "count": len(reactions)}


def _dedupe(values: List[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


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
