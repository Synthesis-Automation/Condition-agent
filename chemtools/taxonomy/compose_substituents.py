#!/usr/bin/env python3
"""Compose substituent groups from fragment definitions.

Usage:
    python -m chemtools.taxonomy.compose_substituents
    python -m chemtools.taxonomy.compose_substituents --write
    python -m chemtools.taxonomy.compose_substituents --write --inplace
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Set

from chemtools.taxonomy.substituent_composer import load_organic_groups_with_compositions


def _group_ids(groups: List[Dict[str, object]]) -> Set[str]:
    return {
        str(entry.get("id") or "").strip()
        for entry in groups
        if isinstance(entry, dict) and entry.get("id")
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compose substituent groups from taxonomy fragment rules.",
    )
    parser.add_argument(
        "--taxonomy-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "data",
        help="Path to taxonomy data directory.",
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="Write a materialized groups file with generated entries appended.",
    )
    parser.add_argument(
        "--inplace",
        action="store_true",
        help="When used with --write, overwrite organic_groups.v1.3.json.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output path for materialized groups JSON (ignored unless --write).",
    )
    args = parser.parse_args()

    taxonomy_dir = args.taxonomy_dir.resolve()
    groups_path = taxonomy_dir / "organic_groups.v1.3.json"
    if not groups_path.exists():
        print(f"ERROR: groups file not found: {groups_path}")
        return 1

    with groups_path.open("r", encoding="utf-8") as handle:
        base_payload = json.load(handle)
    base_groups = base_payload.get("groups", []) or []
    if not isinstance(base_groups, list):
        print("ERROR: groups payload malformed (groups must be a list).")
        return 1

    merged_payload = load_organic_groups_with_compositions(groups_path)
    merged_groups = merged_payload.get("groups", []) or []
    composed_meta = merged_payload.get("composed_groups", {}) or {}
    errors = composed_meta.get("errors", []) or []

    base_ids = _group_ids(base_groups)
    generated = [
        entry
        for entry in merged_groups
        if isinstance(entry, dict) and str(entry.get("id") or "").strip() not in base_ids
    ]

    print(f"Base groups: {len(base_groups)}")
    print(f"Generated groups: {len(generated)}")
    if generated:
        print("Generated IDs:")
        for entry in sorted(generated, key=lambda x: str(x.get("id") or "")):
            print(f"  - {entry.get('id')}")
    if errors:
        print("Composition errors:")
        for msg in errors:
            print(f"  - {msg}")

    if not args.write:
        return 0

    output_path = args.output
    if output_path is None:
        if args.inplace:
            output_path = groups_path
        else:
            output_path = groups_path.with_name("organic_groups.v1.3.composed.json")

    output_payload = dict(base_payload)
    output_payload["groups"] = list(base_groups) + generated
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(output_payload, handle, indent=2, ensure_ascii=False)

    print(f"Wrote materialized groups file: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
