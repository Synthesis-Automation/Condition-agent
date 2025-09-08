from __future__ import annotations

import argparse
import json
import os
import sys
from typing import List, Dict, Any

# Ensure project root on path when running from source
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from chemtools.registry import search


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="List registry items filtered by role and optional query filters."
    )
    p.add_argument(
        "--role",
        choices=["CATALYST", "LIGAND", "BASE", "SOLVENT", "ADDITIVE"],
        help="Filter by role (case-insensitive)",
    )
    p.add_argument("--q", help="Optional substring filter against name/token/abbreviation/generic_core/uid")
    p.add_argument("--compound-type", dest="compound_type", help="Substring filter against compound_type")
    p.add_argument("--limit", type=int, default=100000, help="Maximum number of items to return (default: 100000)")
    p.add_argument("--sort-by", choices=["name", "uid"], default="name", help="Sort key (default: name)")
    p.add_argument("--jsonl", action="store_true", help="Emit JSON lines (one object per line)")
    p.add_argument("--pretty", action="store_true", help="Pretty-print a JSON array instead of lines")
    p.add_argument("--top", type=int, default=20, help="When not using JSON output, show the first N rows (default: 20)")

    args = p.parse_args(argv)

    items: List[Dict[str, Any]] = search(
        q=args.q,
        role=(args.role.upper() if args.role else None),
        compound_type=args.compound_type,
        limit=args.limit,
    )

    # Sort
    key = (args.sort_by or "name")
    items.sort(key=lambda r: str(r.get(key) or "").lower())

    # Output
    if args.jsonl:
        if args.pretty:
            print(json.dumps(items, ensure_ascii=False, indent=2))
        else:
            for obj in items:
                print(json.dumps(obj, ensure_ascii=False))
        return 0

    # Human-readable summary
    print(f"Found {len(items)} records" + (f" for role {args.role}" if args.role else ""))
    print("Showing top", min(args.top, len(items)))
    for i, rec in enumerate(items[: args.top], 1):
        uid = rec.get("uid", "")
        role = rec.get("role", "")
        name = rec.get("name", "")
        ctype = rec.get("compound_type", "")
        print(f"{i:3d}. {uid}\t{role}\t{name}\t{ctype}")

    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
