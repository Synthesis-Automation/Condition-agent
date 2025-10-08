#!/usr/bin/env python3
from __future__ import annotations

"""
Find reactions by condition core.

Examples:
  # Exact/substring match on Pd/XPhos across all families
  python scripts/find_by_core.py "Pd/XPhos" --limit 10 --pretty

  # Ligand-only fuzzy search (find any family containing 'XPhos')
  python scripts/find_by_core.py "XPhos" --pretty

  # Restrict to a family label (normalized like chemtools.precedent)
  python scripts/find_by_core.py "Pd/XPhos" --family "Suzuki_CC" --limit 5 --jsonl
"""

import argparse
import json
import os
import sys
from typing import List


# Ensure project root on path when running from source
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from chemtools.precedent import find_reactions_by_core  # type: ignore


def main(argv: List[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Find reactions by condition core (e.g., Pd/XPhos)")
    p.add_argument("core", help="Core query like 'Pd/XPhos', 'Pd', or 'XPhos'")
    p.add_argument("--family", help="Optional family filter (e.g., 'Ullmann C–N', 'Suzuki_CC')")
    p.add_argument("--no-fuzzy", dest="fuzzy", action="store_false", help="Disable fuzzy ligand name matching")
    p.set_defaults(fuzzy=True)
    p.add_argument("--limit", type=int, default=25, help="Max results (default: 25)")
    p.add_argument("--jsonl", action="store_true", help="Emit compact JSON lines")
    p.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    args = p.parse_args(argv)

    rows = find_reactions_by_core(args.core, family=args.family, fuzzy=bool(args.fuzzy), limit=int(args.limit))
    payload = {
        "query": {"core": args.core, "family": args.family, "fuzzy": bool(args.fuzzy)},
        "count": len(rows),
        "results": rows,
    }
    if args.jsonl and not args.pretty:
        print(json.dumps(payload, ensure_ascii=False))
    else:
        print(json.dumps(payload, ensure_ascii=False, indent=(2 if args.pretty else None)))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

