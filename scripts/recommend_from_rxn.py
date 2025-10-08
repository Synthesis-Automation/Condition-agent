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

from chemtools.recommend import recommend_from_reaction


def main(argv: List[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Recommend conditions from a reaction SMILES")
    p.add_argument("reaction", help="Reaction SMILES, e.g. 'Brc1ccccc1.NCCOC>>COCCNc1ccccc1'")
    p.add_argument("--k", type=int, default=25, help="Neighbors for prototype (default: 25)")
    p.add_argument("--jsonl", action="store_true", help="Emit compact JSON (one line)")
    p.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    # Constraint knobs (simple flags)
    p.add_argument("--no-hmpa", dest="no_hmpa", action="store_true", help="Disallow HMPA")
    p.add_argument("--no-chloro", dest="no_chloro", action="store_true", help="Disallow chlorinated solvents")
    p.add_argument("--aqueous-only", dest="aqueous_only", action="store_true", help="Allow only water as solvent")
    p.add_argument("--min-bp", dest="min_bp", type=float, default=None, help="Minimum solvent boiling point (C)")
    args = p.parse_args(argv)

    relax: Dict[str, Any] = {
        "use_drfp": True,
        "precompute_drfp": True,
        "precompute_scope": "candidates",
    }
    constraints: Dict[str, Any] = {
        "no_HMPA": bool(args.no_hmpa),
        "no_chlorinated": bool(args.no_chloro),
        "aqueous_only": bool(args.aqueous_only),
    }
    if args.min_bp is not None:
        constraints["min_bp_C"] = float(args.min_bp)

    out = recommend_from_reaction(args.reaction, k=int(args.k), relax=relax, constraint_rules=constraints)
    if args.jsonl and not args.pretty:
        print(json.dumps(out, ensure_ascii=False))
    else:
        print(json.dumps(out, ensure_ascii=False, indent=(2 if args.pretty else None)))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

