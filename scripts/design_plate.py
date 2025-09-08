from __future__ import annotations

import argparse
import os
import sys
from typing import List, Dict, Any


# Ensure project root on path when running from source
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from chemtools.recommend import design_plate_from_reaction


def main(argv: List[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Design a 24-well plate across cores for a reaction SMILES")
    p.add_argument("reaction", help="Reaction SMILES, e.g. 'Brc1ccccc1.NCCOC>>COCCNc1ccccc1'")
    p.add_argument("--plate-size", dest="plate_size", type=int, default=24, help="Number of wells (default: 24)")
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

    out = design_plate_from_reaction(args.reaction, plate_size=int(args.plate_size), relax=relax, constraint_rules=constraints)
    # Print CSV to stdout
    sys.stdout.write(out.get("csv") or "")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

