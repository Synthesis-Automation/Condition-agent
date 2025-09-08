from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Dict, Any, List, Tuple

# Ensure project root on path when running from source
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from chemtools.smiles import normalize_reaction
from chemtools.router import detect_family
from chemtools.featurizers import ullmann as feat_ullmann
from chemtools.precedent import knn


def _pick_electrophile_nucleophile(reactants: List[str]) -> Tuple[str, str]:
    def is_electrophile(s: str) -> bool:
        t = (s or "").lower()
        return (
            ("br" in t) or ("cl" in t) or (" i" in t)
            or ("os(=o)(=o)c(f)(f)f" in t) or ("otf" in t)
        )
    if not reactants:
        return "", ""
    if len(reactants) == 1:
        return reactants[0], ""
    r0, r1 = reactants[0], reactants[1]
    if is_electrophile(r0):
        return r0, r1
    if is_electrophile(r1):
        return r1, r0
    return r0, r1


def main(argv: List[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Run precedent (kNN) search from a reaction SMILES")
    p.add_argument("reaction", help="Reaction SMILES, e.g. 'Brc1ccccc1.NCCOC>>COCCNc1ccccc1'")
    p.add_argument("--k", type=int, default=10, help="Number of neighbors (default: 10)")
    p.add_argument("--strict-bin", dest="strict_bin", action="store_true", help="Require exact bin matches before fallbacks")
    p.add_argument("--no-strict-bin", dest="strict_bin", action="store_false", help="Allow fallbacks when exact bin scarce")
    p.set_defaults(strict_bin=False)
    p.add_argument("--min-candidates", type=int, default=5, help="Minimum candidates to gather with fallbacks")
    p.add_argument("--jsonl", action="store_true", help="Emit a compact JSON object")
    p.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    args = p.parse_args(argv)

    # 1) Normalize reaction and extract reactant SMILES
    norm = normalize_reaction(args.reaction)
    reactants = [
        (r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or "")
        for r in (norm.get("reactants") or [])
    ]

    # 2) Detect family
    fam = detect_family(reactants).get("family") or "Unknown"

    # 3) Pick electrophile/nucleophile and featurize (Ullmann C–N supported)
    elec, nuc = _pick_electrophile_nucleophile(reactants)
    features: Dict[str, Any] = {}
    if fam == "Ullmann_CN":
        features = feat_ullmann.featurize(elec, nuc)
    else:
        # Best-effort: still use Ullmann featurizer heuristics; router may misclassify or multiple families may apply.
        features = feat_ullmann.featurize(elec, nuc)

    # 4) Run kNN
    relax = {"strict_bin": bool(args.strict_bin), "min_candidates": int(args.min_candidates)}
    result = knn(family=fam, features=features, k=int(args.k), relax=relax)

    payload = {
        "input_reaction": args.reaction,
        "family": fam,
        "electrophile": elec,
        "nucleophile": nuc,
        "features": features,
        "result": result,
    }
    if args.jsonl and not args.pretty:
        print(json.dumps(payload, ensure_ascii=False))
    else:
        print(json.dumps(payload, ensure_ascii=False, indent=(2 if args.pretty else None)))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

