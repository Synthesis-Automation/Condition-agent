from __future__ import annotations

import argparse
import json
import os
import sys

# Ensure project root on sys.path for local execution
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from chemtools import smiles, router
try:
    from chemtools.reaction_type_detector import detect_reaction_type as rxn_detect_type  # type: ignore
    HAS_RXN_INSIGHT = True
except Exception:
    HAS_RXN_INSIGHT = False


def main() -> int:
    p = argparse.ArgumentParser(description="Detect reaction type using rxn-insight (if available) with router fallback.")
    p.add_argument("reaction", help="Reaction SMILES: reactants>>products")
    p.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    args = p.parse_args()

    rxn = args.reaction
    norm = smiles.normalize_reaction(rxn)
    reactants = [
        (r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or "")
        for r in (norm.get("reactants") or [])
    ]
    fallback = router.detect_family(reactants)
    auto = None
    if HAS_RXN_INSIGHT:
        try:
            auto = rxn_detect_type(norm.get("normalized") or rxn)  # type: ignore[misc]
        except Exception:
            auto = None
    selected = None
    if isinstance(auto, dict) and (auto.get("mapped_family") or auto.get("success")):
        selected = auto.get("mapped_family") or fallback.get("family")
    else:
        selected = fallback.get("family")
    out = {
        "input": {"reaction_smiles": norm.get("normalized") or rxn},
        "rxn_insight_available": bool(HAS_RXN_INSIGHT),
        "rxn_insight": auto,
        "router_fallback": fallback,
        "selected_family": selected,
    }
    print(json.dumps(out, ensure_ascii=False, indent=2 if args.pretty else None))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

