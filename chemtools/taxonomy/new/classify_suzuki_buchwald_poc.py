#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
POC classifier for Suzuki + Buchwald using:
- atomic SMARTS features
- derived boolean features
- simple reaction-type rules

Requires RDKit.
Usage:
  python classify_suzuki_buchwald_poc.py --reactant "Brc1ccccc1" --reactant "OB(O)c1ccccc1"
"""
import argparse
import json
import sys
from pathlib import Path

def _ensure_repo_root_on_path() -> None:
    try:
        __import__("chemtools")
        return
    except ModuleNotFoundError:
        repo_root = Path(__file__).resolve().parents[3]
        sys.path.insert(0, str(repo_root))

def mixture_or(feature_dicts):
    out = {}
    for fd in feature_dicts:
        for k, v in fd.items():
            out[k] = out.get(k, False) or bool(v)
    return out

def main():
    _ensure_repo_root_on_path()

    ap = argparse.ArgumentParser()
    ap.add_argument("--reactant", action="append", required=True, help="Reactant SMILES (repeatable)")
    ap.add_argument("--atomic", default="calculable_features.atomic.suzuki_buchwald.poc.json")
    ap.add_argument("--derived", default="calculable_features.derived.suzuki_buchwald.poc.json")
    ap.add_argument("--reaction-types", default="reaction_types.suzuki_buchwald.poc.json")
    args = ap.parse_args()

    from chemtools.taxonomy.new.feature_engine import (
        classify_reaction_smiles,
        load_feature_definitions,
        validate_features,
    )

    folder = Path(__file__).resolve().parent
    atomic_path = (folder / args.atomic) if not Path(args.atomic).is_absolute() else Path(args.atomic)
    derived_path = (folder / args.derived) if not Path(args.derived).is_absolute() else Path(args.derived)
    reaction_types_path = (
        (folder / args.reaction_types)
        if not Path(args.reaction_types).is_absolute()
        else Path(args.reaction_types)
    )

    atomic_defs, derived_defs = load_feature_definitions(atomic_path=atomic_path, derived_path=derived_path)
    reaction_types = json.loads(reaction_types_path.read_text(encoding="utf-8"))
    validate_features(atomic_defs, derived_defs, reaction_types=reaction_types)

    result = classify_reaction_smiles(
        args.reactant,
        atomic_defs=atomic_defs,
        derived_defs=derived_defs,
        reaction_types=reaction_types,
    )

    hits = result.get("reaction_type_matches") or []
    print("Matched reaction types:", [m.id for m in hits])
    # Print a few key tokens
    keys = ["sp2_electrophile_present","suzuki_boron_partner_present","buchwald_amine_partner_present"]
    mix = result.get("reaction_has") or {}
    print("Key mixture tokens:", {k: mix.get(k, False) for k in keys})

if __name__ == "__main__":
    main()
