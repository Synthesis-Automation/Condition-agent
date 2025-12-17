#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CLI demo for Suzuki + Buchwald POC classification.

Examples (from repo root):
  python -m chemtools.taxonomy.new.rxn_classifier --reactant "Brc1ccccc1" --reactant "OB(O)c1ccccc1"

Example (from this folder):
  python -m rxn_classifier --reactant "Brc1ccccc1" --reactant "OB(O)c1ccccc1"
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any


def _ensure_repo_root_on_path() -> None:
    """
    Allow running `python -m rxn_classifier` from this directory.
    """
    try:
        __import__("chemtools")
        return
    except ModuleNotFoundError:
        repo_root = Path(__file__).resolve().parents[3]
        sys.path.insert(0, str(repo_root))


def main(argv: list[str] | None = None) -> int:
    _ensure_repo_root_on_path()

    from chemtools.taxonomy.new.feature_engine import (
        classify_reaction_smiles,
        load_feature_definitions,
        validate_features,
    )
    from chemtools.util.rdkit_helpers import rdkit_available

    ap = argparse.ArgumentParser()
    ap.add_argument("--reactant", action="append", required=True, help="Reactant SMILES (repeatable)")
    ap.add_argument("--atomic", default="calculable_features.atomic.suzuki_buchwald.poc.json")
    ap.add_argument("--derived", default="calculable_features.derived.suzuki_buchwald.poc.json")
    ap.add_argument("--reactant-types", default="reactant_types.suzuki_buchwald.poc.json")
    ap.add_argument("--reaction-types", default="reaction_types.suzuki_buchwald.poc.json")
    ap.add_argument("--no-validate", action="store_true", help="Skip SMARTS + schema validation")
    ap.add_argument("--json", action="store_true", help="Output full result JSON")
    args = ap.parse_args(argv)

    if not rdkit_available():
        print("RDKit is not available (or disabled via CHEMTOOLS_DISABLE_RDKIT=1).", file=sys.stderr)
        return 2

    folder = Path(__file__).resolve().parent
    atomic_path = (folder / args.atomic) if not Path(args.atomic).is_absolute() else Path(args.atomic)
    derived_path = (folder / args.derived) if not Path(args.derived).is_absolute() else Path(args.derived)
    reactant_types_path = (
        (folder / args.reactant_types) if not Path(args.reactant_types).is_absolute() else Path(args.reactant_types)
    )
    reaction_types_path = (
        (folder / args.reaction_types) if not Path(args.reaction_types).is_absolute() else Path(args.reaction_types)
    )

    atomic_defs, derived_defs = load_feature_definitions(atomic_path=atomic_path, derived_path=derived_path)
    reactant_types = json.loads(reactant_types_path.read_text(encoding="utf-8"))
    reaction_types = json.loads(reaction_types_path.read_text(encoding="utf-8"))

    if not args.no_validate:
        warnings = validate_features(
            atomic_defs,
            derived_defs,
            reactant_types=reactant_types,
            reaction_types=reaction_types,
        )
        for w in warnings:
            print(f"warning: {w}", file=sys.stderr)

    result: dict[str, Any] = classify_reaction_smiles(
        args.reactant,
        atomic_defs=atomic_defs,
        derived_defs=derived_defs,
        reactant_types=reactant_types,
        reaction_types=reaction_types,
    )

    if args.json:
        def _default(obj: Any) -> Any:
            if hasattr(obj, "to_dict"):
                return obj.to_dict()
            if hasattr(obj, "__dict__"):
                return obj.__dict__
            return str(obj)

        print(json.dumps(result, ensure_ascii=False, indent=2, default=_default))
        return 0

    matches = result.get("reaction_type_matches") or []
    if matches:
        print("Matched reaction types:")
        for m in matches:
            print(f"- {m.id}: {m.name}")
            why = sorted(set(m.why_all_of) | set(m.why_any_of))
            if why:
                print(f"  why: {', '.join(why)}")
    else:
        print("Matched reaction types: []")

    key_tokens = [
        "sp2_electrophile_present",
        "suzuki_boron_partner_present",
        "buchwald_amine_partner_present",
    ]
    reaction_has = result.get("reaction_has") or {}
    key_summary = {k: bool(reaction_has.get(k, False)) for k in key_tokens}
    print("Key mixture tokens:", key_summary)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
