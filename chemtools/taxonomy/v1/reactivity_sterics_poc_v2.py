#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CLI: reactivity sterics v2 (Ar- anchor, ortho bulk -> 0–10 score).

Example:
  python reactivity_sterics_poc_v2.py "Cc1ccccc1Br" "Brc1ccccc1" "O=Cc1ccccc1"
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

try:
    from .ar_context_sterics_v2 import analyze_smiles_reactivity_sterics_v2
except Exception:
    from ar_context_sterics_v2 import analyze_smiles_reactivity_sterics_v2


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("smiles", nargs="+", help="One or more SMILES to analyze.")
    ap.add_argument("--compiled", default="calculable_features.compiled.v1.json", help="Compiled features JSON.")
    ap.add_argument("--groups", default="organic_groups.v1.json", help="Groups JSON.")
    args = ap.parse_args()

    compiled = Path(args.compiled)
    groups = Path(args.groups)

    results = [
        analyze_smiles_reactivity_sterics_v2(s, compiled_features_path=compiled, groups_path=groups)
        for s in args.smiles
    ]
    print(json.dumps({"results": results}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

