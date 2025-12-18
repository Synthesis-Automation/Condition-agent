#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CLI: reactivity electronics POC v1 (Ar- electron-poor score via Gasteiger delta charge).

Example:
  python reactivity_electronics_poc_v1.py "Brc1ccccc1" "COc1ccc(Br)cc1" "O=[N+]([O-])c1ccc(Br)cc1"
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

try:
    from .ar_context_electronics_v1 import analyze_smiles_reactivity_electronics_v1
except Exception:
    from ar_context_electronics_v1 import analyze_smiles_reactivity_electronics_v1


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("smiles", nargs="+", help="One or more SMILES to analyze.")
    ap.add_argument("--compiled", default="specs/calculable_features.compiled.v1.json", help="Compiled features JSON.")
    ap.add_argument("--groups", default="specs/organic_groups.v1.json", help="Groups JSON.")
    ap.add_argument("--k", default=25.0, type=float, help="Scaling factor for score mapping (raw = 5 + k*delta_q).")
    args = ap.parse_args()

    compiled = Path(args.compiled)
    groups = Path(args.groups)

    results = [
        analyze_smiles_reactivity_electronics_v1(
            s, compiled_features_path=compiled, groups_path=groups, k=float(args.k)
        )
        for s in args.smiles
    ]
    print(json.dumps({"results": results}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
