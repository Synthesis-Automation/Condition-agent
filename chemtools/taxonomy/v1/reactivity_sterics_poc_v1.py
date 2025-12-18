#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CLI: steric reactivity feature layer (POC) — ortho substitution counts around Ar-* sites.

This uses:
- calculable_features.compiled.v1.json (templated atomic features provide Ar-* anchor sites)
- organic_groups.v1.json (chemist-style labels, including Ar-)

Example:
  python reactivity_sterics_poc_v1.py "BrC1=CC=CC=C1C"
  python reactivity_sterics_poc_v1.py "O=Cc1ccccc1"  # Ar-CHO
"""

from __future__ import annotations
import argparse
import json
from pathlib import Path

from ar_context_sterics_v1 import analyze_smiles_ortho_counts

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("smiles", nargs="+", help="One or more SMILES to analyze.")
    ap.add_argument("--compiled", default="calculable_features.compiled.v1.json", help="Compiled features JSON.")
    ap.add_argument("--groups", default="organic_groups.v1.json", help="Groups JSON.")
    args = ap.parse_args()

    compiled = Path(args.compiled)
    groups = Path(args.groups)

    results = [analyze_smiles_ortho_counts(s, compiled, groups) for s in args.smiles]
    print(json.dumps({"results": results}, indent=2, ensure_ascii=False))

if __name__ == "__main__":
    main()
