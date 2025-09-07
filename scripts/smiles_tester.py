#!/usr/bin/env python3
"""
Simple tester for the SMILES and reaction SMILES normalizer.

Usage examples:
  python scripts/smiles_tester.py "Brc1ccc(F)cc1.O"
  python scripts/smiles_tester.py --reaction "Brc1ccc(F)cc1.Nc1ccccc1>>"
  python scripts/smiles_tester.py --demo
  python scripts/smiles_tester.py --file smiles.txt

If no inputs are provided, enters interactive mode.
"""
from __future__ import annotations
import argparse
import json
import os
import sys
from typing import Iterable


# Ensure we can import the local package without installation
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from chemtools.smiles import normalize, normalize_reaction  # type: ignore


def iter_inputs(args) -> Iterable[str]:
    if args.file:
        with open(args.file, "r", encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if s:
                    yield s
        return
    if args.demo:
        for s in [
            "Brc1ccc(F)cc1.O",            # aryl bromide + water
            "Nc1ccccc1.Cl",                # anilinium chloride salt
            "O=C([O-])c1ccccc1",          # benzoate -> benzoic acid
            "C[N+](C)(C)C",               # tetramethylammonium
            "notasmiles",                  # invalid when RDKit available
            "Brc1ccc(F)cc1.Nc1ccccc1>>",  # reaction SMILES (empty agents)
        ]:
            yield s
        return
    if args.smiles:
        for s in args.smiles:
            yield s
        return
    # Interactive
    print("Enter SMILES (blank line to quit):")
    while True:
        try:
            s = input("> ").strip()
        except EOFError:
            break
        if not s:
            break
        yield s


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Test the SMILES normalizer")
    p.add_argument("smiles", nargs="*", help="One or more SMILES/reaction SMILES strings")
    p.add_argument("--file", "-f", help="File with one SMILES per line")
    p.add_argument("--demo", action="store_true", help="Run on a small demo set")
    p.add_argument("--reaction", action="store_true", help="Treat inputs as reaction SMILES (split by '>')")
    args = p.parse_args(argv)

    any_done = False
    for s in iter_inputs(args):
        any_done = True
        if args.reaction or ('>' in s):
            res = normalize_reaction(s)
        else:
            res = normalize(s)
        print(json.dumps({"input": s, **res}, indent=2, ensure_ascii=False))
    if not any_done:
        p.print_help()
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

