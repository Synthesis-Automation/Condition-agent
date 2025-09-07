#!/usr/bin/env python3
"""
Reaction-family Router Tester
=============================

Probe reaction-family detection from reactant SMILES or reaction SMILES.

Examples:
  # Reactants as reaction SMILES (preferred)
  python scripts/router_tester.py --reaction "Brc1ccc(F)cc1.Nc1ccccc1>>"
  # Reactants as dot-separated SMILES
  python scripts/router_tester.py "Brc1ccc(F)cc1.Nc1ccccc1"
  # From file
  python scripts/router_tester.py -f reactions.txt --reaction --jsonl
  # Demo set
  python scripts/router_tester.py --demo --pretty
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Iterable, List

# Ensure local package import
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from chemtools.router import detect_family  # type: ignore
from chemtools.smiles import normalize_reaction  # type: ignore


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
            # Ullmann C–N
            "Brc1ccc(F)cc1.Nc1ccccc1>>",
            # Suzuki
            "Brc1ccccc1.c1ccc(B(O)O)cc1>>",
            # Sonogashira
            "Brc1ccccc1.C#CC>>",
            # Amide coupling
            "c1ccc(C(=O)O)cc1.NCC>>",
        ]:
            yield s
        return
    if args.queries:
        for s in args.queries:
            yield s
        return
    # Interactive
    print("Enter reactants (dot-separated) or reaction SMILES, blank to quit:")
    while True:
        try:
            s = input("> ").strip()
        except EOFError:
            break
        if not s:
            break
        yield s


def to_reactant_list(s: str, treat_as_reaction: bool) -> List[str]:
    if treat_as_reaction or (">" in s):
        res = normalize_reaction(s)
        out: List[str] = []
        for it in res.get("reactants", []):
            smi = it.get("smiles_norm") or it.get("largest_smiles") or it.get("input") or ""
            if smi:
                out.append(smi)
        return out
    # Dot-separated reactants, or space-separated fallback
    if "." in s:
        parts = [p for p in s.split(".") if p]
    else:
        parts = [p for p in s.split() if p]
    return parts


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Test reaction-family detection from SMILES")
    p.add_argument("queries", nargs="*", help="Dot-separated reactants or reaction SMILES")
    p.add_argument("--file", "-f", help="File with one query per line")
    p.add_argument("--demo", action="store_true", help="Run a small demo set")
    p.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    p.add_argument("--jsonl", action="store_true", help="Emit one compact JSON object per line")
    p.add_argument("--reaction", action="store_true", help="Treat inputs as reaction SMILES")
    p.add_argument("--field", choices=["family", "confidence"], help="Print only a single field per query")
    args = p.parse_args(argv)

    any_done = False
    for s in iter_inputs(args):
        any_done = True
        reactants = to_reactant_list(s, args.reaction)
        result = detect_family(reactants)
        if args.field:
            print(result.get(args.field, ""))
            continue
        payload = {
            "query": s,
            "reactants": reactants,
            **result,
        }
        if args.jsonl and not args.pretty:
            print(json.dumps(payload, ensure_ascii=False))
        else:
            print(json.dumps(payload, indent=(2 if args.pretty else None), ensure_ascii=False))

    if not any_done:
        p.print_help()
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

