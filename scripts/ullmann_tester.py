#!/usr/bin/env python3
"""
Ullmann Featurizer Tester
=========================

Inspect features extracted for electrophile/nucleophile SMILES pairs.

Examples:
  python scripts/ullmann_tester.py "Brc1ccc(F)cc1" "Nc1ccccc1"
  python scripts/ullmann_tester.py -f pairs.txt --jsonl
  python scripts/ullmann_tester.py --demo --pretty --field bin

File format (pairs.txt): one pair per line, separated by tab, comma, pipe, or whitespace.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Iterable, Tuple, List

# Ensure local package import without installation
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from chemtools.featurizers.ullmann import featurize  # type: ignore


def _split_pair(line: str) -> Tuple[str, str] | None:
    s = line.strip()
    if not s:
        return None
    for sep in ("\t", ",", "|"):
        if sep in s:
            a, b = [p.strip() for p in s.split(sep, 1)]
            return (a, b)
    parts = s.split()
    if len(parts) >= 2:
        return (parts[0], parts[1])
    return None


def iter_pairs(args) -> Iterable[Tuple[str, str]]:
    if args.file:
        with open(args.file, "r", encoding="utf-8") as f:
            for line in f:
                pr = _split_pair(line)
                if pr:
                    yield pr
        return
    if args.demo:
        demo: List[Tuple[str, str]] = [
            ("Brc1ccc(F)cc1", "Nc1ccccc1"),      # aryl bromide + aniline
            ("Clc1ccncc1", "Nc1ccccc1"),         # 4-chloropyridine + aniline
            ("OS(=O)(=O)C(F)(F)F-c1ccccc1", "Nc1ccccc1"),  # aryl triflate + aniline (illustrative)
            ("Brc1ccccc1", "c1ccc(\"O\")cc1"),         # aryl bromide + phenol (string contains quotes; demo only)
        ]
        for p in demo:
            yield p
        return
    if args.smiles:
        if len(args.smiles) % 2 != 0:
            raise SystemExit("Provide an even number of SMILES: electrophile nucleophile [electrophile nucleophile]...")
        it = iter(args.smiles)
        for e, n in zip(it, it):  # type: ignore
            yield (e, n)
        return
    # Interactive
    print("Enter electrophile and nucleophile (separated by space). Blank line to quit:")
    while True:
        try:
            line = input("> ")
        except EOFError:
            break
        pr = _split_pair(line)
        if not pr:
            break
        yield pr


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Test Ullmann featurizer on SMILES pairs")
    p.add_argument("smiles", nargs="*", help="Electrophile and nucleophile SMILES in pairs")
    p.add_argument("--file", "-f", help="File with one SMILES pair per line")
    p.add_argument("--demo", action="store_true", help="Run on a small demo set")
    p.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    p.add_argument("--jsonl", action="store_true", help="Emit one compact JSON object per line")
    p.add_argument("--field", choices=[
        "LG", "elec_class", "ortho_count", "para_EWG", "heteroaryl",
        "nuc_class", "n_basicity", "steric_alpha", "bin"
    ], help="Print only a single feature per pair")
    args = p.parse_args(argv)

    any_done = False
    for e, n in iter_pairs(args):
        any_done = True
        feats = featurize(e, n)
        if args.field:
            print(feats.get(args.field, ""))
            continue
        payload = {
            "electrophile": e,
            "nucleophile": n,
            **feats,
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

