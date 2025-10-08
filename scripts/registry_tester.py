#!/usr/bin/env python3
"""
Simple tester for the Registry resolver.

Usage examples:
  python scripts/registry_tester.py "108-88-3" "Toluene" PhMe
  python scripts/registry_tester.py --file queries.txt --jsonl
  python scripts/registry_tester.py --demo --pretty
  python scripts/registry_tester.py --data data/registry_sample.jsonl --demo

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

from chemtools.registry import resolve  # type: ignore


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
            "108-88-3",   # toluene (CAS)
            "Toluene",    # name
            "PhMe",       # abbreviation/alias
            "7778-53-2",  # K3PO4 (base)
            "K3PO4",      # token
            "7681-65-4",  # CuI
            "CuI",        # token
            "p-Xylene",   # hyphenated name
            "pxylene",    # punctuation-insensitive alias
        ]:
            yield s
        return
    if args.queries:
        for s in args.queries:
            yield s
        return
    # Interactive
    print("Enter query (CAS/name/alias), blank line to quit:")
    while True:
        try:
            s = input("> ").strip()
        except EOFError:
            break
        if not s:
            break
        yield s


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Test the registry resolver")
    p.add_argument("queries", nargs="*", help="One or more queries (CAS/name/alias)")
    p.add_argument("--file", "-f", help="File with one query per line")
    p.add_argument("--demo", action="store_true", help="Run on a small demo set")
    p.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    p.add_argument("--jsonl", action="store_true", help="Emit one compact JSON object per line")
    p.add_argument("--field", choices=["uid", "role", "name"], help="Print only a single field per line")
    p.add_argument("--data", help="Path to a registry JSONL file (overrides default via CHEMTOOLS_REGISTRY_PATH)")
    p.add_argument("--fail-on-miss", action="store_true", help="Exit nonzero if any NOT_FOUND results")
    args = p.parse_args(argv)

    # Optional registry path override
    if args.data:
        os.environ["CHEMTOOLS_REGISTRY_PATH"] = os.path.abspath(args.data)

    any_done = False
    any_missing = False
    for q in iter_inputs(args):
        any_done = True
        res = resolve(q)
        if args.field:
            if isinstance(res, dict) and "error" not in res:
                print(res.get(args.field, ""))
            else:
                any_missing = True
                print("NOT_FOUND")
            continue
        payload = {"query": q, **res} if isinstance(res, dict) else {"query": q, "result": res}
        if args.jsonl and not args.pretty:
            print(json.dumps(payload, ensure_ascii=False))
        else:
            print(json.dumps(payload, indent=(2 if args.pretty else None), ensure_ascii=False))
        if isinstance(res, dict) and res.get("error"):
            any_missing = True

    if not any_done:
        p.print_help()
    return 1 if (any_missing and args.fail_on_miss) else 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

