#!/usr/bin/env python3
"""
ConditionCore Tester
====================

Quickly parse a condition core from a list of reagents or dataset rows.

Examples:
  # Minimal Cu/DMEDA
  python scripts/condition_core_tester.py "142-71-2:CATALYST" "110-70-3:LIGAND"

  # Using names/aliases (role inferred via registry if omitted)
  python scripts/condition_core_tester.py CuI LIGAND:DMEDA BASE:K3PO4

  # From JSON file (array or object with {reagents: [...]})
  python scripts/condition_core_tester.py -f reagents.json --pretty

  # From Ullmann dataset JSONL (first 3 entries)
  python scripts/condition_core_tester.py --dataset data/reaction_dataset/Ullman-C-N.jsonl --limit 3 --jsonl

If no inputs are provided, enters interactive mode.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Dict, Any, Iterable, List, Tuple

# Ensure local package import without installation
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from chemtools.condition_core import parse  # type: ignore
from chemtools.registry import resolve as reg_resolve  # type: ignore


def _norm_role(s: str | None) -> str | None:
    if not s:
        return None
    s = s.strip().upper()
    m = {
        "CAT": "CATALYST",
        "METAL": "CATALYST",
        "CATALYST": "CATALYST",
        "LIG": "LIGAND",
        "LIGAND": "LIGAND",
        "BASE": "BASE",
        "SOLV": "SOLVENT",
        "SOLVENT": "SOLVENT",
        "ADDITIVE": "ADDITIVE",
    }
    return m.get(s, s)


def _reagent_from_token(token: str) -> Dict[str, Any]:
    # token forms:
    #   "value:ROLE" or "ROLE:value" or just "value"
    raw = token.strip()
    role = None
    value = raw
    if ":" in raw:
        a, b = raw.split(":", 1)
        # decide which side is role by presence in role set
        ra = _norm_role(a)
        rb = _norm_role(b)
        if ra in {"CATALYST", "LIGAND", "BASE", "SOLVENT", "ADDITIVE"}:
            role = ra
            value = b
        elif rb in {"CATALYST", "LIGAND", "BASE", "SOLVENT", "ADDITIVE"}:
            role = rb
            value = a
        else:
            value = raw
    value = value.strip()

    info = reg_resolve(value)
    uid = info.get("uid") if isinstance(info, dict) else None
    name = info.get("name") if isinstance(info, dict) else None
    if not role and isinstance(info, dict):
        role = info.get("role")
    return {
        "uid": uid or value,
        "role": role or "ADDITIVE",
        "name": name or value,
    }


def _load_json(path: str) -> Any:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _iter_dataset_rows(path: str, limit: int | None) -> Iterable[Dict[str, Any]]:
    count = 0
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            try:
                obj = json.loads(line)
            except Exception:
                continue
            yield obj
            count += 1
            if limit and count >= limit:
                break


def _reagents_from_dataset_row(obj: Dict[str, Any]) -> List[Dict[str, Any]]:
    reagents: List[Dict[str, Any]] = []
    # Catalyst block
    cat = (obj.get("catalyst") or {})
    for lst_key in ("full_system", "core"):
        for it in (cat.get(lst_key) or []):
            uid = (it.get("cas") or it.get("uid") or "").strip()
            name = (it.get("name") or "").strip()
            info = reg_resolve(uid or name)
            role = info.get("role") if isinstance(info, dict) else None
            reagents.append({"uid": uid or name, "name": name or uid, "role": role or "CATALYST"})
    # Other reagents
    for it in (obj.get("reagents") or []):
        reagents.append({"uid": (it.get("cas") or "").strip(), "name": it.get("name"), "role": (it.get("role") or "").upper() or "ADDITIVE"})
    for it in (obj.get("solvents") or []):
        reagents.append({"uid": (it.get("cas") or "").strip(), "name": it.get("name"), "role": "SOLVENT"})
    return reagents


def iter_inputs(args) -> Iterable[Tuple[List[Dict[str, Any]], str]]:
    # Returns sequence of (reagents, text)
    if args.dataset:
        for row in _iter_dataset_rows(args.dataset, args.limit):
            yield _reagents_from_dataset_row(row), (args.text or "")
        return
    if args.file:
        data = _load_json(args.file)
        if isinstance(data, dict) and "reagents" in data:
            yield data["reagents"], (data.get("text") or args.text or "")
        elif isinstance(data, list):
            yield data, (args.text or "")
        else:
            raise SystemExit("Unsupported JSON format. Use an array of reagents or {reagents:[...]}")
        return
    if args.json:
        raw = args.json
        # Be robust to Windows PowerShell: single quotes are not special to the OS argv parser,
        # so they may be included literally. Try several trims before parsing.
        s = raw.strip()
        try:
            data = json.loads(s)
        except Exception:
            s2 = s
            if s2 and s2[0] in "'\"":
                s2 = s2[1:]
            if s2 and s2[-1] in "'\"":
                s2 = s2[:-1]
            try:
                data = json.loads(s2)
            except Exception:
                # Last resort: unescape common quote-escape artifacts
                s3 = s2.replace('\"', '"').replace("\'", "'")
                data = json.loads(s3)
        if isinstance(data, dict) and "reagents" in data:
            yield data["reagents"], (data.get("text") or args.text or "")
        elif isinstance(data, list):
            yield data, (args.text or "")
        else:
            raise SystemExit("Unsupported JSON format for --json input")
        return
    if args.items:
        reagents = [_reagent_from_token(tok) for tok in args.items]
        yield reagents, (args.text or "")
        return
    # Interactive
    print("Enter reagents as 'value:ROLE' tokens separated by spaces. Blank line to quit.")
    print("Example: 142-71-2:CATALYST 110-70-3:LIGAND 7778-53-2:BASE")
    while True:
        try:
            line = input("> ").strip()
        except EOFError:
            break
        if not line:
            break
        tokens = [t for t in line.split() if t]
        reagents = [_reagent_from_token(tok) for tok in tokens]
        yield reagents, (args.text or "")


def main(argv: List[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Parse condition core from reagents")
    p.add_argument("items", nargs="*", help="Tokens of the form value:ROLE or ROLE:value or just value (role inferred)")
    p.add_argument("--file", "-f", help="Path to JSON file (array of reagents or {reagents:[...]})")
    p.add_argument("--json", help="Inline JSON (array or {reagents:[...]})")
    p.add_argument("--dataset", help="Path to Ullmann dataset JSONL to sample from")
    p.add_argument("--limit", type=int, default=0, help="Limit dataset rows processed (0=all)")
    p.add_argument("--text", help="Optional free-text hint")
    p.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    p.add_argument("--jsonl", action="store_true", help="Emit one compact JSON object per line")
    p.add_argument("--field", choices=["core", "metal_source_uid", "ligand_uid", "precatalyst"], help="Print only this field")
    args = p.parse_args(argv)

    lim = args.limit or None
    any_done = False
    for reagents, text in iter_inputs(args):
        any_done = True
        res = parse(reagents, text)
        if args.field:
            print(res.get(args.field, ""))
            continue
        payload = {"reagents": reagents, "text": text, **res}
        if args.jsonl and not args.pretty:
            print(json.dumps(payload, ensure_ascii=False))
        else:
            print(json.dumps(payload, indent=(2 if args.pretty else None), ensure_ascii=False))
    if not any_done:
        p.print_help()
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
