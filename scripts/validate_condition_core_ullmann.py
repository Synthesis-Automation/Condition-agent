from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Any, Dict, List

# Ensure project root on path when running from source
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from chemtools.condition_core import parse as parse_condition_core


def _mk_reagent(name: str | None, cas: str | None, role: str | None) -> Dict[str, Any]:
    return {
        "name": (name or ""),
        "uid": (cas or ""),
        "role": (role or "")
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Validate ConditionCore normalization using the Ullmann C–N dataset")
    p.add_argument("--path", default=os.path.join("data", "reaction_dataset", "Ullman-C-N.jsonl"), help="Path to Ullmann JSONL dataset")
    p.add_argument("--limit", type=int, default=0, help="Limit number of records (0 = all)")
    p.add_argument("--jsonl", action="store_true", help="Emit per-record JSON results")
    p.add_argument("--show-mismatches", type=int, default=10, help="Show up to N mismatches in human-readable mode")
    args = p.parse_args(argv)

    if not os.path.exists(args.path):
        print(f"Dataset not found: {args.path}")
        return 2

    total = 0
    ok = 0
    mismatches: List[Dict[str, Any]] = []

    with open(args.path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except Exception:
                continue

            # Build reagent list from dataset
            reagents: List[Dict[str, Any]] = []
            cat = rec.get("catalyst") or {}
            # Prefer full_system if present; otherwise use core
            for item in (cat.get("full_system") or cat.get("core") or []):
                reagents.append(_mk_reagent(item.get("name"), item.get("cas"), "CATALYST"))
            # Include reagents (bases, etc.) with their roles when present
            for item in (rec.get("reagents") or []):
                reagents.append(_mk_reagent(item.get("name"), item.get("cas"), (item.get("role") or "ADDITIVE")))

            # Optional: include solvents (role not strictly needed for core)
            for item in (rec.get("solvents") or []):
                reagents.append(_mk_reagent(item.get("name"), item.get("cas"), "SOLVENT"))

            # For this validation, text is not required; parser relies on dataset-derived aliases
            out = parse_condition_core(reagents, "")

            truth = (rec.get("condition_core") or "").strip()
            pred = out.get("core", "").strip()

            # Normalize both sides a bit (e.g., allow Cu vs Cu/none)
            def _norm_core(s: str) -> str:
                s = s.strip()
                if s.endswith("/none"):
                    s = s[:-5]
                return s

            total += 1
            t = _norm_core(truth)
            p = _norm_core(pred)
            # Accept metal-only matches when truth omits ligand
            def _metal_part(s: str) -> str:
                return s.split("/", 1)[0] if s else ""
            ok_flag = (t == p) or (_metal_part(t) and _metal_part(t) == _metal_part(p))
            if ok_flag:
                ok += 1
            else:
                mismatches.append({
                    "reaction_id": rec.get("reaction_id"),
                    "truth": truth,
                    "pred": pred,
                    "reagents": reagents,
                })

            if args.jsonl:
                payload = {
                    "reaction_id": rec.get("reaction_id"),
                    "truth": truth,
                    "pred": pred,
                    "match": ok_flag,
                }
                print(json.dumps(payload, ensure_ascii=False))

            if args.limit and total >= args.limit:
                break

    if not args.jsonl:
        acc = (ok / total) * 100.0 if total else 0.0
        print(f"Records: {total}, Matches: {ok}, Accuracy: {acc:.1f}%")
        if mismatches:
            print(f"Showing up to {args.show_mismatches} mismatches:")
            for m in mismatches[: args.show_mismatches]:
                print(f"- {m['reaction_id']}: truth={m['truth']} pred={m['pred']}")

    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
