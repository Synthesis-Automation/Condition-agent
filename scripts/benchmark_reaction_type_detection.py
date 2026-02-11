#!/usr/bin/env python3
"""
Benchmark legacy vs specificity-aware reaction type validation.

The script runs reaction featurization to produce CRK keys, then evaluates:
- Legacy matcher (first-hit taxonomy order)
- New matcher (specificity-aware ranked taxonomy match)

Usage:
    python scripts/benchmark_reaction_type_detection.py
    python scripts/benchmark_reaction_type_detection.py --input data/reaction_dataset/C_N_Coupling.csv
    python scripts/benchmark_reaction_type_detection.py --limit 200 --out-csv results/benchmark_rows.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from chemtools.featurizers.formatters.detection_validation import (
    validate_detection_with_crk_key,
)
from chemtools.featurizers.unified import featurize_reaction
from chemtools.taxonomy.reaction_catalog import (
    load_reaction_catalog,
)
from chemtools.util.ingestion_routing import classify_reaction_for_taxonomy_benchmark


REACTION_COL_CANDIDATES = (
    "reaction_smiles",
    "reaction",
    "rxn_smiles",
    "smiles",
)
TRUTH_COL_CANDIDATES = (
    "reaction_type_standardized",
    "reaction_type",
    "rxn_type",
    "reaction_family",
)


def _normalize_type(label: Any) -> str:
    text = str(label or "").strip()
    if not text:
        return ""
    definitions, alias_map = load_reaction_catalog()
    if text in definitions:
        return text
    return alias_map.get(text.lower(), text)


def _pick_column(fieldnames: Iterable[str], candidates: Tuple[str, ...]) -> Optional[str]:
    names = {str(n).strip(): str(n).strip() for n in (fieldnames or [])}
    lowered = {k.lower(): k for k in names}
    for candidate in candidates:
        col = lowered.get(candidate.lower())
        if col:
            return col
    return None


def _iter_input_files(input_value: Optional[str]) -> List[Path]:
    if input_value:
        p = Path(input_value)
        if p.is_file():
            return [p]
        if p.is_dir():
            return sorted(x for x in p.rglob("*.csv") if x.is_file())
        return []
    default_dir = PROJECT_ROOT / "data" / "reaction_dataset"
    if default_dir.exists():
        return sorted(x for x in default_dir.rglob("*.csv") if x.is_file())
    return []


def _extract_crk_raw(payload: Dict[str, Any]) -> str:
    detection = payload.get("detection")
    if isinstance(detection, dict):
        validation = detection.get("validation")
        if isinstance(validation, dict):
            raw = str(validation.get("reaction_key_raw") or "").strip()
            if raw:
                return raw
    return str(payload.get("reaction_key") or "").strip()


def _safe_featurize(reaction_smiles: str) -> Optional[Dict[str, Any]]:
    try:
        return featurize_reaction(
            reaction_smiles,
            options={"detailed": True, "confirm_coupling_products": True},
        )
    except Exception:
        return None


def _benchmark_file(
    path: Path,
    *,
    limit: int,
    routing_policy: str,
) -> Dict[str, Any]:
    rows_out: List[Dict[str, Any]] = []
    total = 0
    usable = 0
    legacy_correct = 0
    new_correct = 0
    routed_excluded = 0
    route_counts: Counter[str] = Counter()

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle)
        reaction_col = _pick_column(reader.fieldnames or [], REACTION_COL_CANDIDATES)
        truth_col = _pick_column(reader.fieldnames or [], TRUTH_COL_CANDIDATES)
        if not reaction_col or not truth_col:
            return {
                "file": str(path),
                "error": f"missing required columns reaction_col={reaction_col}, truth_col={truth_col}",
                "rows": rows_out,
                "total": 0,
                "usable": 0,
                "legacy_correct": 0,
                "new_correct": 0,
            }

        for idx, row in enumerate(reader, start=1):
            if limit > 0 and usable >= limit:
                break
            total += 1
            reaction_smiles = str(row.get(reaction_col) or "").strip()
            truth = _normalize_type(row.get(truth_col))
            if not reaction_smiles or not truth:
                continue
            route = classify_reaction_for_taxonomy_benchmark(reaction_smiles)
            route_name = str(route.get("route") or "")
            if route_name:
                route_counts[route_name] += 1
            if routing_policy == "exclude_complex" and bool(route.get("excluded")):
                routed_excluded += 1
                continue
            payload = _safe_featurize(reaction_smiles)
            if not payload:
                continue
            crk_raw = _extract_crk_raw(payload)
            if not crk_raw:
                continue

            legacy = validate_detection_with_crk_key(
                initial_detection="Unknown",
                initial_confidence=0.0,
                reaction_key=crk_raw,
                use_legacy=True,
                include_evidence=False,
            ).get("reaction_type")
            new = validate_detection_with_crk_key(
                initial_detection="Unknown",
                initial_confidence=0.0,
                reaction_key=crk_raw,
                use_legacy=False,
                include_evidence=False,
            ).get("reaction_type")

            legacy_norm = _normalize_type(legacy)
            new_norm = _normalize_type(new)
            legacy_hit = legacy_norm == truth
            new_hit = new_norm == truth

            usable += 1
            if legacy_hit:
                legacy_correct += 1
            if new_hit:
                new_correct += 1

            rows_out.append(
                {
                    "file": str(path),
                    "row_index": idx,
                    "reaction_smiles": reaction_smiles,
                    "truth": truth,
                    "legacy_prediction": legacy_norm,
                    "new_prediction": new_norm,
                    "legacy_correct": legacy_hit,
                    "new_correct": new_hit,
                }
            )

    return {
        "file": str(path),
        "rows": rows_out,
        "total": total,
        "usable": usable,
        "routed_excluded": routed_excluded,
        "routing_route_counts": route_counts.most_common(),
        "legacy_correct": legacy_correct,
        "new_correct": new_correct,
    }


def _ratio(num: int, den: int) -> float:
    return float(num) / float(den) if den > 0 else 0.0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Benchmark legacy vs specificity-aware reaction type matching.",
    )
    parser.add_argument(
        "--input",
        help="CSV file or directory (default: data/reaction_dataset).",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Max usable rows per file (0 = no limit).",
    )
    parser.add_argument(
        "--out-csv",
        help="Optional output CSV for per-row predictions.",
    )
    parser.add_argument(
        "--routing-policy",
        default="exclude_complex",
        choices=["none", "exclude_complex"],
        help="Ingestion routing policy before taxonomy benchmarking.",
    )
    args = parser.parse_args()

    files = _iter_input_files(args.input)
    if not files:
        print("No input CSV files found.")
        return 1

    all_rows: List[Dict[str, Any]] = []
    grand_total = 0
    grand_usable = 0
    grand_routed_excluded = 0
    grand_legacy = 0
    grand_new = 0

    print(f"Benchmarking {len(files)} file(s)...")
    for path in files:
        result = _benchmark_file(
            path,
            limit=max(0, int(args.limit)),
            routing_policy=str(args.routing_policy),
        )
        if result.get("error"):
            print(f"- {path.name}: {result['error']}")
            continue

        usable = int(result["usable"])
        routed_excluded = int(result.get("routed_excluded") or 0)
        legacy_correct = int(result["legacy_correct"])
        new_correct = int(result["new_correct"])
        legacy_acc = _ratio(legacy_correct, usable)
        new_acc = _ratio(new_correct, usable)
        delta = new_acc - legacy_acc

        print(
            f"- {path.name}: usable={usable}, routed_excluded={routed_excluded}, "
            f"legacy={legacy_acc:.3f}, new={new_acc:.3f}, delta={delta:+.3f}"
        )

        grand_total += int(result["total"])
        grand_usable += usable
        grand_routed_excluded += routed_excluded
        grand_legacy += legacy_correct
        grand_new += new_correct
        all_rows.extend(result["rows"])

    if grand_usable == 0:
        print("No usable rows were evaluated.")
        return 1

    grand_legacy_acc = _ratio(grand_legacy, grand_usable)
    grand_new_acc = _ratio(grand_new, grand_usable)
    print("\nOverall:")
    print(f"  total rows scanned: {grand_total}")
    print(f"  routed excluded rows: {grand_routed_excluded}")
    print(f"  usable rows: {grand_usable}")
    print(f"  legacy accuracy: {grand_legacy_acc:.4f}")
    print(f"  new accuracy:    {grand_new_acc:.4f}")
    print(f"  delta:           {grand_new_acc - grand_legacy_acc:+.4f}")

    if args.out_csv:
        out_path = Path(args.out_csv)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = [
            "file",
            "row_index",
            "reaction_smiles",
            "truth",
            "legacy_prediction",
            "new_prediction",
            "legacy_correct",
            "new_correct",
        ]
        with out_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for row in all_rows:
                writer.writerow(row)
        print(f"Wrote per-row output: {out_path}")

    # Emit machine-readable summary line
    print(
        json.dumps(
            {
                "routing_policy": str(args.routing_policy),
                "usable": grand_usable,
                "routed_excluded": grand_routed_excluded,
                "legacy_accuracy": grand_legacy_acc,
                "new_accuracy": grand_new_acc,
                "delta": grand_new_acc - grand_legacy_acc,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
