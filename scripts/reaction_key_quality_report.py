#!/usr/bin/env python3
"""
Dataset-level diagnostics for reaction key reliability.

Usage examples:
  python scripts/reaction_key_quality_report.py --input data/reactions.txt
  python scripts/reaction_key_quality_report.py --input data/reactions.csv --column reaction_smiles
  python scripts/reaction_key_quality_report.py --input data/reactions.jsonl --column rxn
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.featurizers.unified import featurize_reaction


def _detect_format(path: Path, explicit: str) -> str:
    if explicit != "auto":
        return explicit
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return "csv"
    if suffix in {".jsonl", ".ndjson"}:
        return "jsonl"
    return "txt"


def _iter_reaction_smiles(
    path: Path,
    *,
    input_format: str,
    column: str,
    delimiter: str,
) -> Iterator[str]:
    if input_format == "csv":
        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            for row in reader:
                value = str((row or {}).get(column, "")).strip()
                if value:
                    yield value
        return
    if input_format == "jsonl":
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                text = line.strip()
                if not text:
                    continue
                try:
                    payload = json.loads(text)
                except Exception:
                    continue
                if isinstance(payload, str):
                    value = payload.strip()
                elif isinstance(payload, dict):
                    value = str(payload.get(column, "")).strip()
                else:
                    value = ""
                if value:
                    yield value
        return
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            yield text


def _reaction_type_label(result: Dict[str, Any]) -> str:
    value = result.get("reaction_type")
    if isinstance(value, dict):
        return str(value.get("reaction_type") or value.get("name") or "Unknown")
    text = str(value or "").strip()
    return text if text else "Unknown"


def _extract_key_quality(result: Dict[str, Any]) -> Dict[str, Any]:
    events = result.get("reaction_events")
    if isinstance(events, dict):
        quality = events.get("reaction_key_quality")
        if isinstance(quality, dict) and quality:
            return quality
    detection = result.get("detection")
    if isinstance(detection, dict):
        quality = detection.get("reaction_key_quality")
        if isinstance(quality, dict) and quality:
            return quality
    return {}


def _motif_signature(result: Dict[str, Any]) -> str:
    aggregates = result.get("aggregates") or {}
    reacted = sorted(str(m) for m in (aggregates.get("reacted_motifs") or []) if m)
    formed = sorted(str(m) for m in (aggregates.get("formed_motifs") or []) if m)
    reacted_sig = "|".join(reacted[:4]) if reacted else "none"
    formed_sig = "|".join(formed[:4]) if formed else "none"
    return f"{reacted_sig} -> {formed_sig}"


def build_report(reactions: Iterable[str], *, limit: Optional[int] = None) -> Dict[str, Any]:
    total = 0
    unknown_type = 0
    empty_key = 0
    low_quality = 0
    error_count = 0
    quality_levels: Counter[str] = Counter()
    anomaly_counts: Counter[str] = Counter()
    unresolved_signatures: Counter[str] = Counter()

    run_options = {
        "detailed": False,
        "include_reaction_type": True,
        "confirm_coupling_products": True,
    }

    for idx, reaction_smiles in enumerate(reactions):
        if limit is not None and idx >= max(0, int(limit)):
            break
        rxn = str(reaction_smiles or "").strip()
        if not rxn:
            continue
        total += 1
        try:
            result = featurize_reaction(rxn, options=run_options)
        except Exception:
            error_count += 1
            continue

        reaction_key = str(result.get("reaction_key") or "").strip()
        if not reaction_key:
            empty_key += 1

        reaction_type = _reaction_type_label(result)
        if reaction_type.lower() == "unknown":
            unknown_type += 1

        quality = _extract_key_quality(result)
        level = str(quality.get("level") or "unknown").strip().lower()
        quality_levels[level] += 1
        score = float(quality.get("score_0_1") or 0.0)
        if level == "low" or score < 0.45:
            low_quality += 1

        events = result.get("reaction_events")
        anomalies = (events or {}).get("anomalies") if isinstance(events, dict) else []
        for anomaly in anomalies or []:
            label = str(anomaly).strip()
            if label:
                anomaly_counts[label] += 1

        if reaction_type.lower() == "unknown" or level in {"low", "unknown"}:
            unresolved_signatures[_motif_signature(result)] += 1

    denominator = max(total, 1)
    return {
        "total_reactions": total,
        "error_count": error_count,
        "unknown_reaction_type": {
            "count": unknown_type,
            "rate": round(unknown_type / denominator, 4),
        },
        "empty_reaction_key": {
            "count": empty_key,
            "rate": round(empty_key / denominator, 4),
        },
        "low_reaction_key_quality": {
            "count": low_quality,
            "rate": round(low_quality / denominator, 4),
        },
        "quality_levels": dict(sorted(quality_levels.items())),
        "top_anomalies": anomaly_counts.most_common(10),
        "top_unresolved_motif_signatures": unresolved_signatures.most_common(15),
    }


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Reaction key quality diagnostics report")
    parser.add_argument("--input", required=True, help="Input file (txt, csv, jsonl)")
    parser.add_argument(
        "--format",
        default="auto",
        choices=["auto", "txt", "csv", "jsonl"],
        help="Input format (default: auto by extension)",
    )
    parser.add_argument(
        "--column",
        default="reaction_smiles",
        help="Column/key used for CSV or JSONL dict payloads",
    )
    parser.add_argument(
        "--delimiter",
        default=",",
        help="CSV delimiter (default: ',')",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional max number of reactions to process",
    )
    parser.add_argument(
        "--output",
        default="",
        help="Optional JSON output path",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    input_path = Path(args.input).expanduser().resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    input_format = _detect_format(input_path, str(args.format))
    reactions = list(
        _iter_reaction_smiles(
            input_path,
            input_format=input_format,
            column=str(args.column),
            delimiter=str(args.delimiter),
        )
    )
    report = build_report(reactions, limit=args.limit)
    output = json.dumps(report, indent=2, sort_keys=True)
    print(output)

    output_path = str(args.output or "").strip()
    if output_path:
        target = Path(output_path).expanduser().resolve()
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(output, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
