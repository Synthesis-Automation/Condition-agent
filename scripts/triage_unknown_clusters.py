#!/usr/bin/env python3
"""
Triage unknown reaction-type predictions into motif-signature clusters.

Example:
  python scripts/triage_unknown_clusters.py \
    --input examples/test_reactions_random_sampled.csv \
    --output-json results/unknown_clusters.random_sampled.json \
    --output-samples-csv results/unknown_cluster_samples.random_sampled.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.featurizers.unified import featurize_reaction


def _reaction_type_label(result: Dict[str, Any]) -> str:
    value = result.get("reaction_type")
    if isinstance(value, dict):
        return str(value.get("reaction_type") or value.get("name") or "Unknown")
    text = str(value or "").strip()
    return text if text else "Unknown"


def _motif_signature(result: Dict[str, Any]) -> str:
    aggregates = result.get("aggregates") or {}
    reacted = sorted(str(m) for m in (aggregates.get("reacted_motifs") or []) if m)
    formed = sorted(str(m) for m in (aggregates.get("formed_motifs") or []) if m)
    reacted_sig = "|".join(reacted[:6]) if reacted else "none"
    formed_sig = "|".join(formed[:6]) if formed else "none"
    return f"{reacted_sig} -> {formed_sig}"


def _event_signature(result: Dict[str, Any]) -> str:
    events = (result.get("reaction_events") or {}).get("events") if isinstance(
        result.get("reaction_events"),
        dict,
    ) else []
    kinds = sorted(
        str(ev.get("kind"))
        for ev in (events or [])
        if isinstance(ev, dict) and ev.get("kind")
    )
    return "|".join(kinds[:6]) if kinds else "none"


def _quality_level(result: Dict[str, Any]) -> str:
    events = result.get("reaction_events") or {}
    quality = events.get("reaction_key_quality") if isinstance(events, dict) else {}
    if isinstance(quality, dict):
        level = str(quality.get("level") or "").strip().lower()
        if level:
            return level
    return "unknown"


def _quality_score(result: Dict[str, Any]) -> Optional[float]:
    events = result.get("reaction_events") or {}
    quality = events.get("reaction_key_quality") if isinstance(events, dict) else {}
    if isinstance(quality, dict):
        try:
            return float(quality.get("score_0_1"))
        except Exception:
            return None
    return None


def _cluster_key(result: Dict[str, Any]) -> Tuple[str, str]:
    return _motif_signature(result), _event_signature(result)


def _iter_csv_rows(path: Path, reaction_col: str) -> Iterable[Tuple[int, Dict[str, str]]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        for row_num, row in enumerate(reader, start=2):
            if not isinstance(row, dict):
                continue
            rxn = str(row.get(reaction_col) or "").strip()
            if not rxn:
                continue
            yield row_num, row


def _reservoir_insert(
    reservoir: List[Dict[str, Any]],
    candidate: Dict[str, Any],
    *,
    seen_count: int,
    capacity: int,
    rng: random.Random,
) -> None:
    if capacity <= 0:
        return
    if len(reservoir) < capacity:
        reservoir.append(candidate)
        return
    pick = rng.randint(1, max(1, seen_count))
    if pick <= capacity:
        reservoir[pick - 1] = candidate


def build_unknown_clusters(
    rows: Iterable[Tuple[int, Dict[str, str]]],
    *,
    sample_per_cluster: int = 5,
    seed: int = 20260211,
    max_rows: Optional[int] = None,
) -> Dict[str, Any]:
    rng = random.Random(seed)
    run_options = {
        "detailed": False,
        "include_reaction_type": True,
        "confirm_coupling_products": True,
    }

    processed = 0
    unknown = 0
    errors = 0
    clusters: Dict[str, Dict[str, Any]] = {}

    for idx, (row_num, row) in enumerate(rows):
        if max_rows is not None and idx >= max(0, int(max_rows)):
            break
        reaction_smiles = str(row.get("reaction_smiles") or "").strip()
        if not reaction_smiles:
            continue
        processed += 1
        try:
            result = featurize_reaction(reaction_smiles, options=run_options)
        except Exception:
            errors += 1
            continue

        predicted = _reaction_type_label(result)
        if predicted.lower() != "unknown":
            continue
        unknown += 1

        motif_sig, event_sig = _cluster_key(result)
        cluster_id = f"{motif_sig} || events:{event_sig}"
        cluster = clusters.setdefault(
            cluster_id,
            {
                "cluster_id": cluster_id,
                "motif_signature": motif_sig,
                "event_signature": event_sig,
                "count": 0,
                "quality_levels": Counter(),
                "anomalies": Counter(),
                "source_files": Counter(),
                "source_reaction_labels": Counter(),
                "samples": [],
            },
        )
        cluster["count"] += 1
        cluster["quality_levels"][_quality_level(result)] += 1

        anomalies = (result.get("reaction_events") or {}).get("anomalies") if isinstance(
            result.get("reaction_events"),
            dict,
        ) else []
        for anomaly in anomalies or []:
            text = str(anomaly).strip()
            if text:
                cluster["anomalies"][text] += 1

        source_file = str(row.get("source_file") or row.get("source_reaction_family") or "").strip()
        source_label = str(row.get("reaction_type") or "").strip()
        if source_file:
            cluster["source_files"][source_file] += 1
        if source_label:
            cluster["source_reaction_labels"][source_label] += 1

        sample_payload = {
            "original_row_number": row_num,
            "reaction_id": str(row.get("reaction_id") or "").strip(),
            "source_file": source_file,
            "source_reaction_label": source_label,
            "reaction_smiles": reaction_smiles,
            "reaction_key": str(result.get("reaction_key") or ""),
            "quality_level": _quality_level(result),
            "quality_score": _quality_score(result),
        }
        _reservoir_insert(
            cluster["samples"],
            sample_payload,
            seen_count=cluster["count"],
            capacity=sample_per_cluster,
            rng=rng,
        )

    ranked = sorted(clusters.values(), key=lambda entry: int(entry.get("count", 0)), reverse=True)
    out_clusters: List[Dict[str, Any]] = []
    for entry in ranked:
        out_clusters.append(
            {
                "cluster_id": entry["cluster_id"],
                "motif_signature": entry["motif_signature"],
                "event_signature": entry["event_signature"],
                "count": int(entry["count"]),
                "quality_levels": dict(sorted(entry["quality_levels"].items())),
                "anomalies": entry["anomalies"].most_common(10),
                "source_files": entry["source_files"].most_common(10),
                "source_reaction_labels": entry["source_reaction_labels"].most_common(10),
                "samples": entry["samples"],
            }
        )

    return {
        "processed_reactions": processed,
        "unknown_reactions": unknown,
        "unknown_rate": round(unknown / max(1, processed), 4),
        "errors": errors,
        "cluster_count": len(out_clusters),
        "clusters": out_clusters,
    }


def write_samples_csv(path: Path, clusters: List[Dict[str, Any]]) -> None:
    fields = [
        "cluster_rank",
        "cluster_count",
        "cluster_id",
        "motif_signature",
        "event_signature",
        "original_row_number",
        "reaction_id",
        "source_file",
        "source_reaction_label",
        "quality_level",
        "quality_score",
        "reaction_key",
        "reaction_smiles",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for rank, cluster in enumerate(clusters, start=1):
            for sample in cluster.get("samples") or []:
                row = {
                    "cluster_rank": rank,
                    "cluster_count": cluster.get("count"),
                    "cluster_id": cluster.get("cluster_id"),
                    "motif_signature": cluster.get("motif_signature"),
                    "event_signature": cluster.get("event_signature"),
                    "original_row_number": sample.get("original_row_number"),
                    "reaction_id": sample.get("reaction_id"),
                    "source_file": sample.get("source_file"),
                    "source_reaction_label": sample.get("source_reaction_label"),
                    "quality_level": sample.get("quality_level"),
                    "quality_score": sample.get("quality_score"),
                    "reaction_key": sample.get("reaction_key"),
                    "reaction_smiles": sample.get("reaction_smiles"),
                }
                writer.writerow(row)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Cluster unknown reaction-type predictions")
    parser.add_argument("--input", required=True, help="Input CSV")
    parser.add_argument("--reaction-column", default="reaction_smiles", help="Reaction SMILES column")
    parser.add_argument("--output-json", default="results/unknown_clusters.json", help="Output JSON path")
    parser.add_argument(
        "--output-samples-csv",
        default="results/unknown_cluster_samples.csv",
        help="Output sample CSV path",
    )
    parser.add_argument("--top-clusters", type=int, default=50, help="Number of clusters to keep in outputs")
    parser.add_argument("--sample-per-cluster", type=int, default=5, help="Samples stored per cluster")
    parser.add_argument("--seed", type=int, default=20260211, help="Random seed for reservoir sampling")
    parser.add_argument("--limit", type=int, default=None, help="Optional max rows to process")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    input_path = Path(args.input).expanduser().resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    data = build_unknown_clusters(
        _iter_csv_rows(input_path, reaction_col=str(args.reaction_column)),
        sample_per_cluster=max(1, int(args.sample_per_cluster)),
        seed=int(args.seed),
        max_rows=args.limit,
    )
    top_n = max(1, int(args.top_clusters))
    data["clusters"] = list(data.get("clusters") or [])[:top_n]

    output_json = Path(args.output_json).expanduser().resolve()
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")

    output_csv = Path(args.output_samples_csv).expanduser().resolve()
    write_samples_csv(output_csv, data["clusters"])

    print(f"Unknown clustering JSON: {output_json}")
    print(f"Unknown cluster samples CSV: {output_csv}")
    print(f"Processed reactions: {data['processed_reactions']}")
    print(f"Unknown reactions: {data['unknown_reactions']} (rate={data['unknown_rate']})")
    print(f"Clusters: {data['cluster_count']}, kept top: {len(data['clusters'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
