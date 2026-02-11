#!/usr/bin/env python3
"""
Discover reaction taxonomy coverage gaps from large reaction datasets.

This script is designed for taxonomy expansion work:
1) Process diverse reaction datasets through deterministic featurization.
2) Quantify unknown reaction types and weak/empty reaction keys.
3) Surface unclassified motifs and motifs outside reaction-type taxonomy slots.
4) Cluster unresolved reactions and suggest closest taxonomy candidates.

Examples:
  python scripts/discover_reaction_coverage_gaps.py --input data-processor/reaction_dataset
  python scripts/discover_reaction_coverage_gaps.py --input data-processor/reaction_dataset --limit 50000
  python scripts/discover_reaction_coverage_gaps.py --input examples/test_reactions.csv --reaction-column reaction_smiles
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, DefaultDict, Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.featurizers.unified import featurize_reaction
from chemtools.taxonomy.reaction_catalog import load_reaction_catalog


REACTION_COL_CANDIDATES: Tuple[str, ...] = (
    "reaction_smiles",
    "reaction",
    "rxn_smiles",
    "rxn",
    "smiles",
)
REACTION_LABEL_CANDIDATES: Tuple[str, ...] = (
    "reaction_type_standardized",
    "reaction_type",
    "reaction_family",
    "source_reaction_family",
)
REACTION_ID_CANDIDATES: Tuple[str, ...] = ("reaction_id", "rxn_id", "id")
DEFAULT_OUTPUT_JSON = "results/reaction_coverage_discovery.json"
DEFAULT_OUTPUT_SAMPLES_CSV = "results/reaction_coverage_discovery_samples.csv"


def _to_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except Exception:
        return default


def _find_column(fieldnames: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    normalized = {str(name).strip().lower(): str(name).strip() for name in fieldnames}
    for candidate in candidates:
        col = normalized.get(str(candidate).strip().lower())
        if col:
            return col
    return None


def _iter_csv_rows(
    path: Path,
    *,
    reaction_column: Optional[str] = None,
    max_rows: Optional[int] = None,
    require_arrow: bool = True,
) -> Iterator[Dict[str, Any]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            return
        reaction_col = reaction_column or _find_column(fieldnames, REACTION_COL_CANDIDATES)
        if not reaction_col:
            return
        label_col = _find_column(fieldnames, REACTION_LABEL_CANDIDATES)
        rid_col = _find_column(fieldnames, REACTION_ID_CANDIDATES)

        emitted = 0
        for row_num, row in enumerate(reader, start=2):
            if max_rows is not None and emitted >= max(0, int(max_rows)):
                break
            if not isinstance(row, dict):
                continue
            rxn = str(row.get(reaction_col) or "").strip()
            if not rxn:
                continue
            if require_arrow and ">>" not in rxn:
                continue
            emitted += 1
            yield {
                "reaction_smiles": rxn,
                "source_file": path.name,
                "source_path": str(path),
                "row_number": row_num,
                "reaction_id": str(row.get(rid_col) or "").strip() if rid_col else "",
                "source_reaction_label": str(row.get(label_col) or "").strip() if label_col else "",
            }


def _iter_input_rows(
    input_path: Path,
    *,
    reaction_column: Optional[str] = None,
    max_rows_per_file: Optional[int] = None,
    require_arrow: bool = True,
) -> Iterator[Dict[str, Any]]:
    if input_path.is_file():
        yield from _iter_csv_rows(
            input_path,
            reaction_column=reaction_column,
            max_rows=max_rows_per_file,
            require_arrow=require_arrow,
        )
        return

    files = sorted(path for path in input_path.rglob("*.csv") if path.is_file())
    for path in files:
        yield from _iter_csv_rows(
            path,
            reaction_column=reaction_column,
            max_rows=max_rows_per_file,
            require_arrow=require_arrow,
        )


def _reaction_type_label(result: Dict[str, Any]) -> str:
    value = result.get("reaction_type")
    if isinstance(value, dict):
        text = str(value.get("reaction_type") or value.get("name") or "").strip()
        return text or "Unknown"
    text = str(value or "").strip()
    return text or "Unknown"


def _extract_quality(result: Dict[str, Any]) -> Dict[str, Any]:
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


def _extract_events(result: Dict[str, Any]) -> List[Dict[str, Any]]:
    events = (result.get("reaction_events") or {}).get("events")
    if isinstance(events, list):
        return [event for event in events if isinstance(event, dict)]
    return []


def _extract_anomalies(result: Dict[str, Any]) -> List[str]:
    anomalies = (result.get("reaction_events") or {}).get("anomalies")
    if not isinstance(anomalies, list):
        return []
    values: List[str] = []
    for anomaly in anomalies:
        text = str(anomaly).strip()
        if text:
            values.append(text)
    return values


def _extract_aggregates(result: Dict[str, Any]) -> Dict[str, Any]:
    aggregates = result.get("aggregates")
    return aggregates if isinstance(aggregates, dict) else {}


def _motif_signature(result: Dict[str, Any], *, max_terms: int = 8) -> str:
    aggregates = _extract_aggregates(result)
    reacted = sorted(str(m) for m in (aggregates.get("reacted_motifs") or []) if str(m).strip())
    formed = sorted(
        str(m)
        for m in (
            aggregates.get("formed_motifs_all")
            or aggregates.get("formed_motifs")
            or []
        )
        if str(m).strip()
    )
    left = "|".join(reacted[:max_terms]) if reacted else "none"
    right = "|".join(formed[:max_terms]) if formed else "none"
    return f"{left} -> {right}"


def _event_signature(result: Dict[str, Any], *, max_terms: int = 8) -> str:
    kinds = sorted(
        str(event.get("kind")).strip()
        for event in _extract_events(result)
        if str(event.get("kind") or "").strip()
    )
    return "|".join(kinds[:max_terms]) if kinds else "none"


def _is_unclassified_motif_id(value: Any) -> bool:
    text = str(value or "").strip()
    if not text:
        return True
    low = text.lower()
    if low == "unknown":
        return True
    if text.startswith("Unclassified-"):
        return True
    return False


def _collect_bundle_motif_ids(bundles: Any) -> List[str]:
    out: List[str] = []
    if not isinstance(bundles, list):
        return out
    for bundle in bundles:
        if not isinstance(bundle, dict):
            continue
        for motif in bundle.get("motifs") or []:
            if not isinstance(motif, dict):
                continue
            mid = str(motif.get("id") or motif.get("compound_id") or "").strip()
            if mid:
                out.append(mid)
    return out


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


def _reaction_taxonomy_sets() -> Dict[str, Dict[str, Any]]:
    definitions, _ = load_reaction_catalog()
    payload: Dict[str, Dict[str, Any]] = {}
    for rid, defn in definitions.items():
        reactants: Set[str] = set()
        products: Set[str] = set()
        for slot in defn.reactants.values():
            reactants.update(str(m) for m in (slot.allowed or []) if str(m).strip())
        for slot in defn.products.values():
            products.update(str(m) for m in (slot.allowed or []) if str(m).strip())
        payload[rid] = {
            "reaction_id": rid,
            "name": str(defn.name or rid),
            "category": str(defn.category or ""),
            "reactants": reactants,
            "products": products,
        }
    return payload


def _overlap_ratio(query: Set[str], target: Set[str]) -> float:
    if not query:
        return 0.0
    return len(query & target) / max(1, len(query))


def _split_signature_side(text: str) -> Set[str]:
    value = str(text or "").strip()
    if not value or value == "none":
        return set()
    return {token.strip() for token in value.split("|") if token.strip() and token.strip() != "none"}


def _parse_motif_signature(signature: str) -> Tuple[Set[str], Set[str]]:
    value = str(signature or "").strip()
    if "->" not in value:
        return set(), set()
    left, right = value.split("->", 1)
    return _split_signature_side(left), _split_signature_side(right)


def _score_taxonomy_candidate(
    *,
    reacted: Set[str],
    formed: Set[str],
    candidate: Dict[str, Any],
) -> Dict[str, Any]:
    reactant_set = set(candidate.get("reactants") or set())
    product_set = set(candidate.get("products") or set())
    reactant_overlap = _overlap_ratio(reacted, reactant_set)
    product_overlap = _overlap_ratio(formed, product_set) if formed else 0.0
    score = (0.65 * reactant_overlap + 0.35 * product_overlap) if formed else (0.85 * reactant_overlap)
    return {
        "reaction_id": candidate.get("reaction_id"),
        "name": candidate.get("name"),
        "category": candidate.get("category"),
        "score": round(score, 4),
        "reactant_overlap": round(reactant_overlap, 4),
        "product_overlap": round(product_overlap, 4),
        "matched_reacted_motifs": sorted(list(reacted & reactant_set)),
        "matched_formed_motifs": sorted(list(formed & product_set)),
    }


def _recommended_action(top_score: float, cluster_count: int) -> str:
    if top_score >= 0.75 and cluster_count >= 3:
        return "expand_existing_slots_or_aliases"
    if top_score >= 0.45:
        return "tighten_or_expand_slot_constraints"
    return "consider_new_reaction_family_or_multistep_bucket"


def build_discovery_report(
    rows: Iterable[Dict[str, Any]],
    *,
    limit: Optional[int] = None,
    sample_per_cluster: int = 5,
    seed: int = 20260211,
    progress_every: int = 0,
) -> Dict[str, Any]:
    rng = random.Random(seed)
    run_options = {
        "detailed": False,
        "include_reaction_type": True,
        "confirm_coupling_products": True,
    }
    taxonomy = _reaction_taxonomy_sets()
    taxonomy_reacted_universe: Set[str] = set()
    taxonomy_formed_universe: Set[str] = set()
    for item in taxonomy.values():
        taxonomy_reacted_universe.update(item.get("reactants") or set())
        taxonomy_formed_universe.update(item.get("products") or set())
    taxonomy_union = taxonomy_reacted_universe | taxonomy_formed_universe

    processed = 0
    errors = 0
    unknown_reaction_type = 0
    empty_reaction_key = 0
    low_reaction_key_quality = 0
    missing_formed_motifs = 0
    unresolved_total = 0
    with_mapping_warning = 0

    quality_levels: Counter[str] = Counter()
    anomalies: Counter[str] = Counter()
    predicted_types: Counter[str] = Counter()
    unresolved_reasons: Counter[str] = Counter()
    unclassified_motifs: Counter[str] = Counter()
    outside_taxonomy_motifs: Counter[str] = Counter()
    unknown_reacted_motifs: Counter[str] = Counter()
    unknown_formed_motifs: Counter[str] = Counter()

    per_file_stats: DefaultDict[str, Dict[str, int]] = defaultdict(
        lambda: {
            "processed": 0,
            "errors": 0,
            "unknown_reaction_type": 0,
            "empty_reaction_key": 0,
            "low_reaction_key_quality": 0,
            "unresolved": 0,
            "with_unclassified_motif": 0,
            "with_outside_taxonomy_motif": 0,
        }
    )

    clusters: Dict[str, Dict[str, Any]] = {}

    for idx, row in enumerate(rows):
        if limit is not None and idx >= max(0, int(limit)):
            break
        reaction_smiles = str(row.get("reaction_smiles") or "").strip()
        if not reaction_smiles:
            continue
        source_file = str(row.get("source_file") or "unknown_source").strip() or "unknown_source"
        source_label = str(row.get("source_reaction_label") or "").strip()

        processed += 1
        per_file_stats[source_file]["processed"] += 1
        if progress_every and processed % max(1, int(progress_every)) == 0:
            print(f"[progress] processed={processed} unknown={unknown_reaction_type} unresolved={unresolved_total}")

        try:
            result = featurize_reaction(reaction_smiles, options=run_options)
        except Exception:
            errors += 1
            per_file_stats[source_file]["errors"] += 1
            continue

        reaction_type = _reaction_type_label(result)
        reaction_key = str(result.get("reaction_key") or "").strip()
        quality = _extract_quality(result)
        quality_level = str(quality.get("level") or "unknown").strip().lower()
        quality_score = _to_float(quality.get("score_0_1"), default=0.0)
        aggregates = _extract_aggregates(result)
        reacted_motifs = sorted(str(m) for m in (aggregates.get("reacted_motifs") or []) if str(m).strip())
        formed_motifs = sorted(
            str(m)
            for m in (
                aggregates.get("formed_motifs_all")
                or aggregates.get("formed_motifs")
                or []
            )
            if str(m).strip()
        )

        predicted_types[reaction_type] += 1
        quality_levels[quality_level] += 1
        for anomaly in _extract_anomalies(result):
            anomalies[anomaly] += 1

        row_reasons: Set[str] = set()
        if reaction_type.lower() == "unknown":
            unknown_reaction_type += 1
            per_file_stats[source_file]["unknown_reaction_type"] += 1
            row_reasons.add("unknown_reaction_type")
            for motif in reacted_motifs:
                unknown_reacted_motifs[motif] += 1
            for motif in formed_motifs:
                unknown_formed_motifs[motif] += 1

        if not reaction_key:
            empty_reaction_key += 1
            per_file_stats[source_file]["empty_reaction_key"] += 1
            row_reasons.add("empty_reaction_key")

        if quality_level == "low" or quality_score < 0.45:
            low_reaction_key_quality += 1
            per_file_stats[source_file]["low_reaction_key_quality"] += 1
            row_reasons.add("low_reaction_key_quality")

        if not formed_motifs:
            missing_formed_motifs += 1
            row_reasons.add("missing_formed_motifs")

        detection = result.get("detection")
        if isinstance(detection, dict) and isinstance(detection.get("mapping_warning"), dict):
            with_mapping_warning += 1
            row_reasons.add("mapping_warning")

        reactant_bundle_motifs = _collect_bundle_motif_ids(result.get("reactants"))
        product_bundle_motifs = _collect_bundle_motif_ids(result.get("products"))
        row_unclassified = sorted(
            {
                motif
                for motif in (reactant_bundle_motifs + product_bundle_motifs)
                if _is_unclassified_motif_id(motif)
            }
        )
        if row_unclassified:
            per_file_stats[source_file]["with_unclassified_motif"] += 1
            row_reasons.add("unclassified_motif")
            for motif in row_unclassified:
                unclassified_motifs[motif] += 1

        row_outside_taxonomy = sorted(
            {
                motif
                for motif in (reacted_motifs + formed_motifs)
                if motif not in taxonomy_union and not _is_unclassified_motif_id(motif)
            }
        )
        if row_outside_taxonomy:
            per_file_stats[source_file]["with_outside_taxonomy_motif"] += 1
            row_reasons.add("motif_outside_reaction_taxonomy")
            for motif in row_outside_taxonomy:
                outside_taxonomy_motifs[motif] += 1

        cluster_reason_keys = {
            "unknown_reaction_type",
            "empty_reaction_key",
            "low_reaction_key_quality",
            "unclassified_motif",
            "motif_outside_reaction_taxonomy",
        }
        unresolved = any(reason in cluster_reason_keys for reason in row_reasons)
        if unresolved:
            unresolved_total += 1
            per_file_stats[source_file]["unresolved"] += 1
            for reason in row_reasons:
                unresolved_reasons[reason] += 1

            motif_sig = _motif_signature(result)
            event_sig = _event_signature(result)
            cluster_id = f"{motif_sig} || events:{event_sig}"
            cluster = clusters.setdefault(
                cluster_id,
                {
                    "cluster_id": cluster_id,
                    "motif_signature": motif_sig,
                    "event_signature": event_sig,
                    "count": 0,
                    "reasons": Counter(),
                    "quality_levels": Counter(),
                    "source_files": Counter(),
                    "source_labels": Counter(),
                    "unclassified_motifs": Counter(),
                    "outside_taxonomy_motifs": Counter(),
                    "samples": [],
                },
            )
            cluster["count"] += 1
            for reason in row_reasons:
                cluster["reasons"][reason] += 1
            cluster["quality_levels"][quality_level] += 1
            if source_file:
                cluster["source_files"][source_file] += 1
            if source_label:
                cluster["source_labels"][source_label] += 1
            for motif in row_unclassified:
                cluster["unclassified_motifs"][motif] += 1
            for motif in row_outside_taxonomy:
                cluster["outside_taxonomy_motifs"][motif] += 1

            sample = {
                "source_file": source_file,
                "source_path": row.get("source_path"),
                "row_number": row.get("row_number"),
                "reaction_id": row.get("reaction_id"),
                "source_reaction_label": source_label,
                "reaction_smiles": reaction_smiles,
                "reaction_type": reaction_type,
                "reaction_key": reaction_key,
                "quality_level": quality_level,
                "quality_score": quality_score,
                "reasons": sorted(row_reasons),
                "unclassified_motifs": row_unclassified,
                "outside_taxonomy_motifs": row_outside_taxonomy,
            }
            _reservoir_insert(
                cluster["samples"],
                sample,
                seen_count=int(cluster["count"]),
                capacity=max(1, int(sample_per_cluster)),
                rng=rng,
            )

    ranked_clusters_raw = sorted(
        clusters.values(),
        key=lambda item: int(item.get("count", 0)),
        reverse=True,
    )
    ranked_clusters: List[Dict[str, Any]] = []
    for cluster in ranked_clusters_raw:
        reacted, formed = _parse_motif_signature(str(cluster.get("motif_signature") or ""))
        candidates: List[Dict[str, Any]] = []
        for candidate in taxonomy.values():
            candidates.append(
                _score_taxonomy_candidate(
                    reacted=reacted,
                    formed=formed,
                    candidate=candidate,
                )
            )
        candidates.sort(
            key=lambda item: (
                float(item.get("score") or 0.0),
                float(item.get("reactant_overlap") or 0.0),
                float(item.get("product_overlap") or 0.0),
            ),
            reverse=True,
        )
        top_candidates = candidates[:5]
        top_score = float(top_candidates[0].get("score") or 0.0) if top_candidates else 0.0
        ranked_clusters.append(
            {
                "cluster_id": cluster["cluster_id"],
                "count": int(cluster["count"]),
                "motif_signature": cluster["motif_signature"],
                "event_signature": cluster["event_signature"],
                "reasons": cluster["reasons"].most_common(10),
                "quality_levels": dict(sorted(cluster["quality_levels"].items())),
                "source_files": cluster["source_files"].most_common(10),
                "source_labels": cluster["source_labels"].most_common(10),
                "unclassified_motifs": cluster["unclassified_motifs"].most_common(10),
                "outside_taxonomy_motifs": cluster["outside_taxonomy_motifs"].most_common(10),
                "recommended_action": _recommended_action(top_score, int(cluster["count"])),
                "top_taxonomy_candidates": top_candidates,
                "samples": cluster["samples"],
            }
        )

    per_file_rows: List[Dict[str, Any]] = []
    for source_file, stats in per_file_stats.items():
        processed_n = int(stats.get("processed", 0))
        denom = max(1, processed_n)
        per_file_rows.append(
            {
                "source_file": source_file,
                "processed": processed_n,
                "errors": int(stats.get("errors", 0)),
                "unknown_reaction_type_rate": round(int(stats.get("unknown_reaction_type", 0)) / denom, 4),
                "low_reaction_key_quality_rate": round(int(stats.get("low_reaction_key_quality", 0)) / denom, 4),
                "empty_reaction_key_rate": round(int(stats.get("empty_reaction_key", 0)) / denom, 4),
                "unresolved_rate": round(int(stats.get("unresolved", 0)) / denom, 4),
                "with_unclassified_motif_rate": round(int(stats.get("with_unclassified_motif", 0)) / denom, 4),
                "with_outside_taxonomy_motif_rate": round(
                    int(stats.get("with_outside_taxonomy_motif", 0)) / denom,
                    4,
                ),
            }
        )
    per_file_rows.sort(key=lambda row: (row["unresolved_rate"], row["processed"]), reverse=True)

    denominator = max(1, processed)
    return {
        "summary": {
            "processed_reactions": processed,
            "featurization_errors": errors,
            "unknown_reaction_type": {
                "count": unknown_reaction_type,
                "rate": round(unknown_reaction_type / denominator, 4),
            },
            "empty_reaction_key": {
                "count": empty_reaction_key,
                "rate": round(empty_reaction_key / denominator, 4),
            },
            "low_reaction_key_quality": {
                "count": low_reaction_key_quality,
                "rate": round(low_reaction_key_quality / denominator, 4),
            },
            "missing_formed_motifs": {
                "count": missing_formed_motifs,
                "rate": round(missing_formed_motifs / denominator, 4),
            },
            "with_mapping_warning": {
                "count": with_mapping_warning,
                "rate": round(with_mapping_warning / denominator, 4),
            },
            "unresolved_reactions": {
                "count": unresolved_total,
                "rate": round(unresolved_total / denominator, 4),
            },
        },
        "taxonomy_coverage": {
            "reaction_types_in_taxonomy": len(taxonomy),
            "reacted_slot_motif_count": len(taxonomy_reacted_universe),
            "formed_slot_motif_count": len(taxonomy_formed_universe),
            "total_slot_motif_count": len(taxonomy_union),
            "top_motifs_outside_reaction_taxonomy": outside_taxonomy_motifs.most_common(50),
        },
        "quality_levels": dict(sorted(quality_levels.items())),
        "top_predicted_reaction_types": predicted_types.most_common(30),
        "unresolved_reason_counts": unresolved_reasons.most_common(20),
        "top_anomalies": anomalies.most_common(20),
        "top_unclassified_motifs": unclassified_motifs.most_common(50),
        "top_unknown_reacted_motifs": unknown_reacted_motifs.most_common(50),
        "top_unknown_formed_motifs": unknown_formed_motifs.most_common(50),
        "files": per_file_rows,
        "clusters": ranked_clusters,
    }


def write_samples_csv(path: Path, clusters: List[Dict[str, Any]]) -> None:
    fields = [
        "cluster_rank",
        "cluster_count",
        "cluster_id",
        "motif_signature",
        "event_signature",
        "recommended_action",
        "source_file",
        "source_path",
        "row_number",
        "reaction_id",
        "source_reaction_label",
        "reaction_type",
        "quality_level",
        "quality_score",
        "reasons",
        "unclassified_motifs",
        "outside_taxonomy_motifs",
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
                    "recommended_action": cluster.get("recommended_action"),
                    "source_file": sample.get("source_file"),
                    "source_path": sample.get("source_path"),
                    "row_number": sample.get("row_number"),
                    "reaction_id": sample.get("reaction_id"),
                    "source_reaction_label": sample.get("source_reaction_label"),
                    "reaction_type": sample.get("reaction_type"),
                    "quality_level": sample.get("quality_level"),
                    "quality_score": sample.get("quality_score"),
                    "reasons": "|".join(sample.get("reasons") or []),
                    "unclassified_motifs": "|".join(sample.get("unclassified_motifs") or []),
                    "outside_taxonomy_motifs": "|".join(sample.get("outside_taxonomy_motifs") or []),
                    "reaction_key": sample.get("reaction_key"),
                    "reaction_smiles": sample.get("reaction_smiles"),
                }
                writer.writerow(row)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Discover reaction taxonomy coverage gaps")
    parser.add_argument(
        "--input",
        default="data-processor/reaction_dataset",
        help="Input CSV file or directory of CSV files",
    )
    parser.add_argument(
        "--reaction-column",
        default="",
        help="Optional explicit reaction SMILES column name (CSV input)",
    )
    parser.add_argument(
        "--max-rows-per-file",
        type=int,
        default=None,
        help="Optional row limit per file before featurization",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional global limit of processed reactions",
    )
    parser.add_argument(
        "--sample-per-cluster",
        type=int,
        default=5,
        help="Reservoir sample size stored for each unresolved cluster",
    )
    parser.add_argument(
        "--top-clusters",
        type=int,
        default=60,
        help="Number of largest unresolved clusters to keep in output",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=500,
        help="Print progress every N processed reactions (0 to disable)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260211,
        help="Random seed for deterministic sampling",
    )
    parser.add_argument(
        "--allow-non-reaction-smiles",
        action="store_true",
        help="Allow rows without '>>' in reaction SMILES column",
    )
    parser.add_argument(
        "--output-json",
        default=DEFAULT_OUTPUT_JSON,
        help=f"Output JSON path (default: {DEFAULT_OUTPUT_JSON})",
    )
    parser.add_argument(
        "--output-samples-csv",
        default=DEFAULT_OUTPUT_SAMPLES_CSV,
        help=f"Output sample CSV path (default: {DEFAULT_OUTPUT_SAMPLES_CSV})",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    input_path = Path(args.input).expanduser().resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input path not found: {input_path}")

    rows = _iter_input_rows(
        input_path,
        reaction_column=str(args.reaction_column).strip() or None,
        max_rows_per_file=args.max_rows_per_file,
        require_arrow=not bool(args.allow_non_reaction_smiles),
    )
    report = build_discovery_report(
        rows,
        limit=args.limit,
        sample_per_cluster=max(1, int(args.sample_per_cluster)),
        seed=int(args.seed),
        progress_every=max(0, int(args.progress_every)),
    )

    top_clusters = max(1, int(args.top_clusters))
    report["clusters"] = list(report.get("clusters") or [])[:top_clusters]

    output_json = Path(args.output_json).expanduser().resolve()
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")

    output_samples = Path(args.output_samples_csv).expanduser().resolve()
    write_samples_csv(output_samples, report.get("clusters") or [])

    summary = report.get("summary") or {}
    processed = summary.get("processed_reactions", 0)
    unknown = (summary.get("unknown_reaction_type") or {}).get("count", 0)
    unresolved = (summary.get("unresolved_reactions") or {}).get("count", 0)
    print(f"Discovery JSON: {output_json}")
    print(f"Discovery samples CSV: {output_samples}")
    print(f"Processed reactions: {processed}")
    print(f"Unknown reaction types: {unknown}")
    print(f"Unresolved reactions: {unresolved}")
    print(f"Clusters kept: {len(report.get('clusters') or [])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
