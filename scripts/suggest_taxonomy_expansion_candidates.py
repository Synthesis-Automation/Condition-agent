#!/usr/bin/env python3
"""
Suggest taxonomy expansion candidates from unknown-cluster motifs.

Input is expected from scripts/triage_unknown_clusters.py output JSON.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Set, Tuple

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.taxonomy.reaction_catalog import load_reaction_catalog


def _split_signature_side(text: str) -> Set[str]:
    value = str(text or "").strip()
    if not value or value == "none":
        return set()
    return {tok.strip() for tok in value.split("|") if tok.strip() and tok.strip() != "none"}


def _parse_motif_signature(signature: str) -> Tuple[Set[str], Set[str]]:
    text = str(signature or "").strip()
    if "->" not in text:
        return set(), set()
    left, right = text.split("->", 1)
    return _split_signature_side(left), _split_signature_side(right)


def _taxonomy_reaction_motif_sets() -> Dict[str, Dict[str, Any]]:
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
            "name": defn.name,
            "category": defn.category,
            "reactants": reactants,
            "products": products,
            "reactant_slot_count": len(defn.reactants),
            "product_slot_count": len(defn.products),
        }
    return payload


def _overlap_ratio(query: Set[str], target: Set[str]) -> float:
    if not query:
        return 0.0
    return len(query & target) / max(1, len(query))


def _score_candidate(
    reacted: Set[str],
    formed: Set[str],
    candidate: Dict[str, Any],
) -> Dict[str, Any]:
    reactant_ratio = _overlap_ratio(reacted, set(candidate.get("reactants") or set()))
    product_target = set(candidate.get("products") or set())
    product_ratio = _overlap_ratio(formed, product_target) if formed else 0.0

    if formed:
        score = 0.65 * reactant_ratio + 0.35 * product_ratio
    else:
        score = 0.85 * reactant_ratio

    return {
        "reaction_id": candidate.get("reaction_id"),
        "name": candidate.get("name"),
        "category": candidate.get("category"),
        "score": round(score, 4),
        "reactant_overlap": round(reactant_ratio, 4),
        "product_overlap": round(product_ratio, 4),
        "matched_reacted_motifs": sorted(list(reacted & set(candidate.get("reactants") or set()))),
        "matched_formed_motifs": sorted(list(formed & product_target)),
    }


def _recommended_action(top_score: float, cluster_count: int) -> str:
    if top_score >= 0.75 and cluster_count >= 3:
        return "expand_existing_taxonomy_motif_sets_or_aliases"
    if top_score >= 0.45:
        return "review_slot_constraints_and_missing_product_patterns"
    return "consider_new_reaction_family_or_multistep_record_bucket"


def suggest_candidates(
    clusters_payload: Dict[str, Any],
    *,
    top_candidates: int = 5,
) -> Dict[str, Any]:
    taxonomy = _taxonomy_reaction_motif_sets()
    clusters = clusters_payload.get("clusters") or []
    output_clusters: List[Dict[str, Any]] = []

    for cluster in clusters:
        motif_signature = str(cluster.get("motif_signature") or "")
        reacted, formed = _parse_motif_signature(motif_signature)
        scored: List[Dict[str, Any]] = []
        for candidate in taxonomy.values():
            scored.append(_score_candidate(reacted, formed, candidate))
        scored.sort(key=lambda row: (float(row.get("score") or 0.0), float(row.get("reactant_overlap") or 0.0)), reverse=True)
        top = scored[: max(1, int(top_candidates))]
        top_score = float(top[0].get("score") or 0.0) if top else 0.0
        count = int(cluster.get("count") or 0)

        output_clusters.append(
            {
                "cluster_id": cluster.get("cluster_id"),
                "count": count,
                "motif_signature": motif_signature,
                "event_signature": cluster.get("event_signature"),
                "source_reaction_labels": cluster.get("source_reaction_labels") or [],
                "recommended_action": _recommended_action(top_score, count),
                "top_taxonomy_candidates": top,
            }
        )

    return {
        "input_cluster_count": len(clusters),
        "taxonomy_reaction_type_count": len(taxonomy),
        "clusters": output_clusters,
    }


def write_flat_csv(path: Path, payload: Dict[str, Any]) -> None:
    fields = [
        "cluster_id",
        "cluster_count",
        "motif_signature",
        "event_signature",
        "recommended_action",
        "candidate_rank",
        "candidate_reaction_id",
        "candidate_name",
        "candidate_category",
        "score",
        "reactant_overlap",
        "product_overlap",
        "matched_reacted_motifs",
        "matched_formed_motifs",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for cluster in payload.get("clusters") or []:
            for rank, candidate in enumerate(cluster.get("top_taxonomy_candidates") or [], start=1):
                writer.writerow(
                    {
                        "cluster_id": cluster.get("cluster_id"),
                        "cluster_count": cluster.get("count"),
                        "motif_signature": cluster.get("motif_signature"),
                        "event_signature": cluster.get("event_signature"),
                        "recommended_action": cluster.get("recommended_action"),
                        "candidate_rank": rank,
                        "candidate_reaction_id": candidate.get("reaction_id"),
                        "candidate_name": candidate.get("name"),
                        "candidate_category": candidate.get("category"),
                        "score": candidate.get("score"),
                        "reactant_overlap": candidate.get("reactant_overlap"),
                        "product_overlap": candidate.get("product_overlap"),
                        "matched_reacted_motifs": "|".join(candidate.get("matched_reacted_motifs") or []),
                        "matched_formed_motifs": "|".join(candidate.get("matched_formed_motifs") or []),
                    }
                )


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Suggest taxonomy expansion candidates from unknown clusters")
    parser.add_argument("--input", required=True, help="Unknown cluster JSON")
    parser.add_argument("--output-json", default="results/taxonomy_expansion_candidates.json", help="Output JSON path")
    parser.add_argument("--output-csv", default="results/taxonomy_expansion_candidates.csv", help="Output CSV path")
    parser.add_argument("--top-candidates", type=int, default=5, help="Top taxonomy candidates per cluster")
    parser.add_argument("--top-clusters", type=int, default=50, help="Only evaluate first N clusters from input")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    input_path = Path(args.input).expanduser().resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    if isinstance(payload.get("clusters"), list):
        payload["clusters"] = payload["clusters"][: max(1, int(args.top_clusters))]

    out = suggest_candidates(payload, top_candidates=max(1, int(args.top_candidates)))

    output_json = Path(args.output_json).expanduser().resolve()
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")

    output_csv = Path(args.output_csv).expanduser().resolve()
    write_flat_csv(output_csv, out)

    print(f"Taxonomy candidate JSON: {output_json}")
    print(f"Taxonomy candidate CSV: {output_csv}")
    print(f"Clusters evaluated: {out['input_cluster_count']}")
    print(f"Taxonomy reaction types considered: {out['taxonomy_reaction_type_count']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
