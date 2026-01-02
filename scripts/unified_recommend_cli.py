#!/usr/bin/env python
"""
CLI for the unified recommendation engine.

Searches across reaction datasets, protocols, and HTE sources.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools.recommend import UnifiedRecommender


def format_recommendation(result, show_details: bool = False) -> str:
    icon = "P" if result.source_type == "protocol" else ("D" if result.source_type == "dataset" else "H")
    lines = [
        f"{icon} [{result.rank}] {result.name}",
        f"    Family: {result.family}",
        f"    Similarity: {result.similarity:.3f}",
        f"    Type: {result.source_type}",
    ]
    if result.drfp_similarity is not None:
        lines.append(f"    DRFP: {result.drfp_similarity:.3f}")
    if result.feature_similarity is not None:
        lines.append(f"    Features: {result.feature_similarity:.3f}")
    if show_details:
        lines.append(f"    ID: {result.id}")
        lines.append(f"    Source: {result.source_file}")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Find reaction condition recommendations using unified DRFP + feature similarity",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("reaction", help="Reaction SMILES string or path to file containing SMILES")
    parser.add_argument("-k", "--top-k", type=int, default=5, help="Number of recommendations to return (default: 5)")
    parser.add_argument("--min-sim", type=float, default=0.0, help="Minimum similarity threshold 0.0-1.0")
    parser.add_argument(
        "--type",
        choices=["protocol", "dataset", "hte", "all"],
        default="all",
        help="Filter by source type (default: all)",
    )
    parser.add_argument("--index-dir", type=str, default=None, help="Path to custom index directory")
    parser.add_argument("--json", action="store_true", help="Output results as JSON")
    parser.add_argument("--details", action="store_true", help="Show IDs and file paths")
    parser.add_argument("--stats", action="store_true", help="Show index statistics before recommendations")

    args = parser.parse_args()

    reaction_smiles = args.reaction
    if Path(reaction_smiles).exists():
        reaction_smiles = Path(reaction_smiles).read_text(encoding="utf-8").strip()

    try:
        recommender = UnifiedRecommender(index_dir=args.index_dir) if args.index_dir else UnifiedRecommender()
    except Exception as exc:
        print(f"Error loading recommender: {exc}", file=sys.stderr)
        return 1

    if args.stats:
        stats = recommender.get_statistics()
        print("=" * 70)
        print("Index Statistics")
        print("=" * 70)
        print(f"Total entries: {stats.get('total_entries')}")
        by_source = stats.get("by_source", {})
        if by_source:
            print(f"Protocols: {by_source.get('protocol', 0)}")
            print(f"Datasets: {by_source.get('dataset', 0)}")
            print(f"HTE: {by_source.get('hte', 0)}")
        drfp = stats.get("drfp", {})
        if drfp:
            print(f"Fingerprints: {drfp.get('computed', 0)}")
        print()

    source_types = None if args.type == "all" else [args.type]
    try:
        results = recommender.recommend(
            reaction_smiles=reaction_smiles,
            top_k=args.top_k,
            min_similarity=args.min_sim,
            source_types=source_types,
        )
    except Exception as exc:
        print(f"Error getting recommendations: {exc}", file=sys.stderr)
        return 1

    if args.json:
        output = {
            "query": reaction_smiles,
            "count": len(results),
            "filters": {
                "top_k": args.top_k,
                "min_similarity": args.min_sim,
                "source_type": args.type,
            },
            "recommendations": [
                {
                    "rank": r.rank,
                    "id": r.id,
                    "name": r.name,
                    "source_type": r.source_type,
                    "family": r.family,
                    "similarity": round(r.similarity, 3),
                    "drfp_similarity": r.drfp_similarity,
                    "feature_similarity": r.feature_similarity,
                    "source_file": r.source_file,
                }
                for r in results
            ],
        }
        print(json.dumps(output, indent=2))
    else:
        print("=" * 70)
        print("Unified Condition Recommender")
        print("=" * 70)
        print(f"Query: {reaction_smiles}")
        print(f"Filters: top_k={args.top_k}, min_sim={args.min_sim}, type={args.type}")
        print()
        if not results:
            print("No recommendations found matching criteria.")
        else:
            print(f"Found {len(results)} recommendation(s):")
            print()
            for result in results:
                print(format_recommendation(result, show_details=args.details))
                print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
