#!/usr/bin/env python
"""
Command-line interface for UnifiedRecommender.

Quick tool for testing reaction condition recommendations using the unified
protocol + rule index.

Usage:
    # Basic usage
    python scripts/unified_recommend_cli.py "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    # With options
    python scripts/unified_recommend_cli.py "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" --top-k 3 --min-sim 0.5
    
    # Filter by protocol or rule
    python scripts/unified_recommend_cli.py "reaction.smi" --type protocol
    python scripts/unified_recommend_cli.py "reaction.smi" --type rule
"""

import argparse
import sys
import json
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools.recommend import UnifiedRecommender


def format_recommendation(result, show_details=False):
    """Format a single recommendation for CLI display."""
    icon = "📋" if result.source_type == "protocol" else "📖"
    
    lines = [
        f"{icon} [{result.rank}] {result.name}",
        f"    Family: {result.family}",
        f"    Similarity: {result.similarity:.3f}",
        f"    Type: {result.source_type}",
        f"    Version: {result.version}",
    ]
    
    if result.tags:
        tags_str = ", ".join(result.tags[:5])  # Limit to 5 tags
        if len(result.tags) > 5:
            tags_str += f" (+{len(result.tags) - 5} more)"
        lines.append(f"    Tags: {tags_str}")
    
    if show_details:
        lines.append(f"    ID: {result.id}")
        lines.append(f"    Source: {result.source_file}")
    
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Find reaction condition recommendations using DRFP similarity",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic Buchwald-Hartwig recommendation
  %(prog)s "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
  
  # Top 3 protocols only with minimum similarity 0.5
  %(prog)s "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" -k 3 --min-sim 0.5 --type protocol
  
  # Read from file
  %(prog)s reaction.smi -k 5
  
  # JSON output for scripting
  %(prog)s "reaction_smiles" --json
        """
    )
    
    parser.add_argument(
        "reaction",
        help="Reaction SMILES string or path to file containing SMILES"
    )
    parser.add_argument(
        "-k", "--top-k",
        type=int,
        default=5,
        help="Number of recommendations to return (default: 5)"
    )
    parser.add_argument(
        "--min-sim",
        type=float,
        default=0.0,
        help="Minimum similarity threshold 0.0-1.0 (default: 0.0)"
    )
    parser.add_argument(
        "--type",
        choices=["protocol", "rule", "all"],
        default="all",
        help="Filter by source type (default: all)"
    )
    parser.add_argument(
        "--index-dir",
        type=str,
        default=None,
        help="Path to custom unified index directory (optional)"
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output results as JSON"
    )
    parser.add_argument(
        "--details",
        action="store_true",
        help="Show additional details (IDs, file paths)"
    )
    parser.add_argument(
        "--stats",
        action="store_true",
        help="Show index statistics before recommendations"
    )
    
    args = parser.parse_args()
    
    # Read reaction SMILES
    reaction_smiles = args.reaction
    if Path(reaction_smiles).exists():
        with open(reaction_smiles, 'r') as f:
            reaction_smiles = f.read().strip()
    
    # Initialize recommender
    try:
        if args.index_dir:
            recommender = UnifiedRecommender(index_dir=args.index_dir)
        else:
            recommender = UnifiedRecommender()
    except Exception as e:
        print(f"Error loading recommender: {e}", file=sys.stderr)
        return 1
    
    # Show statistics if requested
    if args.stats:
        stats = recommender.get_statistics()
        print("=" * 70)
        print("Index Statistics")
        print("=" * 70)
        print(f"Version: {stats['build_info']['version']}")
        print(f"Protocols: {stats['protocols']['count']}")
        print(f"Rules: {stats['rules']['count']}")
        print(f"Fingerprints: {stats['drfp']['computed']}")
        print()
    
    # Get recommendations
    source_types = None
    if args.type != "all":
        source_types = [args.type]
    
    try:
        results = recommender.recommend(
            reaction_smiles=reaction_smiles,
            top_k=args.top_k,
            min_similarity=args.min_sim,
            source_types=source_types
        )
    except Exception as e:
        print(f"Error getting recommendations: {e}", file=sys.stderr)
        return 1
    
    # Output results
    if args.json:
        # JSON output for scripting
        output = {
            "query": reaction_smiles,
            "count": len(results),
            "filters": {
                "top_k": args.top_k,
                "min_similarity": args.min_sim,
                "source_type": args.type
            },
            "recommendations": [
                {
                    "rank": r.rank,
                    "id": r.id,
                    "name": r.name,
                    "source_type": r.source_type,
                    "family": r.family,
                    "similarity": round(r.similarity, 3),
                    "tags": r.tags,
                    "version": r.version,
                    "source_file": r.source_file
                }
                for r in results
            ]
        }
        print(json.dumps(output, indent=2))
    else:
        # Human-readable output
        print("=" * 70)
        print("Unified Condition Recommender")
        print("=" * 70)
        print(f"Query: {reaction_smiles}")
        print(f"Filters: top_k={args.top_k}, min_sim={args.min_sim}, type={args.type}")
        print()
        
        if not results:
            print("No recommendations found matching criteria.")
            print("Try lowering min_similarity or removing type filter.")
        else:
            print(f"Found {len(results)} recommendation(s):")
            print()
            for result in results:
                print(format_recommendation(result, show_details=args.details))
                print()
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
