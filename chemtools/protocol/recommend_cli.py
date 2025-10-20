#!/usr/bin/env python3
"""
Protocol Recommendation CLI

Find matching protocols for a reaction using DRFP similarity.

Usage:
    # Recommend protocols for a reaction
    python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1"
    
    # Specify number of results
    python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --k 5
    
    # Filter by reaction family
    python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --family Suzuki
    
    # Filter by tags
    python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --tags coupling,Pd
    
    # Disable SMARTS filtering
    python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --no-smarts-filter
    
    # Show legacy format (for debugging)
    python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --legacy-format
    
    # Export to JSON file
    python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --output results.json
    
    # Pretty print JSON
    python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --pretty
"""

import sys
import json
import argparse
import logging
from pathlib import Path
from typing import Optional, List

from .recommend import ProtocolRecommender

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s'
)
logger = logging.getLogger(__name__)


def print_recommendations(results: dict, pretty: bool = False):
    """Print recommendations in human-readable format"""
    meta = results.get('meta', {})
    detection = results.get('detection', {})
    recommendations = results.get('recommended_conditions', [])
    
    print("=" * 70)
    print("Protocol Recommendations")
    print("=" * 70)
    print()
    
    # Meta info
    print(f"Model: {meta.get('model_type', 'Unknown')}")
    print(f"Status: {meta.get('status', 'Unknown')}")
    print(f"Processing time: {meta.get('processing_time_ms', 0):.1f} ms")
    print()
    
    # Detection info
    print(f"Detected type: {detection.get('detected_type', 'Unknown')}")
    if detection.get('confidence'):
        print(f"Confidence: {detection['confidence']:.3f}")
    print()
    
    # Input
    input_data = results.get('input', {})
    if input_data.get('reaction_smiles'):
        print(f"Reaction: {input_data['reaction_smiles']}")
        print()
    
    # Recommendations
    if not recommendations:
        print("No matching protocols found.")
        return
    
    print(f"Found {len(recommendations)} matching protocol(s):")
    print()
    
    for rec in recommendations:
        rank = rec.get('rank', '?')
        confidence = rec.get('confidence', 0)
        similarity = rec.get('similarity', confidence)
        metadata = rec.get('protocol_metadata', {})
        
        print(f"{'=' * 70}")
        print(f"Rank {rank} - Similarity: {similarity:.3f}")
        print(f"{'=' * 70}")
        
        # Protocol info
        title = metadata.get('title', 'Untitled')
        journal = metadata.get('journal', '')
        year = metadata.get('year', '')
        
        print(f"Title: {title}")
        if journal:
            print(f"Journal: {journal}", end='')
            if year:
                print(f" ({year})", end='')
            print()
        
        family = metadata.get('reaction_family', '')
        if family:
            print(f"Family: {family}")
        
        # Conditions
        conditions = rec.get('conditions', {})
        if conditions:
            print()
            print("Conditions:")
            if conditions.get('catalyst'):
                print(f"  Catalyst: {conditions['catalyst']}")
            if conditions.get('ligand'):
                print(f"  Ligand: {conditions['ligand']}")
            if conditions.get('base'):
                print(f"  Base: {conditions['base']}")
            if conditions.get('solvent'):
                print(f"  Solvent: {conditions['solvent']}")
            if conditions.get('temperature_C'):
                print(f"  Temperature: {conditions['temperature_C']} °C")
            if conditions.get('time_h'):
                print(f"  Time: {conditions['time_h']} h")
            if conditions.get('atmosphere'):
                print(f"  Atmosphere: {conditions['atmosphere']}")
        
        # File reference
        filename = metadata.get('filename', '')
        if filename:
            print()
            print(f"Source file: {filename}")
        
        print()


def main(argv: Optional[List[str]] = None) -> int:
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Find matching protocols for a reaction using DRFP similarity",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        'reaction_smiles',
        help="Reaction SMILES string (e.g., 'CCBr.c1ccccc1B(O)O>>CCc1ccccc1')"
    )
    
    parser.add_argument(
        '--k', '-n',
        type=int,
        default=5,
        help="Number of top recommendations to return (default: 5)"
    )
    
    parser.add_argument(
        '--family', '-f',
        help="Filter by reaction family (e.g., 'Suzuki', 'Buchwald')"
    )
    
    parser.add_argument(
        '--tags', '-t',
        help="Filter by tags (comma-separated, e.g., 'coupling,Pd')"
    )
    
    parser.add_argument(
        '--min-similarity',
        type=float,
        default=0.0,
        help="Minimum similarity threshold (0.0-1.0, default: 0.0)"
    )
    
    parser.add_argument(
        '--no-smarts-filter',
        action='store_true',
        help="Disable SMARTS structural pre-filtering"
    )
    
    parser.add_argument(
        '--legacy-format',
        action='store_true',
        help="Use legacy output format instead of standard format"
    )
    
    parser.add_argument(
        '--index-path',
        type=Path,
        help="Path to protocol index file (default: auto-detect)"
    )
    
    parser.add_argument(
        '--output', '-o',
        type=Path,
        help="Save results to JSON file"
    )
    
    parser.add_argument(
        '--pretty', '-p',
        action='store_true',
        help="Pretty print JSON output"
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help="Show verbose output"
    )
    
    args = parser.parse_args(argv)
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Parse tags
    tags = None
    if args.tags:
        tags = [t.strip() for t in args.tags.split(',')]
    
    # Create recommender
    try:
        recommender = ProtocolRecommender(index_path=args.index_path)
    except Exception as e:
        print(f"❌ Failed to load protocol recommender: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1
    
    # Get recommendations
    try:
        results = recommender.recommend(
            reaction_smiles=args.reaction_smiles,
            k=args.k,
            reaction_family=args.family,
            tags=tags,
            min_similarity=args.min_similarity,
            use_standard_format=not args.legacy_format,
            use_smarts_filter=not args.no_smarts_filter
        )
    except Exception as e:
        print(f"❌ Recommendation failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1
    
    # Save to file if requested
    if args.output:
        try:
            with open(args.output, 'w', encoding='utf-8') as f:
                json.dump(results, f, indent=(2 if args.pretty else None), ensure_ascii=False)
            print(f"✅ Results saved to: {args.output}")
            print()
        except Exception as e:
            print(f"❌ Failed to save results: {e}")
            return 1
    
    # Print results
    if args.pretty or not args.output:
        if args.legacy_format:
            # Print JSON for legacy format
            print(json.dumps(results, indent=2, ensure_ascii=False))
        else:
            # Print human-readable format for standard format
            print_recommendations(results, args.pretty)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
