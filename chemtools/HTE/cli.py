"""
Command-line interface for HTE-based condition recommendation
"""
import argparse
import sys
from pathlib import Path

from chemtools.HTE import HTERecommender, format_result


def main():
    parser = argparse.ArgumentParser(
        description="HTE-based condition recommendation from reactant SMILES",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # C-N Coupling
  python -m chemtools.HTE.cli -a "c1ccc(Br)cc1" -b "CCN" -k 5
  
  # Suzuki Coupling with palladium catalysts only
  python -m chemtools.HTE.cli -a "c1ccc(Cl)cc1" -b "c1ccc(B(O)O)cc1" --reaction Suzuki --catalyst Pd
  
  # C-N Coupling with copper catalysts only
  python -m chemtools.HTE.cli -a "c1ccc(Br)cc1" -b "c1ccc(N)cc1" --catalyst copper
  
  # Single reactant
  python -m chemtools.HTE.cli -a "c1ccc(Br)cc1"
  
  # Batch processing from file
  python -m chemtools.HTE.cli --batch queries.txt --output results.txt
        """
    )
    
    # Input options
    parser.add_argument(
        '-a', '--reactant-a',
        type=str,
        help='SMILES string for first reactant'
    )
    parser.add_argument(
        '-b', '--reactant-b',
        type=str,
        help='SMILES string for second reactant (optional)'
    )
    parser.add_argument(
        '--batch',
        type=str,
        help='File with reactant pairs (format: SMILES_A SMILES_B per line)'
    )
    
    # Options
    parser.add_argument(
        '-k', '--top-k',
        type=int,
        default=5,
        help='Number of recommendations to return (default: 5)'
    )
    parser.add_argument(
        '--min-exp',
        type=int,
        default=2,
        help='Minimum experiments required per condition (default: 2)'
    )
    parser.add_argument(
        '--reaction',
        type=str,
        help='Filter by reaction type (e.g., Suzuki, C_N_Coupling)'
    )
    parser.add_argument(
        '--catalyst',
        type=str,
        help='Filter by catalyst metal type (e.g., Pd, Cu, Ni, palladium, copper)'
    )
    parser.add_argument(
        '--db-path',
        type=str,
        default='data/HTE_db/HTE_0.jsonl',
        help='Path to HTE database JSONL (default: data/HTE_db/HTE_0.jsonl)'
    )
    
    # Output options
    parser.add_argument(
        '-o', '--output',
        type=str,
        help='Output file path (default: stdout)'
    )
    parser.add_argument(
        '--json',
        action='store_true',
        help='Output as JSON format'
    )
    parser.add_argument(
        '--compact',
        action='store_true',
        help='Compact output (top condition only)'
    )
    parser.add_argument(
        '--stats',
        action='store_true',
        help='Show database statistics only'
    )
    
    args = parser.parse_args()
    
    # Initialize recommender
    try:
        recommender = HTERecommender(args.db_path)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    
    # Show statistics and exit
    if args.stats:
        stats = recommender.get_statistics()
        print("\n📊 HTE DATABASE STATISTICS")
        print("="*60)
        for key, value in stats.items():
            print(f"{key:30s}: {value}")
        return 0
    
    # Determine output stream
    output = sys.stdout
    if args.output:
        output = open(args.output, 'w')
    
    try:
        # Batch processing
        if args.batch:
            return process_batch(
                recommender, args.batch, output, 
                args.top_k, args.min_exp, args.reaction, args.catalyst,
                args.json, args.compact
            )
        
        # Single query
        if not args.reactant_a:
            parser.print_help()
            return 1
        
        result = recommender.recommend(
            reactant_a_smiles=args.reactant_a,
            reactant_b_smiles=args.reactant_b,
            top_k=args.top_k,
            min_experiments=args.min_exp,
            reaction_type_filter=args.reaction,
            catalyst_filter=args.catalyst
        )
        
        if args.json:
            output.write(result_to_json(result))
        elif args.compact:
            output.write(format_compact(result))
        else:
            output.write(format_result(result))
        
        return 0
    
    finally:
        if args.output:
            output.close()


def process_batch(
    recommender, batch_file, output, 
    top_k, min_exp, reaction_filter, catalyst_filter,
    json_format, compact
):
    """Process batch file with reactant pairs"""
    batch_path = Path(batch_file)
    if not batch_path.exists():
        print(f"Error: Batch file not found: {batch_file}", file=sys.stderr)
        return 1
    
    with open(batch_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    
    results = []
    for i, line in enumerate(lines, 1):
        parts = line.split()
        if len(parts) == 0:
            continue
        
        reactant_a = parts[0]
        reactant_b = parts[1] if len(parts) > 1 else None
        
        print(f"Processing {i}/{len(lines)}: {reactant_a} + {reactant_b or 'None'}", 
              file=sys.stderr)
        
        result = recommender.recommend(
            reactant_a_smiles=reactant_a,
            reactant_b_smiles=reactant_b,
            top_k=top_k,
            min_experiments=min_exp,
            reaction_type_filter=reaction_filter,
            catalyst_filter=catalyst_filter
        )
        
        results.append(result)
        
        if not json_format:
            output.write("\n" + "="*80 + "\n")
            output.write(f"QUERY {i}: {reactant_a} + {reactant_b or 'None'}\n")
            output.write("="*80 + "\n")
            
            if compact:
                output.write(format_compact(result))
            else:
                output.write(format_result(result))
            output.write("\n")
    
    if json_format:
        import json
        output.write(json.dumps([result_to_dict(r) for r in results], indent=2))
    
    return 0


def format_compact(result) -> str:
    """Format result in compact form (top condition only)"""
    lines = [
        f"Reactant A: {result.reactant_a_type} ({result.reactant_a_category})"
    ]
    
    if result.reactant_b_smiles:
        lines.append(f"Reactant B: {result.reactant_b_type} ({result.reactant_b_category})")
    
    lines.append(f"Predicted: {result.predicted_reaction_type} ({result.reaction_type_confidence*100:.0f}%)")
    lines.append(f"Matches: {result.total_matching_experiments} experiments")
    
    if result.recommendations:
        rec = result.recommendations[0]
        lines.extend([
            "",
            f"TOP RECOMMENDATION (Z-Score: {rec.avg_z_score:.2f}, Confidence: {rec.confidence_score:.1f}/100)",
            f"  Catalyst: {rec.catalyst}",
            f"  Ligand: {rec.ligand}",
            f"  Base: {rec.base}",
            f"  Solvent: {rec.solvent}",
            f"  Success: {rec.success_rate:.1f}% ({rec.num_experiments} exp, avg {rec.avg_yield:.1f}%)"
        ])
    else:
        lines.append("\nNo recommendations found")
    
    return "\n".join(lines)


def result_to_dict(result) -> dict:
    """Convert result to dictionary for JSON export"""
    return {
        'reactant_a': {
            'smiles': result.reactant_a_smiles,
            'type': result.reactant_a_type,
            'category': result.reactant_a_category
        },
        'reactant_b': {
            'smiles': result.reactant_b_smiles,
            'type': result.reactant_b_type,
            'category': result.reactant_b_category
        } if result.reactant_b_smiles else None,
        'predicted_reaction': {
            'type': result.predicted_reaction_type,
            'confidence': result.reaction_type_confidence
        },
        'database_match': {
            'experiments': result.total_matching_experiments,
            'coverage_pct': result.database_coverage
        },
        'recommendations': [
            {
                'rank': i + 1,
                'avg_z_score': rec.avg_z_score,
                'confidence_score': rec.confidence_score,
                'success_rate': rec.success_rate,
                'avg_yield': rec.avg_yield,
                'median_yield': rec.median_yield,
                'num_experiments': rec.num_experiments,
                'conditions': {
                    'catalyst': rec.catalyst,
                    'ligand': rec.ligand,
                    'base': rec.base,
                    'solvent': rec.solvent,
                    'secondary_solvent': rec.secondary_solvent,
                    'additive': rec.additive,
                    'coupling_reagent': rec.coupling_reagent
                },
                'metadata': {
                    'reaction_type': rec.reaction_type,
                    'reactant_types': rec.reactant_types,
                    'z_score_range': rec.z_score_range
                }
            }
            for i, rec in enumerate(result.recommendations)
        ]
    }


def result_to_json(result) -> str:
    """Convert result to JSON string"""
    import json
    return json.dumps(result_to_dict(result), indent=2)


if __name__ == '__main__':
    sys.exit(main())
