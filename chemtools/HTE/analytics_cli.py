"""
CLI for HTE Analytics Tools

Provides command-line access to HTE database analytics:
- List reactant pairs
- Analyze catalysts
- View reaction type summaries
- Export filtered datasets
"""

import argparse
import sys
from pathlib import Path
from chemtools.HTE.analytics import HTEAnalytics


def cmd_list_pairs(args):
    """List reactant pairs command"""
    analytics = HTEAnalytics(args.db_path)
    
    print("\n" + "="*80)
    print("📋 REACTANT PAIR ANALYSIS")
    print("="*80)
    
    if args.reaction:
        print(f"Reaction Type: {args.reaction}")
    if args.catalyst:
        print(f"Catalyst Filter: {args.catalyst}")
    if args.min_experiments > 1:
        print(f"Min Experiments: {args.min_experiments}")
    print()
    
    df = analytics.list_reactant_pairs(
        reaction_type=args.reaction,
        catalyst_filter=args.catalyst,
        min_experiments=args.min_experiments,
        sort_by=args.sort
    )
    
    if len(df) == 0:
        print("❌ No matching reactant pairs found")
        return
    
    print(f"Found {len(df)} reactant pair combinations\n")
    
    # Display results
    if args.compact:
        # Compact format
        for i, row in df.head(args.top).iterrows():
            print(f"{i+1}. {row['Reactant_A_Type']} + {row['Reactant_B_Type']}")
            print(f"   Reaction: {row['Reaction_Type']}")
            print(f"   Experiments: {row['Num_Experiments']}, "
                  f"Avg Yield: {row['Avg_Yield']:.1f}%, "
                  f"Success Rate: {row['Success_Rate']:.1f}%")
            print(f"   Top Catalyst: {row['Top_Catalyst']}")
            print()
    else:
        # Full table format
        pd_options = {
            'display.max_rows': args.top,
            'display.max_columns': None,
            'display.width': None,
            'display.max_colwidth': 50
        }
        
        import pandas as pd
        with pd.option_context(*[item for pair in pd_options.items() for item in pair]):
            print(df.head(args.top).to_string(index=False))
    
    if args.output:
        df.to_csv(args.output, index=False)
        print(f"\n✅ Saved results to {args.output}")


def cmd_catalysts(args):
    """Analyze catalysts command"""
    analytics = HTEAnalytics(args.db_path)
    
    print("\n" + "="*80)
    print("⚗️  CATALYST ANALYSIS")
    print("="*80)
    
    if args.reaction:
        print(f"Reaction Type: {args.reaction}")
    if args.reactant_a:
        print(f"Reactant A Type: {args.reactant_a}")
    if args.reactant_b:
        print(f"Reactant B Type: {args.reactant_b}")
    print()
    
    df = analytics.get_catalyst_stats(
        reaction_type=args.reaction,
        reactant_a_type=args.reactant_a,
        reactant_b_type=args.reactant_b
    )
    
    if len(df) == 0:
        print("❌ No catalysts found matching criteria")
        return
    
    print(f"Found {len(df)} catalysts\n")
    
    # Display results
    if args.compact:
        for i, row in df.head(args.top).iterrows():
            print(f"{i+1}. {row['Catalyst']} [{row['Metal']}]")
            print(f"   Experiments: {row['Num_Experiments']}, "
                  f"Avg Yield: {row['Avg_Yield']:.1f}%, "
                  f"Success Rate: {row['Success_Rate']:.1f}%")
            if pd.notna(row['Reaction_Types']) and len(row['Reaction_Types']) < 100:
                print(f"   Reactions: {row['Reaction_Types']}")
            print()
    else:
        import pandas as pd
        pd_options = {
            'display.max_rows': args.top,
            'display.max_columns': None,
            'display.width': None,
            'display.max_colwidth': 60
        }
        with pd.option_context(*[item for pair in pd_options.items() for item in pair]):
            print(df.head(args.top).to_string(index=False))
    
    if args.output:
        df.to_csv(args.output, index=False)
        print(f"\n✅ Saved results to {args.output}")


def cmd_reactions(args):
    """Analyze reaction types command"""
    analytics = HTEAnalytics(args.db_path)
    
    print("\n" + "="*80)
    print("🧪 REACTION TYPE SUMMARY")
    print("="*80)
    print()
    
    df = analytics.get_reaction_type_summary()
    
    print(f"Found {len(df)} reaction types in database\n")
    
    # Display results
    if args.compact:
        for i, row in df.head(args.top).iterrows():
            print(f"{i+1}. {row['Reaction_Type']}")
            print(f"   Experiments: {row['Num_Experiments']:,}, "
                  f"Pairs: {row['Num_Reactant_Pairs']}, "
                  f"Catalysts: {row['Num_Catalysts']}")
            print(f"   Avg Yield: {row['Avg_Yield']:.1f}%, "
                  f"Success Rate: {row['Success_Rate']:.1f}%")
            print(f"   Top Catalyst: {row['Top_Catalyst']}")
            print(f"   Top Pair: {row['Top_Reactant_Pair']}")
            print()
    else:
        import pandas as pd
        pd_options = {
            'display.max_rows': args.top,
            'display.max_columns': None,
            'display.width': None,
            'display.max_colwidth': 40
        }
        with pd.option_context(*[item for pair in pd_options.items() for item in pair]):
            print(df.head(args.top).to_string(index=False))
    
    if args.output:
        df.to_csv(args.output, index=False)
        print(f"\n✅ Saved results to {args.output}")


def cmd_metals(args):
    """Analyze metal usage command"""
    analytics = HTEAnalytics(args.db_path)
    
    print("\n" + "="*80)
    print("🔬 METAL USAGE ANALYSIS")
    print("="*80)
    print()
    
    result = analytics.analyze_metal_usage()
    
    print(f"Total Experiments: {result['total_experiments']:,}\n")
    
    print("Metal Distribution:")
    print("-" * 50)
    
    import pandas as pd
    df = result['metal_distribution']
    
    for _, row in df.iterrows():
        metal = row['Metal']
        count = row['Num_Experiments']
        pct = row['Percentage']
        bar = '█' * int(pct / 2)
        print(f"{metal:>4}: {bar:<35} {count:>6,} ({pct:>5.1f}%)")
    
    if args.detailed:
        print("\n\nMetal Usage by Reaction Type:")
        print("-" * 50)
        
        for metal, reactions in sorted(result['by_reaction_type'].items(), 
                                      key=lambda x: sum(x[1].values()), reverse=True):
            if metal and metal != 'Other':
                print(f"\n{metal}:")
                for rxn, count in sorted(reactions.items(), key=lambda x: x[1], reverse=True)[:5]:
                    print(f"  {rxn}: {count:,}")
    
    if args.output:
        df.to_csv(args.output, index=False)
        print(f"\n✅ Saved metal distribution to {args.output}")


def cmd_export(args):
    """Export filtered dataset command"""
    analytics = HTEAnalytics(args.db_path)
    
    print("\n" + "="*80)
    print("💾 EXPORT FILTERED DATASET")
    print("="*80)
    
    if args.reaction:
        print(f"Reaction Type: {args.reaction}")
    if args.catalyst:
        print(f"Catalyst Filter: {args.catalyst}")
    if args.reactant_a:
        print(f"Reactant A Type: {args.reactant_a}")
    if args.reactant_b:
        print(f"Reactant B Type: {args.reactant_b}")
    if args.min_yield:
        print(f"Min Yield: {args.min_yield}%")
    print()
    
    count = analytics.export_subset(
        output_path=args.output,
        reaction_type=args.reaction,
        catalyst_filter=args.catalyst,
        reactant_a_type=args.reactant_a,
        reactant_b_type=args.reactant_b,
        min_yield=args.min_yield
    )
    
    print(f"\n✅ Export complete: {count:,} experiments")


def main():
    parser = argparse.ArgumentParser(
        description="HTE Database Analytics Tools",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # List all Suzuki reactant pairs with Pd catalysts
  python -m chemtools.HTE.analytics pairs --reaction Suzuki --catalyst Pd --top 10
  
  # Analyze Cu catalysts in C-N coupling
  python -m chemtools.HTE.analytics catalysts --reaction "C-N" --catalyst Cu --compact
  
  # View reaction type summary
  python -m chemtools.HTE.analytics reactions --top 20
  
  # Analyze metal usage
  python -m chemtools.HTE.analytics metals --detailed
  
  # Export Pd-catalyzed Suzuki data
  python -m chemtools.HTE.analytics export suzuki_pd.csv --reaction Suzuki --catalyst Pd
        """
    )
    
    parser.add_argument(
        '--db-path',
        default='data/HTE_db',
        help='Path to HTE database CSV/JSONL or directory (default: data/HTE_db)'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Analytics command')
    
    # List pairs command
    pairs_parser = subparsers.add_parser('pairs', help='List reactant pairs')
    pairs_parser.add_argument('--reaction', help='Filter by reaction type')
    pairs_parser.add_argument('--catalyst', help='Filter by catalyst metal (e.g., Pd, Cu)')
    pairs_parser.add_argument('--min-experiments', type=int, default=1,
                             help='Minimum number of experiments (default: 1)')
    pairs_parser.add_argument('--sort', choices=['count', 'success_rate'], default='count',
                             help='Sort by count or success_rate (default: count)')
    pairs_parser.add_argument('--top', type=int, default=20,
                             help='Number of results to show (default: 20)')
    pairs_parser.add_argument('--compact', action='store_true',
                             help='Use compact output format')
    pairs_parser.add_argument('-o', '--output', help='Save results to CSV')
    
    # Catalysts command
    cat_parser = subparsers.add_parser('catalysts', help='Analyze catalysts')
    cat_parser.add_argument('--reaction', help='Filter by reaction type')
    cat_parser.add_argument('--reactant-a', help='Filter by reactant A type')
    cat_parser.add_argument('--reactant-b', help='Filter by reactant B type')
    cat_parser.add_argument('--top', type=int, default=20,
                           help='Number of results to show (default: 20)')
    cat_parser.add_argument('--compact', action='store_true',
                           help='Use compact output format')
    cat_parser.add_argument('-o', '--output', help='Save results to CSV')
    
    # Reactions command
    rxn_parser = subparsers.add_parser('reactions', help='Analyze reaction types')
    rxn_parser.add_argument('--top', type=int, default=20,
                           help='Number of results to show (default: 20)')
    rxn_parser.add_argument('--compact', action='store_true',
                           help='Use compact output format')
    rxn_parser.add_argument('-o', '--output', help='Save results to CSV')
    
    # Metals command
    metals_parser = subparsers.add_parser('metals', help='Analyze metal usage')
    metals_parser.add_argument('--detailed', action='store_true',
                              help='Show detailed breakdown by reaction type')
    metals_parser.add_argument('-o', '--output', help='Save results to CSV')
    
    # Export command
    export_parser = subparsers.add_parser('export', help='Export filtered dataset')
    export_parser.add_argument('output', help='Output CSV file path')
    export_parser.add_argument('--reaction', help='Filter by reaction type')
    export_parser.add_argument('--catalyst', help='Filter by catalyst metal')
    export_parser.add_argument('--reactant-a', help='Filter by reactant A type')
    export_parser.add_argument('--reactant-b', help='Filter by reactant B type')
    export_parser.add_argument('--min-yield', type=float, help='Minimum yield threshold')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    try:
        if args.command == 'pairs':
            cmd_list_pairs(args)
        elif args.command == 'catalysts':
            cmd_catalysts(args)
        elif args.command == 'reactions':
            cmd_reactions(args)
        elif args.command == 'metals':
            cmd_metals(args)
        elif args.command == 'export':
            cmd_export(args)
        return 0
    except Exception as e:
        print(f"\n❌ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
