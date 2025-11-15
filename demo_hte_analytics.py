"""
Demo: HTE Analytics Tools

Demonstrates the new analytics capabilities for exploring the HTE database.
"""

from chemtools.HTE import HTEAnalytics


def demo_reactant_pairs():
    """Demo: List reactant pairs by reaction type and catalyst"""
    print("\n" + "="*80)
    print("DEMO 1: Reactant Pair Analysis")
    print("="*80)
    
    analytics = HTEAnalytics()
    
    # Example 1: All Suzuki pairs with Pd
    print("\n1️⃣  Suzuki Coupling - Palladium Catalysts")
    print("-" * 80)
    df = analytics.list_reactant_pairs(
        reaction_type="Suzuki",
        catalyst_filter="Pd",
        min_experiments=50,
        sort_by="count"
    )
    
    print(f"Found {len(df)} reactant pairs (≥50 experiments)\n")
    for i, row in df.head(5).iterrows():
        print(f"{i+1}. {row['Reactant_A_Type']} + {row['Reactant_B_Type']}")
        print(f"   Experiments: {row['Num_Experiments']}, "
              f"Success: {row['Success_Rate']:.1f}%, "
              f"Avg Yield: {row['Avg_Yield']:.1f}%")
        print(f"   Top Catalyst: {row['Top_Catalyst']}")
    
    # Example 2: C-N Coupling with Cu
    print("\n2️⃣  C-N Coupling - Copper Catalysts")
    print("-" * 80)
    df = analytics.list_reactant_pairs(
        reaction_type="C-N",
        catalyst_filter="Cu",
        min_experiments=20,
        sort_by="count"
    )
    
    print(f"Found {len(df)} reactant pairs (≥20 experiments)\n")
    for i, row in df.head(5).iterrows():
        print(f"{i+1}. {row['Reactant_A_Type']} + {row['Reactant_B_Type']}")
        print(f"   Experiments: {row['Num_Experiments']}, "
              f"Success: {row['Success_Rate']:.1f}%, "
              f"Avg Yield: {row['Avg_Yield']:.1f}%")
        print(f"   Top Catalyst: {row['Top_Catalyst']}")


def demo_catalyst_analysis():
    """Demo: Analyze catalyst statistics"""
    print("\n" + "="*80)
    print("DEMO 2: Catalyst Analysis")
    print("="*80)
    
    analytics = HTEAnalytics()
    
    # Example 1: Top Pd catalysts for Suzuki
    print("\n1️⃣  Top Palladium Catalysts for Suzuki Coupling")
    print("-" * 80)
    df = analytics.get_catalyst_stats(reaction_type="Suzuki")
    df = df[df['Metal'] == 'Pd']
    
    print(f"Found {len(df)} Pd catalysts for Suzuki\n")
    for i, row in df.head(10).iterrows():
        print(f"{i+1}. {row['Catalyst']}")
        print(f"   Experiments: {row['Num_Experiments']}, "
              f"Avg Yield: {row['Avg_Yield']:.1f}%, "
              f"Success: {row['Success_Rate']:.1f}%")
    
    # Example 2: Top Cu catalysts overall
    print("\n2️⃣  Top Copper Catalysts (All Reactions)")
    print("-" * 80)
    df = analytics.get_catalyst_stats()
    df = df[df['Metal'] == 'Cu']
    
    print(f"Found {len(df)} Cu catalysts\n")
    for i, row in df.head(10).iterrows():
        print(f"{i+1}. {row['Catalyst']}")
        print(f"   Experiments: {row['Num_Experiments']}, "
              f"Avg Yield: {row['Avg_Yield']:.1f}%")
        # Show first few reaction types
        rxns = row['Reaction_Types'].split(', ')[:3]
        print(f"   Reactions: {', '.join(rxns)}")


def demo_reaction_summary():
    """Demo: Reaction type summary"""
    print("\n" + "="*80)
    print("DEMO 3: Reaction Type Summary")
    print("="*80)
    
    analytics = HTEAnalytics()
    
    df = analytics.get_reaction_type_summary()
    
    print(f"\nFound {len(df)} reaction types in database\n")
    
    for i, row in df.head(15).iterrows():
        print(f"{i+1}. {row['Reaction_Type']}")
        print(f"   Experiments: {row['Num_Experiments']:,}, "
              f"Pairs: {row['Num_Reactant_Pairs']}, "
              f"Catalysts: {row['Num_Catalysts']}")
        print(f"   Performance: {row['Success_Rate']:.1f}% success, {row['Avg_Yield']:.1f}% avg")
        print(f"   Top: {row['Top_Catalyst']}")
        print(f"   Most Common: {row['Top_Reactant_Pair']}")
        print()


def demo_metal_analysis():
    """Demo: Metal usage analysis"""
    print("\n" + "="*80)
    print("DEMO 4: Metal Usage Analysis")
    print("="*80)
    
    analytics = HTEAnalytics()
    
    result = analytics.analyze_metal_usage()
    
    print(f"\nTotal Experiments: {result['total_experiments']:,}\n")
    print("Metal Distribution:")
    print("-" * 60)
    
    df = result['metal_distribution']
    
    for _, row in df.iterrows():
        metal = row['Metal']
        count = row['Num_Experiments']
        pct = row['Percentage']
        bar = '█' * int(pct / 2)
        print(f"{metal:>5}: {bar:<35} {count:>6,} ({pct:>5.1f}%)")
    
    print("\n\nTop Reaction Types by Metal:")
    print("-" * 60)
    
    # Show top 3 reactions for each major metal
    for metal in ['Pd', 'Cu', 'Ni']:
        if metal in result['by_reaction_type']:
            reactions = result['by_reaction_type'][metal]
            top_3 = sorted(reactions.items(), key=lambda x: x[1], reverse=True)[:3]
            print(f"\n{metal}:")
            for rxn, count in top_3:
                print(f"  {rxn}: {count:,} experiments")


def demo_similar_pairs():
    """Demo: Find similar reactant pairs"""
    print("\n" + "="*80)
    print("DEMO 5: Find Similar Reactant Pairs")
    print("="*80)
    
    analytics = HTEAnalytics()
    
    # Query: ArCl + ArB(OH)2 (typical Suzuki)
    print("\n1️⃣  Query: ArCl + ArB(OH)2 (Suzuki)")
    print("-" * 80)
    print("Finding pairs with same reaction type...\n")
    
    df = analytics.find_similar_pairs(
        reactant_a_type="ArCl",
        reactant_b_type="ArB(OH)2",
        similarity_criteria="reaction_type"
    )
    
    if len(df) > 0:
        print(f"Found {len(df)} similar pairs\n")
        for i, row in df.head(5).iterrows():
            print(f"{i+1}. {row['Reactant_A_Type']} + {row['Reactant_B_Type']}")
            print(f"   Reaction: {row['Reaction_Type']}")
            print(f"   Experiments: {row['Num_Experiments']}, "
                  f"Avg Yield: {row['Avg_Yield']:.1f}%")
    else:
        print("No similar pairs found")
    
    # Query with catalyst similarity
    print("\n2️⃣  Query: ArBr + RNH2-a-branch")
    print("-" * 80)
    print("Finding pairs with same catalyst metal...\n")
    
    df = analytics.find_similar_pairs(
        reactant_a_type="ArBr",
        reactant_b_type="RNH2-a-branch",
        similarity_criteria="catalyst"
    )
    
    if len(df) > 0:
        print(f"Found {len(df)} similar pairs\n")
        for i, row in df.head(5).iterrows():
            print(f"{i+1}. {row['Reactant_A_Type']} + {row['Reactant_B_Type']}")
            print(f"   Reaction: {row['Reaction_Type']}")
            print(f"   Experiments: {row['Num_Experiments']}, "
                  f"Success: {row['Success_Rate']:.1f}%")


def demo_export():
    """Demo: Export filtered subset"""
    print("\n" + "="*80)
    print("DEMO 6: Export Filtered Datasets")
    print("="*80)
    
    analytics = HTEAnalytics()
    
    # Export high-performing Suzuki data
    print("\n1️⃣  Exporting high-performing Suzuki reactions (≥80% yield)...")
    count = analytics.export_subset(
        output_path="suzuki_high_yield.csv",
        reaction_type="Suzuki",
        min_yield=80.0
    )
    print(f"   → Exported {count:,} experiments")
    
    # Export Cu-catalyzed C-N coupling
    print("\n2️⃣  Exporting Cu-catalyzed C-N coupling data...")
    count = analytics.export_subset(
        output_path="cn_copper.csv",
        reaction_type="C-N",
        catalyst_filter="Cu"
    )
    print(f"   → Exported {count:,} experiments")


if __name__ == "__main__":
    print("\n" + "="*80)
    print("HTE ANALYTICS TOOLS - COMPREHENSIVE DEMO")
    print("="*80)
    
    try:
        demo_reactant_pairs()
        demo_catalyst_analysis()
        demo_reaction_summary()
        demo_metal_analysis()
        demo_similar_pairs()
        demo_export()
        
        print("\n" + "="*80)
        print("✅ DEMO COMPLETE")
        print("="*80)
        print("\nFor CLI usage, try:")
        print("  python -m chemtools.HTE.analytics_cli pairs --reaction Suzuki --catalyst Pd")
        print("  python -m chemtools.HTE.analytics_cli catalysts --reaction 'C-N' --compact")
        print("  python -m chemtools.HTE.analytics_cli reactions --top 20")
        print("  python -m chemtools.HTE.analytics_cli metals --detailed")
        print()
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
