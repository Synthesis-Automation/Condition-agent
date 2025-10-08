"""
Demo: Dataset Analytics Module

Demonstrates all capabilities of the new dataset analytics module for
reaction condition analysis, ranking, and plate design.
"""

from chemtools import chem

def print_section(title):
    """Print a formatted section header."""
    print("\n" + "=" * 80)
    print(f" {title}")
    print("=" * 80 + "\n")


def demo_1_available_families():
    """Demo: List all available reaction families."""
    print_section("1. Available Reaction Families")
    
    families = chem.analytics.get_all_families()
    print(f"Found {len(families)} reaction families:\n")
    for family in families:
        print(f"  - {family}")


def demo_2_basic_stats():
    """Demo: Get basic dataset statistics."""
    print_section("2. Basic Dataset Statistics")
    
    for family in ["Suzuki", "Amide_formation"]:
        print(f"[{family}]")
        stats = chem.analytics.get_stats(family)
        
        print(f"  Total reactions: {stats['total_reactions']:,}")
        print(f"  Unique catalysts: {stats['unique_catalysts']}")
        print(f"  Unique ligands: (in catalytic system)")
        print(f"  Unique bases: {stats['unique_bases']}")
        print(f"  Unique solvents: {stats['unique_solvents']}")
        print(f"  Unique condition cores: {stats['unique_condition_cores']}")
        
        if stats['yield_stats']:
            ys = stats['yield_stats']
            print(f"\n  Yield data:")
            print(f"    Coverage: {ys['count']}/{stats['total_reactions']} reactions ({ys['count']/stats['total_reactions']*100:.1f}%)")
            print(f"    Range: {ys['min']:.1f}% - {ys['max']:.1f}%")
            print(f"    Mean: {ys['mean']:.1f}%")
            print(f"    Median: {ys['median']:.1f}%")
        
        if stats['temperature_stats']:
            ts = stats['temperature_stats']
            print(f"\n  Temperature data:")
            print(f"    Coverage: {ts['count']}/{stats['total_reactions']} reactions ({ts['count']/stats['total_reactions']*100:.1f}%)")
            print(f"    Range: {ts['min']:.0f}°C - {ts['max']:.0f}°C")
            print(f"    Mean: {ts['mean']:.0f}°C")
        
        print()


def demo_3_common_catalysts():
    """Demo: Get most common catalysts with yield data."""
    print_section("3. Common Catalysts (Ranked by Frequency)")
    
    print("[C_N_Coupling_Pd - Top 10 Pd Catalysts]\n")
    catalysts = chem.analytics.get_common_catalysts("C_N_Coupling_Pd", top_n=10)
    
    if catalysts:
        for name, count, avg_yield in catalysts:
            yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
            print(f"  {count:>5} reactions | Avg yield: {yield_str:>6} | {name}")
    else:
        print("  No catalyst data found")


def demo_4_common_ligands():
    """Demo: Get most common ligands with yield data."""
    print_section("4. Common Ligands (Ranked by Frequency)")
    
    print("[C_N_Coupling_Pd - Top 10 Ligands]\n")
    ligands = chem.analytics.get_common_ligands("C_N_Coupling_Pd", top_n=10)
    
    if ligands:
        for name, count, avg_yield in ligands:
            yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
            print(f"  {count:>5} reactions | Avg yield: {yield_str:>6} | {name}")
    else:
        print("  No ligand data found")


def demo_5_common_bases():
    """Demo: Get most common bases with yield data."""
    print_section("5. Common Bases (Ranked by Frequency)")
    
    print("[Suzuki - Top 10 Bases]\n")
    bases = chem.analytics.get_common_bases("Suzuki", top_n=10)
    
    if bases:
        for name, count, avg_yield in bases:
            yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
            print(f"  {count:>5} reactions | Avg yield: {yield_str:>6} | {name}")
    else:
        print("  No base data found")


def demo_6_common_solvents():
    """Demo: Get most common solvents with yield data."""
    print_section("6. Common Solvents (Ranked by Frequency)")
    
    print("[Suzuki - Top 10 Solvents]\n")
    solvents = chem.analytics.get_common_solvents("Suzuki", top_n=10)
    
    if solvents:
        for name, count, avg_yield in solvents:
            yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
            print(f"  {count:>5} reactions | Avg yield: {yield_str:>6} | {name}")
    else:
        print("  No solvent data found")


def demo_7_common_reagents():
    """Demo: Get most common reagents (all roles)."""
    print_section("7. Common Reagents (All Roles)")
    
    print("[Amide_formation - Top 10 Coupling Reagents]\n")
    reagents = chem.analytics.get_common_reagents("Amide_formation", role="COUPLING_REAGENT", top_n=10)
    
    if reagents:
        for name, role, count, avg_yield in reagents:
            yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
            print(f"  {count:>5} reactions | Avg yield: {yield_str:>6} | {name}")
    else:
        print("  No coupling reagent data found")
    
    print("\n[Amide_formation - Top 10 All Reagents]\n")
    all_reagents = chem.analytics.get_common_reagents("Amide_formation", top_n=10)
    
    if all_reagents:
        for name, role, count, avg_yield in all_reagents:
            yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
            print(f"  {count:>5} reactions | {role:<20} | Avg yield: {yield_str:>6} | {name}")
    else:
        print("  No reagent data found")


def demo_8_high_yield_filter():
    """Demo: Filter by minimum yield threshold."""
    print_section("8. High-Yield Filtering (min_yield >= 80%)")
    
    print("[Suzuki - High-Yield Catalysts (>= 80% avg yield)]\n")
    
    # Get all catalysts
    all_cats = chem.analytics.get_common_catalysts("Suzuki", top_n=20, min_yield=None)
    high_yield_cats = [c for c in all_cats if c[2] and c[2] >= 80.0]
    
    if high_yield_cats:
        print(f"Found {len(high_yield_cats)} catalysts with >= 80% avg yield:\n")
        for name, count, avg_yield in high_yield_cats[:10]:
            print(f"  {count:>5} reactions | Avg yield: {avg_yield:>6.1f}% | {name}")
    else:
        print("  No high-yield catalysts found")


def demo_9_condition_cores():
    """Demo: Get most common condition core combinations."""
    print_section("9. Common Condition Cores")
    
    print("[C_N_Coupling_Pd - Top 10 Condition Cores]\n")
    cores = chem.analytics.get_condition_cores("C_N_Coupling_Pd", top_n=10)
    
    if cores:
        for i, (core, count, avg_yield) in enumerate(cores, 1):
            yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
            core_display = core[:70] + "..." if len(core) > 70 else core
            print(f"  {i:>2}. {count:>4} reactions | Avg yield: {yield_str:>6} | {core_display}")
    else:
        print("  No condition core data found")


def demo_10_plate_recommendations():
    """Demo: Generate plate recommendations for HTE."""
    print_section("10. Plate Recommendations for High-Throughput Screening")
    
    print("[C_N_Coupling_Pd - 24-well plate, optimized for diversity]\n")
    
    plate = chem.analytics.get_plate_recommendations(
        family="C_N_Coupling_Pd",
        n_conditions=24,
        min_yield=60.0,
        optimize_for='diversity'
    )
    
    print(f"Generated {len(plate)} conditions:\n")
    
    for i, cond in enumerate(plate[:10], 1):  # Show first 10
        cat = cond['catalyst'] or 'None'
        lig = cond['ligand'] or 'None'
        base = cond['base'] or 'None'
        solv = cond['solvent'] or 'None'
        temp = f"{cond['temperature_c']:.0f}°C" if cond['temperature_c'] else 'N/A'
        time = f"{cond['time_h']:.1f}h" if cond['time_h'] else 'N/A'
        yield_str = f"{cond['avg_yield']:.1f}%" if cond['avg_yield'] else 'N/A'
        
        print(f"  {i:>2}. [{cond['condition_id']}]")
        print(f"      Catalyst: {cat[:40]}")
        print(f"      Ligand: {lig[:40]}")
        print(f"      Base: {base[:40]}")
        print(f"      Solvent: {solv[:40]}")
        print(f"      Conditions: {temp}, {time}")
        print(f"      Performance: {yield_str} avg yield, {cond['frequency']} precedents")
        print(f"      Score: {cond['score']:.1f}")
        print()
    
    if len(plate) > 10:
        print(f"  ... and {len(plate) - 10} more conditions")


def demo_11_full_summary():
    """Demo: Print comprehensive analytics summary."""
    print_section("11. Comprehensive Analytics Summary")
    
    print("Using print_summary() for quick overview:\n")
    chem.analytics.print_summary("Suzuki", top_n=5)


def demo_12_comparison():
    """Demo: Compare datasets."""
    print_section("12. Dataset Comparison")
    
    families = ["Suzuki", "C_N_Coupling_Pd", "C_N_Coupling_Cu", "Amide_formation"]
    
    print(f"{'Family':<25} {'Reactions':<12} {'Avg Yield':<12} {'Avg Temp':<12}")
    print("-" * 65)
    
    for family in families:
        try:
            stats = chem.analytics.get_stats(family)
            total = stats['total_reactions']
            
            if stats['yield_stats']:
                avg_yield = f"{stats['yield_stats']['mean']:.1f}%"
            else:
                avg_yield = "N/A"
            
            if stats['temperature_stats']:
                avg_temp = f"{stats['temperature_stats']['mean']:.0f}°C"
            else:
                avg_temp = "N/A"
            
            print(f"{family:<25} {total:<12,} {avg_yield:<12} {avg_temp:<12}")
        except FileNotFoundError:
            print(f"{family:<25} {'Not found':<12} {'N/A':<12} {'N/A':<12}")


if __name__ == "__main__":
    print("=" * 80)
    print(" CHEMTOOLS DATASET ANALYTICS - COMPREHENSIVE DEMO")
    print("=" * 80)
    print()
    print("This demo showcases all capabilities of the new analytics module:")
    print("  - Dataset statistics and distributions")
    print("  - Ranked reagent/catalyst/ligand/base/solvent lists")
    print("  - Yield-based filtering")
    print("  - High-throughput plate design recommendations")
    print()
    
    try:
        demo_1_available_families()
        demo_2_basic_stats()
        demo_3_common_catalysts()
        demo_4_common_ligands()
        demo_5_common_bases()
        demo_6_common_solvents()
        demo_7_common_reagents()
        demo_8_high_yield_filter()
        demo_9_condition_cores()
        demo_10_plate_recommendations()
        demo_11_full_summary()
        demo_12_comparison()
        
        print_section("DEMO COMPLETE")
        print("All analytics features demonstrated successfully!")
        print()
        print("Key capabilities:")
        print("  ✓ Dataset statistics and metadata")
        print("  ✓ Frequency rankings with yield data")
        print("  ✓ High-yield filtering")
        print("  ✓ Condition core analysis")
        print("  ✓ Plate design recommendations")
        print("  ✓ Multi-dataset comparison")
        print()
        print("Use these analytics to:")
        print("  → Inform condition recommendations")
        print("  → Design HTE experiments")
        print("  → Understand dataset coverage")
        print("  → Identify successful patterns")
        print()
        
    except Exception as e:
        print(f"\n[ERROR] {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
