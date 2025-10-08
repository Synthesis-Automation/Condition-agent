#!/usr/bin/env python
"""
Demo: Catalytic System Analysis
Shows how to analyze complete catalytic systems (catalyst + ligand combinations)
rather than individual components.
"""
import sys
from pathlib import Path

# Add parent directory to path
parent_dir = Path(__file__).parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from chemtools import chem


def demo_catalytic_systems():
    """Demonstrate catalytic system analysis vs individual component analysis."""
    
    print("=" * 80)
    print("DEMO: Catalytic System Analysis")
    print("=" * 80)
    print()
    print("In real chemistry, catalytic systems work as complete units:")
    print("e.g., Pd(OAc)₂ + RuPhos is different from Pd(OAc)₂ + XPhos")
    print()
    print("This demo shows how to analyze complete catalytic systems with")
    print("the new get_common_catalytic_systems() function.")
    print()
    
    family = "C_N_Coupling_Pd"
    
    # ========================================================================
    # Part 1: Individual Component Analysis (Old Way)
    # ========================================================================
    print("─" * 80)
    print("PART 1: Individual Component Analysis (Old Way)")
    print("─" * 80)
    print()
    
    catalysts = chem.analytics.get_common_catalysts(family, top_n=5)
    print("Top 5 Individual Catalysts:")
    for i, (name, count, avg_yield) in enumerate(catalysts, 1):
        yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
        print(f"  {i}. {name}")
        print(f"     {count} reactions, avg yield: {yield_str}")
    
    print()
    ligands = chem.analytics.get_common_ligands(family, top_n=5)
    if ligands:
        print("Top 5 Individual Ligands:")
        for i, (name, count, avg_yield) in enumerate(ligands, 1):
            yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
            print(f"  {i}. {name}")
            print(f"     {count} reactions, avg yield: {yield_str}")
    else:
        print("⚠️  No separate ligand data found in this dataset.")
        print("   (Ligands are stored in catalytic_system with catalysts)")
    
    print()
    print("❌ Problem: We lose the important catalyst-ligand pairing!")
    print("   We can't see which ligands work best with which catalysts.")
    
    # ========================================================================
    # Part 2: Complete System Analysis (New Way)
    # ========================================================================
    print()
    print("─" * 80)
    print("PART 2: Complete Catalytic System Analysis (New Way)")
    print("─" * 80)
    print()
    
    systems = chem.analytics.get_common_catalytic_systems(family, top_n=10)
    print("Top 10 Catalytic Systems (Catalyst + Ligand Combinations):")
    print()
    for i, (system_str, count, avg_yield) in enumerate(systems, 1):
        yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
        print(f"  {i}. {system_str}")
        print(f"     {count} reactions, avg yield: {yield_str}")
        print()
    
    print("✅ Advantage: We now see complete catalytic systems!")
    print("   - Which catalyst/ligand pairs are most common?")
    print("   - What's the success rate of each complete system?")
    print("   - Can identify synergistic combinations")
    
    # ========================================================================
    # Part 3: High-Yield Systems
    # ========================================================================
    print()
    print("─" * 80)
    print("PART 3: High-Yield Catalytic Systems (≥85% average yield)")
    print("─" * 80)
    print()
    
    high_yield_systems = chem.analytics.get_common_catalytic_systems(
        family, top_n=20, min_yield=85.0
    )
    
    if high_yield_systems:
        print(f"Found {len(high_yield_systems)} high-yield catalytic systems:")
        print()
        for i, (system_str, count, avg_yield) in enumerate(high_yield_systems, 1):
            print(f"  {i}. {system_str}")
            print(f"     {count} reactions, avg yield: {avg_yield:.1f}%")
            print()
    else:
        print("⚠️  No catalytic systems meet the ≥85% yield threshold")
        print("   (Try lowering the min_yield parameter)")
    
    # ========================================================================
    # Part 4: Comparison with Suzuki (larger dataset)
    # ========================================================================
    print()
    print("─" * 80)
    print("PART 4: Suzuki Coupling - Catalytic Systems")
    print("─" * 80)
    print()
    
    suzuki_systems = chem.analytics.get_common_catalytic_systems("Suzuki", top_n=10)
    print("Top 10 Catalytic Systems for Suzuki Coupling:")
    print()
    for i, (system_str, count, avg_yield) in enumerate(suzuki_systems, 1):
        yield_str = f"{avg_yield:.1f}%" if avg_yield else "N/A"
        # Truncate very long names
        display_str = system_str[:75] + "..." if len(system_str) > 75 else system_str
        print(f"  {i}. {display_str}")
        print(f"     {count} reactions, avg yield: {yield_str}")
        print()
    
    # ========================================================================
    # Summary
    # ========================================================================
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print()
    print("The new get_common_catalytic_systems() function:")
    print("  ✅ Preserves catalyst-ligand relationships")
    print("  ✅ Shows complete systems as they're used in practice")
    print("  ✅ Helps identify successful combinations")
    print("  ✅ Supports yield filtering for high-performance systems")
    print("  ✅ Useful for plate design and condition recommendations")
    print()
    print("Use this for:")
    print("  • Understanding which catalyst/ligand pairs work best")
    print("  • Designing HTE plates with proven combinations")
    print("  • Identifying gaps in your experimental coverage")
    print("  • Learning from successful precedents")
    print()


if __name__ == "__main__":
    demo_catalytic_systems()
