"""
Demo script for reagent database analytics.

Shows how to get statistics about ligands, bases, and other reagents.

Usage:
    python scripts/demo_reagent_analytics.py
"""

import sys
from pathlib import Path

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools.reagent import (
    get_database_statistics,
    get_family_statistics,
    find_reagents_by_family,
    count_reagents_by_type,
    get_all_reagent_types,
    print_database_summary,
    print_role_summary,
)


def demo_basic_counts():
    """Demo: Basic counting functions."""
    print("="*70)
    print("BASIC COUNTS")
    print("="*70)
    
    # Get all available types
    types = get_all_reagent_types()
    print(f"\n✅ Available reagent types: {', '.join(types)}")
    
    # Count each type
    print("\nReagent counts by type:")
    for reagent_type in types:
        count = count_reagents_by_type(reagent_type)
        print(f"  {reagent_type:25s}: {count:4d} reagents")


def demo_statistics():
    """Demo: Detailed statistics."""
    print("\n" + "="*70)
    print("DETAILED STATISTICS")
    print("="*70)
    
    stats = get_database_statistics()
    
    print(f"\nTotal reagents: {stats['total_reagents']}")
    print(f"  With CAS numbers: {stats['total_with_cas']} ({stats.get('percent_with_cas', 0):.1f}%)")
    print(f"  With InChIKeys: {stats['total_with_inchikey']} ({stats.get('percent_with_inchikey', 0):.1f}%)")
    print(f"  With SMILES: {stats['total_with_smiles']} ({stats.get('percent_with_smiles', 0):.1f}%)")
    
    print(f"\nTop 5 families:")
    for family, count in stats['top_families'][:5]:
        print(f"  {family:40s}: {count:3d} reagents")


def demo_family_search():
    """Demo: Find reagents by family."""
    print("\n" + "="*70)
    print("FAMILY SEARCH")
    print("="*70)
    
    # Get ligand family statistics
    ligand_stats = get_family_statistics('ligand')
    
    print(f"\nTotal ligands: {ligand_stats['total_reagents']}")
    print(f"Total ligand families: {ligand_stats['total_families']}")
    
    # Find all phosphines
    print("\n✅ Finding all trialkyl/triaryl phosphines...")
    phosphines = find_reagents_by_family('ligand', 'trialkyl_triaryl_phosphines')
    
    print(f"\nFound {len(phosphines)} phosphine ligands:")
    for ligand in phosphines[:10]:  # Show first 10
        name = ligand.get('name', 'Unknown')
        cas = ligand.get('cas', 'N/A')
        abbr = ligand.get('abbreviation', [])
        abbr_str = f" ({abbr[0]})" if abbr else ""
        print(f"  - {name}{abbr_str}")
        print(f"    CAS: {cas}")
    
    if len(phosphines) > 10:
        print(f"  ... and {len(phosphines) - 10} more")


def demo_base_analysis():
    """Demo: Analyze base families."""
    print("\n" + "="*70)
    print("BASE FAMILY ANALYSIS")
    print("="*70)
    
    base_stats = get_family_statistics('base')
    
    print(f"\nTotal bases: {base_stats['total_reagents']}")
    print(f"Total base families: {base_stats['total_families']}")
    
    print("\nBase families:")
    for family_data in base_stats['families'][:10]:  # Top 10
        name = family_data['name']
        count = family_data['count']
        percent = (count / base_stats['total_reagents'] * 100)
        print(f"  {name:35s}: {count:2d} ({percent:4.1f}%)")


def demo_pretty_summaries():
    """Demo: Pretty-printed summaries."""
    print("\n" + "="*70)
    print("PRETTY SUMMARIES")
    print("="*70)
    
    print("\n[1] Database Summary:")
    print_database_summary()
    
    print("\n[2] Ligand Summary:")
    print_role_summary('ligand')


def show_api_reference():
    """Show API reference."""
    print("\n" + "="*70)
    print("API REFERENCE")
    print("="*70)
    print("""
✅ Basic Counts:
   from chemtools.reagent import (
       get_all_reagent_types,
       count_reagents_by_type,
       get_all_reagents_by_type
   )
   
   types = get_all_reagent_types()
   count = count_reagents_by_type('ligand')
   all_ligands = get_all_reagents_by_type('ligand')

✅ Statistics:
   from chemtools.reagent import get_database_statistics
   
   stats = get_database_statistics()
   # Returns: total_reagents, by_type, families, etc.

✅ Family Analysis:
   from chemtools.reagent import (
       get_family_statistics,
       find_reagents_by_family
   )
   
   stats = get_family_statistics('ligand')
   phosphines = find_reagents_by_family('ligand', 'trialkyl_triaryl_phosphines')

✅ Pretty Print:
   from chemtools.reagent import print_database_summary, print_role_summary
   
   print_database_summary()
   print_role_summary('ligand')

✅ Missing Data:
   from chemtools.reagent import (
       get_missing_data_report,
       print_missing_data_report
   )
   
   report = get_missing_data_report()
   print_missing_data_report()
""")


def main():
    """Run all demos."""
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    
    print("="*70)
    print("REAGENT DATABASE ANALYTICS DEMO")
    print("="*70)
    
    demos = [
        demo_basic_counts,
        demo_statistics,
        demo_family_search,
        demo_base_analysis,
    ]
    
    for demo in demos:
        try:
            demo()
        except Exception as e:
            print(f"\n❌ Error in {demo.__name__}: {e}")
            import traceback
            traceback.print_exc()
    
    show_api_reference()
    
    print("\n" + "="*70)
    print("DEMO COMPLETE")
    print("="*70)
    print("\n💡 Try the CLI:")
    print("   python -m chemtools.reagent.analytics summary")
    print("   python -m chemtools.reagent.analytics role ligand")
    print("   python -m chemtools.reagent.analytics missing\n")


if __name__ == "__main__":
    main()
