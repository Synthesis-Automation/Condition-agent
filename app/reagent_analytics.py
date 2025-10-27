#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reagent Database Analytics and Summary Tool
============================================

This tool demonstrates the reagent module's capabilities by:
1. Providing comprehensive database statistics
2. Demonstrating search by CAS number
3. Demonstrating search by name/abbreviation
4. Showing role-based filtering

Similar to the dataset analysis app, this showcases the reagent module's
search and analytics functionality.

Usage:
    python app/reagent_analytics.py                    # Full analytics + demos
    python app/reagent_analytics.py --search "NaOtBu"  # Search by name
    python app/reagent_analytics.py --cas "865-48-5"   # Search by CAS
    python app/reagent_analytics.py --role base        # Show all bases
"""

import sys
from pathlib import Path
from typing import Dict, List, Any, Optional
import argparse

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools import reagent


def print_section(title: str, char: str = "="):
    """Print a formatted section header."""
    width = 80
    print()
    print(char * width)
    print(title.center(width))
    print(char * width)
    print()


def print_subsection(title: str):
    """Print a formatted subsection header."""
    print()
    print("-" * 80)
    print(f"  {title}")
    print("-" * 80)


def format_reagent(reagent_data: Dict[str, Any], compact: bool = False) -> str:
    """Format reagent information for display."""
    name = reagent_data.get('name', 'Unknown')
    cas = reagent_data.get('cas', 'No CAS')
    abbrev = reagent_data.get('abbreviation', [])
    abbrev_str = f"({abbrev[0]})" if abbrev else ""
    
    if compact:
        return f"{name} {abbrev_str} [CAS: {cas}]"
    
    lines = []
    lines.append(f"  Name: {name}")
    if abbrev:
        lines.append(f"  Abbreviation: {', '.join(abbrev)}")
    lines.append(f"  CAS: {cas}")
    
    if reagent_data.get('smiles'):
        lines.append(f"  SMILES: {reagent_data['smiles']}")
    
    if reagent_data.get('inchi_key'):
        lines.append(f"  InChIKey: {reagent_data['inchi_key']}")
    
    if reagent_data.get('aliases'):
        lines.append(f"  Aliases: {', '.join(reagent_data['aliases'][:3])}")
    
    if reagent_data.get('roles'):
        roles_str = ', '.join(reagent_data['roles'].keys())
        lines.append(f"  Roles: {roles_str}")
    
    return '\n'.join(lines)


def show_database_statistics():
    """Display comprehensive database statistics."""
    print_section("REAGENT DATABASE STATISTICS")
    
    # Get statistics from analytics module  
    # Use the correct directory (data/reagent_db, not data/reagents)
    from chemtools.reagent.lookup import get_data_dir
    reagent_db_dir = get_data_dir() / "reagent_db"
    stats = reagent.get_database_statistics(registry_dir=reagent_db_dir)
    
    print(f"Database Path: {stats.get('registry_dir', 'N/A')}")
    print(f"\nTotal Reagents: {stats.get('total_reagents', 0)}")
    
    # Show distribution by role
    print_subsection("Distribution by Role")
    by_type = stats.get('by_type', {})
    for role, count in sorted(by_type.items(), key=lambda x: x[1], reverse=True):
        pct = (count / stats['total_reagents'] * 100) if stats['total_reagents'] > 0 else 0
        print(f"  {role:20s}: {count:4d} ({pct:5.1f}%)")
    
    # Show coverage statistics
    print_subsection("Data Coverage")
    print(f"  With CAS numbers:     {stats.get('total_with_cas', 0):4d} ({stats.get('percent_with_cas', 0):.1f}%)")
    print(f"  With InChIKey:        {stats.get('total_with_inchikey', 0):4d} ({stats.get('percent_with_inchikey', 0):.1f}%)")
    print(f"  With SMILES:          {stats.get('total_with_smiles', 0):4d} ({stats.get('percent_with_smiles', 0):.1f}%)")
    print(f"  With abbreviations:   {stats.get('total_with_abbreviations', 0):4d} ({stats.get('percent_with_abbreviations', 0):.1f}%)")
    
    # Show top families
    print_subsection("Top 10 Reagent Families")
    top_families = stats.get('top_families', [])
    for family, count in top_families[:10]:
        print(f"  {family:30s}: {count:4d}")
    
    # Show multi-role reagents
    print_subsection("Multi-Role Reagents")
    print(f"  Reagents with multiple roles: {stats.get('multi_role_reagents', 0)}")
    roles_per_reagent = stats.get('roles_per_reagent', {})
    for num_roles, count in sorted(roles_per_reagent.items()):
        print(f"    {num_roles} role(s): {count} reagents")


def demonstrate_cas_search():
    """Demonstrate searching reagents by CAS number."""
    print_section("SEARCH BY CAS NUMBER - DEMONSTRATIONS")
    
    # Test cases: (name, CAS number, expected_role)
    test_cases = [
        ("Sodium tert-butoxide", "865-48-5", "base"),
        ("Toluene", "108-88-3", "solvent"),
        ("Tetrahydrofuran (THF)", "109-99-9", "solvent"),
        ("Palladium(II) acetate", "3375-31-3", "metal_precursor"),
        ("XPhos", "564483-18-7", "ligand"),
        ("tBuXPhos", "564483-18-7", "ligand"),  # Same CAS as XPhos - should find it
        ("Potassium carbonate", "584-08-7", "base"),
        ("Cesium carbonate", "534-17-8", "base"),
        ("Dimethylformamide (DMF)", "68-12-2", "solvent"),
    ]
    
    found_count = 0
    not_found = []
    
    for name, cas, expected_role in test_cases:
        print_subsection(f"Searching: {name} (CAS: {cas})")
        
        # Search in the expected role database
        result = reagent.find_reagent(cas, expected_role)
        
        if result:
            found_count += 1
            print(f"[OK] FOUND in {expected_role} database:")
            print(format_reagent(result))
        else:
            not_found.append((name, cas, expected_role))
            print(f"[X] NOT FOUND in {expected_role} database")
            
            # Try to find in any database
            all_types = reagent.get_all_reagent_types()
            found_elsewhere = False
            for rtype in all_types:
                if rtype != expected_role:
                    result = reagent.find_reagent(cas, rtype)
                    if result:
                        print(f"  -> But found in {rtype} database:")
                        print(format_reagent(result, compact=True))
                        found_elsewhere = True
                        break
            
            if not found_elsewhere:
                print(f"  -> Not found in any database")
    
    # Summary
    print_subsection("CAS Search Summary")
    print(f"  Searched: {len(test_cases)} reagents")
    print(f"  Found:    {found_count} ({found_count/len(test_cases)*100:.1f}%)")
    print(f"  Missing:  {len(not_found)}")
    
    if not_found:
        print("\n  Missing reagents:")
        for name, cas, role in not_found:
            print(f"    - {name} (CAS: {cas}) - expected in {role}")


def demonstrate_name_search():
    """Demonstrate searching reagents by name/abbreviation."""
    print_section("SEARCH BY NAME/ABBREVIATION - DEMONSTRATIONS")
    
    # Test cases: (search_term, expected_role)
    test_cases = [
        ("NaOtBu", "base"),
        ("sodium tert-butoxide", "base"),
        ("THF", "solvent"),
        ("tetrahydrofuran", "solvent"),
        ("Pd(OAc)2", "metal_precursor"),
        ("palladium acetate", "metal_precursor"),
        ("XPhos", "ligand"),
        ("2-dicyclohexylphosphino-2',4',6'-triisopropylbiphenyl", "ligand"),
        ("K2CO3", "base"),
        ("potassium carbonate", "base"),
        ("Cs2CO3", "base"),
        ("cesium carbonate", "base"),
        ("DMF", "solvent"),
        ("dimethylformamide", "solvent"),
        ("DMSO", "solvent"),
        ("toluene", "solvent"),
    ]
    
    found_count = 0
    not_found = []
    
    for search_term, expected_role in test_cases:
        print_subsection(f"Searching: '{search_term}' in {expected_role}")
        
        result = reagent.find_reagent(search_term, expected_role)
        
        if result:
            found_count += 1
            print(f"[OK] FOUND:")
            print(format_reagent(result, compact=True))
        else:
            not_found.append((search_term, expected_role))
            print(f"[X] NOT FOUND in {expected_role} database")
    
    # Summary
    print_subsection("Name Search Summary")
    print(f"  Searched: {len(test_cases)} terms")
    print(f"  Found:    {found_count} ({found_count/len(test_cases)*100:.1f}%)")
    print(f"  Missing:  {len(not_found)}")
    
    if not_found:
        print("\n  Missing reagents:")
        for term, role in not_found:
            print(f"    - '{term}' in {role}")


def demonstrate_role_based_search():
    """Demonstrate role-based filtering."""
    print_section("ROLE-BASED FILTERING - DEMONSTRATIONS")
    
    # Show reagents from each major role
    roles_to_show = ["base", "solvent", "ligand", "metal_precursor", "additive"]
    
    for role in roles_to_show:
        print_subsection(f"All {role.upper()}s in Database")
        
        all_reagents = reagent.get_all_reagents_by_type(role)
        count = len(all_reagents)
        
        print(f"  Total count: {count}")
        
        if count > 0:
            print(f"\n  First 5 examples:")
            for r in all_reagents[:5]:
                name = r.get('name', 'Unknown')
                cas = r.get('cas', 'No CAS')
                abbrev = r.get('abbreviation', [])
                abbrev_str = f" ({abbrev[0]})" if abbrev else ""
                try:
                    print(f"    - {name}{abbrev_str} [CAS: {cas}]")
                except UnicodeEncodeError:
                    # Handle special characters in reagent names on Windows
                    print(f"    - {name.encode('ascii', 'replace').decode('ascii')}{abbrev_str} [CAS: {cas}]")
        
        # Show family breakdown for this role
        family_stats = reagent.get_family_statistics(role)
        families = family_stats.get('families', [])
        
        if len(families) > 0:
            print(f"\n  Family breakdown ({family_stats.get('total_families', 0)} families):")
            for fam in families[:5]:
                print(f"    - {fam['name']:30s}: {fam['count']:3d} reagents")


def search_reagent_interactive(search_term: str, role: Optional[str] = None):
    """Search for a reagent interactively."""
    print_section(f"SEARCHING FOR: {search_term}")
    
    if role:
        # Search in specific role
        print(f"Searching in {role} database...")
        result = reagent.find_reagent(search_term, role)
        
        if result:
            print(f"\n[OK] FOUND:")
            print(format_reagent(result))
        else:
            print(f"\n[X] NOT FOUND in {role} database")
    else:
        # Search in all databases
        print("Searching in all databases...")
        all_types = reagent.get_all_reagent_types()
        found = False
        
        for rtype in all_types:
            result = reagent.find_reagent(search_term, rtype)
            if result:
                print(f"\n[OK] FOUND in {rtype} database:")
                print(format_reagent(result))
                found = True
                break
        
        if not found:
            print(f"\n[X] NOT FOUND in any database")


def export_reagent_inventory(output_file: str):
    """Export complete reagent inventory to JSON file."""
    print_section(f"EXPORTING REAGENT INVENTORY")
    
    print(f"Collecting all reagents...")
    
    inventory = {
        "total_reagents": 0,
        "by_role": {},
        "metadata": {
            "source": "chemtools.reagent database",
            "generated_by": "reagent_analytics.py"
        }
    }
    
    all_types = reagent.get_all_reagent_types()
    
    for rtype in all_types:
        reagents = reagent.get_all_reagents_by_type(rtype)
        inventory["by_role"][rtype] = reagents
        inventory["total_reagents"] += len(reagents)
        print(f"  - {rtype:20s}: {len(reagents):4d} reagents")
    
    # Write to file
    output_path = Path(output_file)
    with open(output_path, 'w', encoding='utf-8') as f:
        import json
        json.dump(inventory, f, indent=2, ensure_ascii=False)
    
    print(f"\n[OK] Exported {inventory['total_reagents']} reagents to: {output_path}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Reagent Database Analytics Tool")
    parser.add_argument("--search", help="Search for a reagent by name")
    parser.add_argument("--cas", help="Search for a reagent by CAS number")
    parser.add_argument("--role", help="Limit search to specific role (base, ligand, solvent, etc.)")
    parser.add_argument("--export", help="Export reagent inventory to JSON file")
    parser.add_argument("--full", action="store_true", help="Run full analytics and demonstrations")
    
    args = parser.parse_args()
    
    # Handle specific operations
    if args.search:
        search_reagent_interactive(args.search, args.role)
        return
    
    if args.cas:
        search_reagent_interactive(args.cas, args.role)
        return
    
    if args.export:
        export_reagent_inventory(args.export)
        return
    
    # Default: run full analytics and demonstrations
    print_section("REAGENT DATABASE ANALYTICS")
    print("\nThis tool demonstrates the chemtools.reagent module's capabilities:")
    print("  1. Database statistics and analysis")
    print("  2. Search by CAS number")
    print("  3. Search by name/abbreviation")
    print("  4. Role-based filtering")
    
    try:
        # Show database statistics
        show_database_statistics()
        
        # Demonstrate CAS search
        demonstrate_cas_search()
        
        # Demonstrate name search
        demonstrate_name_search()
        
        # Demonstrate role-based filtering
        demonstrate_role_based_search()
        
        print_section("ANALYTICS COMPLETE")
        print("\nFor more options:")
        print("  python app/reagent_analytics.py --search <name>")
        print("  python app/reagent_analytics.py --cas <cas_number>")
        print("  python app/reagent_analytics.py --role <role> --search <name>")
        print("  python app/reagent_analytics.py --export reagents.json")
        
    except Exception as e:
        print(f"\n[X] Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
