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

# Try to import colorama for colored output
try:
    from colorama import init, Fore, Back, Style
    init(autoreset=True)  # Auto-reset colors after each print
    COLORS_AVAILABLE = True
except ImportError:
    # Fallback: no colors
    COLORS_AVAILABLE = False
    class Fore:
        CYAN = YELLOW = GREEN = RED = BLUE = MAGENTA = WHITE = LIGHTBLUE_EX = LIGHTGREEN_EX = LIGHTYELLOW_EX = LIGHTMAGENTA_EX = LIGHTBLACK_EX = ""
    class Back:
        BLACK = ""
    class Style:
        BRIGHT = RESET_ALL = ""


def print_section(title: str, color=None):
    """Print a formatted section header with color."""
    width = 80
    if color is None:
        color = Fore.CYAN + Style.BRIGHT
    print()
    print(color + "=" * width)
    print(color + title.center(width))
    print(color + "=" * width + Style.RESET_ALL)
    print()


def print_subsection(title: str, color=None):
    """Print a formatted subsection header with color."""
    if color is None:
        color = Fore.YELLOW + Style.BRIGHT
    print()
    print(color + "-" * 80)
    print(color + f"  {title}")
    print(color + "-" * 80 + Style.RESET_ALL)


def format_reagent(reagent_data: Dict[str, Any], compact: bool = False) -> str:
    """Format reagent information for display with colors."""
    name = reagent_data.get('name', 'Unknown')
    cas = reagent_data.get('cas', 'No CAS')
    abbrev = reagent_data.get('abbreviation', [])
    abbrev_str = f" {Fore.LIGHTBLUE_EX}({abbrev[0]}){Style.RESET_ALL}" if abbrev else ""
    
    if compact:
        return f"{Fore.WHITE}{name}{Style.RESET_ALL}{abbrev_str} {Fore.LIGHTBLACK_EX}[CAS: {cas}]{Style.RESET_ALL}"
    
    lines = []
    lines.append(f"  {Fore.CYAN}Name:{Style.RESET_ALL} {Fore.WHITE + Style.BRIGHT}{name}{Style.RESET_ALL}")
    if abbrev:
        lines.append(f"  {Fore.CYAN}Abbreviation:{Style.RESET_ALL} {Fore.LIGHTBLUE_EX}{', '.join(abbrev)}{Style.RESET_ALL}")
    lines.append(f"  {Fore.CYAN}CAS:{Style.RESET_ALL} {Fore.YELLOW}{cas}{Style.RESET_ALL}")
    
    if reagent_data.get('smiles'):
        lines.append(f"  {Fore.CYAN}SMILES:{Style.RESET_ALL} {Fore.GREEN}{reagent_data['smiles']}{Style.RESET_ALL}")
    
    if reagent_data.get('inchi_key'):
        lines.append(f"  {Fore.CYAN}InChIKey:{Style.RESET_ALL} {Fore.GREEN}{reagent_data['inchi_key']}{Style.RESET_ALL}")
    
    if reagent_data.get('aliases'):
        aliases_str = ', '.join(reagent_data['aliases'][:3])
        lines.append(f"  {Fore.CYAN}Aliases:{Style.RESET_ALL} {Fore.LIGHTBLACK_EX}{aliases_str}{Style.RESET_ALL}")
    
    if reagent_data.get('roles'):
        roles_str = ', '.join(reagent_data['roles'].keys())
        lines.append(f"  {Fore.CYAN}Roles:{Style.RESET_ALL} {Fore.MAGENTA}{roles_str}{Style.RESET_ALL}")
    
    return '\n'.join(lines)


def show_database_statistics():
    """Display comprehensive database statistics with colors."""
    print_section("REAGENT DATABASE STATISTICS")
    
    # Get statistics from analytics module  
    from chemtools.reagent.lookup import get_data_dir
    reagent_db_dir = get_data_dir() / "reagent_db"
    stats = reagent.get_database_statistics(registry_dir=reagent_db_dir)
    
    print(f"{Fore.CYAN}Database Path:{Style.RESET_ALL} {Fore.WHITE}{stats.get('registry_dir', 'N/A')}{Style.RESET_ALL}")
    print(f"\n{Fore.CYAN}Total Reagents:{Style.RESET_ALL} {Fore.GREEN + Style.BRIGHT}{stats.get('total_reagents', 0)}{Style.RESET_ALL}")
    
    # Show distribution by role
    print_subsection("Distribution by Role")
    by_type = stats.get('by_type', {})
    for role, count in sorted(by_type.items(), key=lambda x: x[1], reverse=True):
        pct = (count / stats['total_reagents'] * 100) if stats['total_reagents'] > 0 else 0
        if count > 100:
            color = Fore.GREEN + Style.BRIGHT
        elif count > 50:
            color = Fore.GREEN
        elif count > 20:
            color = Fore.YELLOW
        else:
            color = Fore.WHITE
        print(f"  {Fore.MAGENTA}{role:20s}{Style.RESET_ALL}: {color}{count:4d}{Style.RESET_ALL} {Fore.LIGHTBLACK_EX}({pct:5.1f}%){Style.RESET_ALL}")
    
    # Show coverage statistics
    print_subsection("Data Coverage")
    def coverage_color(pct):
        if pct >= 90: return Fore.GREEN + Style.BRIGHT
        elif pct >= 70: return Fore.GREEN
        elif pct >= 50: return Fore.YELLOW
        else: return Fore.RED
    
    cas_pct = stats.get('percent_with_cas', 0)
    ink_pct = stats.get('percent_with_inchikey', 0)
    smi_pct = stats.get('percent_with_smiles', 0)
    abb_pct = stats.get('percent_with_abbreviations', 0)
    
    print(f"  {Fore.CYAN}With CAS numbers:{Style.RESET_ALL}     {stats.get('total_with_cas', 0):4d} {coverage_color(cas_pct)}({cas_pct:.1f}%){Style.RESET_ALL}")
    print(f"  {Fore.CYAN}With InChIKey:{Style.RESET_ALL}        {stats.get('total_with_inchikey', 0):4d} {coverage_color(ink_pct)}({ink_pct:.1f}%){Style.RESET_ALL}")
    print(f"  {Fore.CYAN}With SMILES:{Style.RESET_ALL}          {stats.get('total_with_smiles', 0):4d} {coverage_color(smi_pct)}({smi_pct:.1f}%){Style.RESET_ALL}")
    print(f"  {Fore.CYAN}With abbreviations:{Style.RESET_ALL}   {stats.get('total_with_abbreviations', 0):4d} {coverage_color(abb_pct)}({abb_pct:.1f}%){Style.RESET_ALL}")
    
    # Show top families
    print_subsection("Top 10 Reagent Families")
    top_families = stats.get('top_families', [])
    for i, (family, count) in enumerate(top_families[:10], 1):
        if i == 1:
            num_color = Fore.YELLOW + Style.BRIGHT
        elif i == 2:
            num_color = Fore.WHITE + Style.BRIGHT
        elif i == 3:
            num_color = Fore.LIGHTYELLOW_EX
        else:
            num_color = Fore.LIGHTBLACK_EX
        print(f"  {num_color}{i:2d}.{Style.RESET_ALL} {Fore.BLUE}{family:30s}{Style.RESET_ALL}: {Fore.GREEN}{count:4d}{Style.RESET_ALL}")
    
    # Show multi-role reagents
    print_subsection("Multi-Role Reagents")
    print(f"  {Fore.CYAN}Reagents with multiple roles:{Style.RESET_ALL} {Fore.YELLOW + Style.BRIGHT}{stats.get('multi_role_reagents', 0)}{Style.RESET_ALL}")
    roles_per_reagent = stats.get('roles_per_reagent', {})
    for num_roles, count in sorted(roles_per_reagent.items()):
        print(f"    {Fore.LIGHTBLACK_EX}{num_roles} role(s):{Style.RESET_ALL} {Fore.WHITE}{count}{Style.RESET_ALL} reagents")


def demonstrate_cas_search():
    """Demonstrate searching reagents by CAS number with colored output."""
    print_section("SEARCH BY CAS NUMBER - DEMONSTRATIONS", Fore.GREEN + Style.BRIGHT)
    
    test_cases = [
        ("Sodium tert-butoxide", "865-48-5", "base"),
        ("Toluene", "108-88-3", "solvent"),
        ("Tetrahydrofuran (THF)", "109-99-9", "solvent"),
        ("Palladium(II) acetate", "3375-31-3", "metal_catalyst"),
        ("XPhos", "564483-18-7", "ligand"),
        ("tBuXPhos", "564483-18-7", "ligand"),
        ("Potassium carbonate", "584-08-7", "base"),
        ("Cesium carbonate", "534-17-8", "base"),
        ("Dimethylformamide (DMF)", "68-12-2", "solvent"),
    ]
    
    found_count = 0
    not_found = []
    
    for name, cas, expected_role in test_cases:
        print_subsection(f"Searching: {name} {Fore.LIGHTBLACK_EX}(CAS: {cas}){Style.RESET_ALL}", Fore.LIGHTBLUE_EX)
        
        result = reagent.find_reagent(cas, expected_role)
        
        if result:
            found_count += 1
            print(f"{Fore.GREEN + Style.BRIGHT}[OK] FOUND{Style.RESET_ALL} in {Fore.MAGENTA}{expected_role}{Style.RESET_ALL} database:")
            print(format_reagent(result))
        else:
            not_found.append((name, cas, expected_role))
            print(f"{Fore.RED + Style.BRIGHT}[X] NOT FOUND{Style.RESET_ALL} in {Fore.MAGENTA}{expected_role}{Style.RESET_ALL} database")
            
            all_types = reagent.get_all_reagent_types()
            found_elsewhere = False
            for rtype in all_types:
                if rtype != expected_role:
                    result = reagent.find_reagent(cas, rtype)
                    if result:
                        print(f"  {Fore.YELLOW}-> But found in {Fore.MAGENTA}{rtype}{Fore.YELLOW} database:{Style.RESET_ALL}")
                        print(f"    {format_reagent(result, compact=True)}")
                        found_elsewhere = True
                        break
            
            if not found_elsewhere:
                print(f"  {Fore.RED}-> Not found in any database{Style.RESET_ALL}")
    
    # Summary
    print_subsection("CAS Search Summary", Fore.CYAN + Style.BRIGHT)
    total = len(test_cases)
    success_rate = (found_count / total * 100) if total > 0 else 0
    success_color = Fore.GREEN + Style.BRIGHT if success_rate == 100 else (Fore.YELLOW if success_rate >= 75 else Fore.RED)
    
    print(f"  {Fore.CYAN}Searched:{Style.RESET_ALL} {Fore.WHITE}{total}{Style.RESET_ALL} reagents")
    print(f"  {Fore.CYAN}Found:{Style.RESET_ALL}    {Fore.GREEN + Style.BRIGHT}{found_count}{Style.RESET_ALL} {success_color}({success_rate:.1f}%){Style.RESET_ALL}")
    print(f"  {Fore.CYAN}Missing:{Style.RESET_ALL}  {Fore.RED if not_found else Fore.GREEN}{len(not_found)}{Style.RESET_ALL}")
    
    if not_found:
        print(f"\n  {Fore.RED + Style.BRIGHT}Missing reagents:{Style.RESET_ALL}")
        for name, cas, role in not_found:
            print(f"    {Fore.RED}*{Style.RESET_ALL} {name} {Fore.LIGHTBLACK_EX}(CAS: {cas}){Style.RESET_ALL} - expected in {Fore.MAGENTA}{role}{Style.RESET_ALL}")


def demonstrate_name_search():
    """Demonstrate searching reagents by name/abbreviation with colored output."""
    print_section("SEARCH BY NAME/ABBREVIATION - DEMONSTRATIONS", Fore.BLUE + Style.BRIGHT)
    
    test_cases = [
        ("NaOtBu", "base"),
        ("sodium tert-butoxide", "base"),
        ("THF", "solvent"),
        ("tetrahydrofuran", "solvent"),
        ("Pd(OAc)2", "metal_catalyst"),
        ("palladium acetate", "metal_catalyst"),
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
        print_subsection(f"Searching: {Fore.WHITE + Style.BRIGHT}'{search_term}'{Style.RESET_ALL} in {Fore.MAGENTA}{expected_role}{Style.RESET_ALL}", Fore.LIGHTBLUE_EX)
        
        result = reagent.find_reagent(search_term, expected_role)
        
        if result:
            found_count += 1
            print(f"{Fore.GREEN + Style.BRIGHT}[OK] FOUND:{Style.RESET_ALL}")
            print(f"  {format_reagent(result, compact=True)}")
        else:
            not_found.append((search_term, expected_role))
            print(f"{Fore.RED + Style.BRIGHT}[X] NOT FOUND{Style.RESET_ALL} in {Fore.MAGENTA}{expected_role}{Style.RESET_ALL} database")
    
    # Summary
    print_subsection("Name Search Summary", Fore.CYAN + Style.BRIGHT)
    total = len(test_cases)
    success_rate = (found_count / total * 100) if total > 0 else 0
    success_color = Fore.GREEN + Style.BRIGHT if success_rate == 100 else (Fore.YELLOW if success_rate >= 75 else Fore.RED)
    
    print(f"  {Fore.CYAN}Searched:{Style.RESET_ALL} {Fore.WHITE}{total}{Style.RESET_ALL} terms")
    print(f"  {Fore.CYAN}Found:{Style.RESET_ALL}    {Fore.GREEN + Style.BRIGHT}{found_count}{Style.RESET_ALL} {success_color}({success_rate:.1f}%){Style.RESET_ALL}")
    print(f"  {Fore.CYAN}Missing:{Style.RESET_ALL}  {Fore.RED if not_found else Fore.GREEN}{len(not_found)}{Style.RESET_ALL}")
    
    if not_found:
        print(f"\n  {Fore.RED + Style.BRIGHT}Missing reagents:{Style.RESET_ALL}")
        for term, role in not_found:
            print(f"    {Fore.RED}*{Style.RESET_ALL} {Fore.WHITE}'{term}'{Style.RESET_ALL} in {Fore.MAGENTA}{role}{Style.RESET_ALL}")


def demonstrate_role_based_search():
    """Demonstrate role-based filtering with colored output."""
    print_section("ROLE-BASED FILTERING - DEMONSTRATIONS", Fore.MAGENTA + Style.BRIGHT)
    
    roles_to_show = ["base", "solvent", "ligand", "metal_catalyst", "additive"]
    
    for role in roles_to_show:
        print_subsection(f"All {Fore.MAGENTA + Style.BRIGHT}{role.upper()}s{Style.RESET_ALL} in Database", Fore.LIGHTMAGENTA_EX)
        
        all_reagents = reagent.get_all_reagents_by_type(role)
        count = len(all_reagents)
        
        print(f"  {Fore.CYAN}Total count:{Style.RESET_ALL} {Fore.GREEN + Style.BRIGHT}{count}{Style.RESET_ALL}")
        
        if count > 0:
            print(f"\n  {Fore.YELLOW}First 5 examples:{Style.RESET_ALL}")
            for i, r in enumerate(all_reagents[:5], 1):
                name = r.get('name', 'Unknown')
                cas = r.get('cas', 'No CAS')
                abbrev = r.get('abbreviation', [])
                abbrev_str = f" {Fore.LIGHTBLUE_EX}({abbrev[0]}){Style.RESET_ALL}" if abbrev else ""
                try:
                    print(f"    {Fore.LIGHTBLACK_EX}{i}.{Style.RESET_ALL} {Fore.WHITE}{name}{Style.RESET_ALL}{abbrev_str} {Fore.LIGHTBLACK_EX}[CAS: {cas}]{Style.RESET_ALL}")
                except UnicodeEncodeError:
                    print(f"    {Fore.LIGHTBLACK_EX}{i}.{Style.RESET_ALL} {Fore.WHITE}{name.encode('ascii', 'replace').decode('ascii')}{Style.RESET_ALL}{abbrev_str} {Fore.LIGHTBLACK_EX}[CAS: {cas}]{Style.RESET_ALL}")
        
        family_stats = reagent.get_family_statistics(role)
        families = family_stats.get('families', [])
        
        if len(families) > 0:
            print(f"\n  {Fore.YELLOW}Family breakdown{Style.RESET_ALL} {Fore.LIGHTBLACK_EX}({family_stats.get('total_families', 0)} families):{Style.RESET_ALL}")
            for i, fam in enumerate(families[:5], 1):
                print(f"    {Fore.LIGHTBLACK_EX}{i}.{Style.RESET_ALL} {Fore.BLUE}{fam['name']:30s}{Style.RESET_ALL}: {Fore.GREEN}{fam['count']:3d}{Style.RESET_ALL} reagents")


def search_reagent_interactive(search_term: str, role: Optional[str] = None):
    """Search for a reagent interactively with colored output."""
    print_section(f"SEARCHING FOR: {Fore.WHITE + Style.BRIGHT}{search_term}{Style.RESET_ALL}", Fore.CYAN + Style.BRIGHT)
    
    if role:
        print(f"Searching in {Fore.MAGENTA}{role}{Style.RESET_ALL} database...")
        result = reagent.find_reagent(search_term, role)
        
        if result:
            print(f"\n{Fore.GREEN + Style.BRIGHT}[OK] FOUND:{Style.RESET_ALL}")
            print(format_reagent(result))
        else:
            print(f"\n{Fore.RED + Style.BRIGHT}[X] NOT FOUND{Style.RESET_ALL} in {Fore.MAGENTA}{role}{Style.RESET_ALL} database")
    else:
        print(f"{Fore.CYAN}Searching in all databases...{Style.RESET_ALL}")
        all_types = reagent.get_all_reagent_types()
        found = False
        
        for rtype in all_types:
            result = reagent.find_reagent(search_term, rtype)
            if result:
                print(f"\n{Fore.GREEN + Style.BRIGHT}[OK] FOUND{Style.RESET_ALL} in {Fore.MAGENTA}{rtype}{Style.RESET_ALL} database:")
                print(format_reagent(result))
                found = True
                break
        
        if not found:
            print(f"\n{Fore.RED + Style.BRIGHT}[X] NOT FOUND{Style.RESET_ALL} in any database")


def export_reagent_inventory(output_file: str):
    """Export complete reagent inventory to JSON file with colored output."""
    print_section(f"EXPORTING REAGENT INVENTORY", Fore.YELLOW + Style.BRIGHT)
    
    print(f"{Fore.CYAN}Collecting all reagents...{Style.RESET_ALL}\n")
    
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
        print(f"  {Fore.MAGENTA}{rtype:20s}{Style.RESET_ALL}: {Fore.GREEN}{len(reagents):4d}{Style.RESET_ALL} reagents")
    
    output_path = Path(output_file)
    with open(output_path, 'w', encoding='utf-8') as f:
        import json
        json.dump(inventory, f, indent=2, ensure_ascii=False)
    
    print(f"\n{Fore.GREEN + Style.BRIGHT}[OK] Exported {inventory['total_reagents']} reagents{Style.RESET_ALL} to: {Fore.CYAN}{output_path}{Style.RESET_ALL}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Reagent Database Analytics Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--search", help="Search for a reagent by name")
    parser.add_argument("--cas", help="Search for a reagent by CAS number")
    parser.add_argument("--role", help="Limit search to specific role (base, ligand, solvent, etc.)")
    parser.add_argument("--export", help="Export reagent inventory to JSON file")
    parser.add_argument("--full", action="store_true", help="Run full analytics and demonstrations")
    
    args = parser.parse_args()
    
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
    print_section("REAGENT DATABASE ANALYTICS", Fore.CYAN + Style.BRIGHT)
    print(f"\n{Fore.WHITE}This tool demonstrates the {Fore.CYAN + Style.BRIGHT}chemtools.reagent{Style.RESET_ALL + Fore.WHITE} module's capabilities:{Style.RESET_ALL}")
    print(f"  {Fore.GREEN}1.{Style.RESET_ALL} Database statistics and analysis")
    print(f"  {Fore.GREEN}2.{Style.RESET_ALL} Search by CAS number")
    print(f"  {Fore.GREEN}3.{Style.RESET_ALL} Search by name/abbreviation")
    print(f"  {Fore.GREEN}4.{Style.RESET_ALL} Role-based filtering")
    
    try:
        show_database_statistics()
        demonstrate_cas_search()
        demonstrate_name_search()
        demonstrate_role_based_search()
        
        print_section("ANALYTICS COMPLETE", Fore.GREEN + Style.BRIGHT)
        print(f"\n{Fore.YELLOW + Style.BRIGHT}For more options:{Style.RESET_ALL}")
        print(f"  {Fore.CYAN}python app/reagent_analytics.py --search <name>{Style.RESET_ALL}")
        print(f"  {Fore.CYAN}python app/reagent_analytics.py --cas <cas_number>{Style.RESET_ALL}")
        print(f"  {Fore.CYAN}python app/reagent_analytics.py --role <role> --search <name>{Style.RESET_ALL}")
        print(f"  {Fore.CYAN}python app/reagent_analytics.py --export reagents.json{Style.RESET_ALL}")
        
    except Exception as e:
        print(f"\n{Fore.RED + Style.BRIGHT}[X] Error during analysis:{Style.RESET_ALL} {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
