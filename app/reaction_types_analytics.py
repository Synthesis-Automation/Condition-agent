#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reaction Types Database Analytics and Summary Tool
===================================================

This tool provides comprehensive analytics and exploration of the reaction types taxonomy database:
1. Providing comprehensive database statistics
2. Listing reaction types by category
3. Showing SMARTS pattern coverage
4. Searching for specific reaction types

Usage:
    python app/reaction_types_analytics.py                    # Full analytics
    python app/reaction_types_analytics.py --search "suzuki"  # Search by name/alias
    python app/reaction_types_analytics.py --category c_c     # Filter by category
    python app/reaction_types_analytics.py --smarts           # Show only types with SMARTS
    python app/reaction_types_analytics.py --id Suzuki_miyaura # Show details for specific ID
"""

import sys
import json
from pathlib import Path
from typing import Dict, List, Any, Optional
import argparse
from collections import Counter

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

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


try:
    from chemtools.taxonomy import reaction_catalog as _reaction_catalog
    DEFAULT_REACTION_TYPES_PATH = _reaction_catalog.REACTION_TYPES_FILE
except Exception:
    _reaction_catalog = None
    DEFAULT_REACTION_TYPES_PATH = (
        Path(__file__).parent.parent / "chemtools" / "taxonomy" / "data" / "reaction_types.v4.0.json"
    )

LEGACY_REACTION_TYPES_PATH = Path(__file__).parent.parent / "chemtools" / "taxonomy" / "data" / "reaction_types.json"
REACTION_TYPES_PATH = DEFAULT_REACTION_TYPES_PATH


def resolve_reaction_types_path() -> Path:
    """Resolve the reaction types taxonomy file."""
    candidates = [DEFAULT_REACTION_TYPES_PATH, LEGACY_REACTION_TYPES_PATH]
    for path in candidates:
        if path and path.exists():
            return path
    tried = ", ".join(str(path) for path in candidates if path)
    raise FileNotFoundError(f"Reaction types file not found. Tried: {tried}")


def load_reaction_types() -> List[Dict[str, Any]]:
    """Load reaction types from JSON file."""
    path = resolve_reaction_types_path()
    global REACTION_TYPES_PATH
    REACTION_TYPES_PATH = path

    with open(path, 'r', encoding='utf-8') as f:
        payload = json.load(f)

    if isinstance(payload, list):
        return payload

    if isinstance(payload, dict):
        reaction_types = payload.get("reaction_types")
        if isinstance(reaction_types, list):
            return reaction_types
        reaction_types = payload.get("reactions")
        if isinstance(reaction_types, list):
            return reaction_types

        flattened: List[Dict[str, Any]] = []
        for category, category_payload in payload.items():
            if not isinstance(category_payload, dict):
                continue
            reactions = category_payload.get("reactions")
            if not isinstance(reactions, list):
                continue
            for entry in reactions:
                if not isinstance(entry, dict):
                    continue
                normalized = dict(entry)
                normalized.setdefault("category", category)
                flattened.append(normalized)
        if flattened:
            return flattened

    raise ValueError(f"Unrecognized reaction types JSON layout in {path}")


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


def format_reaction_type(rtype: Dict[str, Any], compact: bool = False) -> str:
    """Format reaction type information for display with colors."""
    rid = rtype.get('id', 'no-id')
    if rtype.get('name'):
        name = rtype.get('name')
    elif _reaction_catalog is not None:
        name = _reaction_catalog.format_reaction_id_display(rid)
    else:
        name = rid
    category = rtype.get('category', 'uncategorized')
    aliases = rtype.get('aliases', [])
    aliases_str = f" {Fore.LIGHTBLUE_EX}({', '.join(aliases[:2])}){Style.RESET_ALL}" if aliases else ""
    
    if compact:
        has_smarts = "Y" if rtype.get('smarts') else "N"
        smarts_color = Fore.GREEN if rtype.get('smarts') else Fore.RED
        return f"{Fore.WHITE}{name}{Style.RESET_ALL}{aliases_str} {smarts_color}[SMARTS: {has_smarts}]{Style.RESET_ALL}"
    
    lines = []
    lines.append(f"  {Fore.CYAN}ID:{Style.RESET_ALL} {Fore.YELLOW + Style.BRIGHT}{rid}{Style.RESET_ALL}")
    lines.append(f"  {Fore.CYAN}Name:{Style.RESET_ALL} {Fore.WHITE + Style.BRIGHT}{name}{Style.RESET_ALL}")
    lines.append(f"  {Fore.CYAN}Category:{Style.RESET_ALL} {Fore.MAGENTA}{category}{Style.RESET_ALL}")
    
    if aliases:
        lines.append(f"  {Fore.CYAN}Aliases:{Style.RESET_ALL} {Fore.LIGHTBLUE_EX}{', '.join(aliases)}{Style.RESET_ALL}")
    
    if rtype.get('description'):
        lines.append(f"  {Fore.CYAN}Description:{Style.RESET_ALL} {Fore.WHITE}{rtype['description']}{Style.RESET_ALL}")
    
    if rtype.get('smarts'):
        lines.append(f"  {Fore.CYAN}SMARTS:{Style.RESET_ALL} {Fore.GREEN}{rtype['smarts']}{Style.RESET_ALL}")
    else:
        lines.append(f"  {Fore.CYAN}SMARTS:{Style.RESET_ALL} {Fore.RED}Not defined{Style.RESET_ALL}")
    
    if rtype.get('reference_reactions'):
        lines.append(f"  {Fore.CYAN}Reference Reactions:{Style.RESET_ALL}")
        for rxn in rtype['reference_reactions'][:2]:
            lines.append(f"    {Fore.LIGHTBLACK_EX}{rxn}{Style.RESET_ALL}")
    
    reactants = rtype.get('reactants')
    if isinstance(reactants, list) and reactants:
        lines.append(f"  {Fore.CYAN}Reactants:{Style.RESET_ALL} {Fore.WHITE}{', '.join(reactants)}{Style.RESET_ALL}")
    elif isinstance(reactants, dict) and reactants:
        slots = ", ".join(sorted(reactants.keys()))
        if slots:
            lines.append(f"  {Fore.CYAN}Reactants:{Style.RESET_ALL} {Fore.WHITE}{slots}{Style.RESET_ALL}")
    
    if rtype.get('catalysts'):
        lines.append(f"  {Fore.CYAN}Catalysts:{Style.RESET_ALL} {Fore.YELLOW}{', '.join(rtype['catalysts'])}{Style.RESET_ALL}")
    
    if rtype.get('conditions'):
        lines.append(f"  {Fore.CYAN}Conditions:{Style.RESET_ALL} {Fore.LIGHTBLACK_EX}{rtype['conditions']}{Style.RESET_ALL}")
    
    return '\n'.join(lines)


def show_database_statistics(reaction_types: List[Dict[str, Any]]):
    """Display comprehensive database statistics with colors."""
    print_section("REACTION TYPES DATABASE STATISTICS")
    
    total = len(reaction_types)
    print(f"{Fore.CYAN}Database Path:{Style.RESET_ALL} {Fore.WHITE}{REACTION_TYPES_PATH}{Style.RESET_ALL}")
    print(f"\n{Fore.CYAN}Total Reaction Types:{Style.RESET_ALL} {Fore.GREEN + Style.BRIGHT}{total}{Style.RESET_ALL}")
    
    # Category distribution
    print_subsection("Distribution by Category")
    categories = Counter(r.get('category', 'uncategorized') for r in reaction_types)
    for category, count in sorted(categories.items(), key=lambda x: x[1], reverse=True):
        pct = (count / total * 100) if total > 0 else 0
        if count > 10:
            color = Fore.GREEN + Style.BRIGHT
        elif count > 5:
            color = Fore.GREEN
        elif count > 2:
            color = Fore.YELLOW
        else:
            color = Fore.WHITE
        # Format category name nicely
        category_display = category.replace('_', ' ').title()
        print(f"  {Fore.MAGENTA}{category_display:40s}{Style.RESET_ALL}: {color}{count:3d}{Style.RESET_ALL} {Fore.LIGHTBLACK_EX}({pct:5.1f}%){Style.RESET_ALL}")
    
    # SMARTS coverage
    print_subsection("SMARTS Pattern Coverage")
    with_smarts = sum(1 for r in reaction_types if r.get('smarts'))
    without_smarts = total - with_smarts
    smarts_pct = (with_smarts / total * 100) if total > 0 else 0
    
    def coverage_color(pct):
        if pct >= 80: return Fore.GREEN + Style.BRIGHT
        elif pct >= 60: return Fore.GREEN
        elif pct >= 40: return Fore.YELLOW
        else: return Fore.RED
    
    print(f"  {Fore.CYAN}With SMARTS patterns:{Style.RESET_ALL}    {Fore.GREEN + Style.BRIGHT}{with_smarts:3d}{Style.RESET_ALL} {coverage_color(smarts_pct)}({smarts_pct:.1f}%){Style.RESET_ALL}")
    print(f"  {Fore.CYAN}Without SMARTS patterns:{Style.RESET_ALL} {Fore.RED}{without_smarts:3d}{Style.RESET_ALL} {Fore.LIGHTBLACK_EX}({100-smarts_pct:.1f}%){Style.RESET_ALL}")
    
    # Reference reactions coverage
    print_subsection("Reference Reactions Coverage")
    with_ref = sum(1 for r in reaction_types if r.get('reference_reactions'))
    ref_pct = (with_ref / total * 100) if total > 0 else 0
    print(f"  {Fore.CYAN}With reference reactions:{Style.RESET_ALL} {Fore.GREEN}{with_ref:3d}{Style.RESET_ALL} {coverage_color(ref_pct)}({ref_pct:.1f}%){Style.RESET_ALL}")
    
    # Catalyst distribution
    print_subsection("Top Catalysts")
    all_catalysts = []
    for r in reaction_types:
        if r.get('catalysts'):
            all_catalysts.extend(r['catalysts'])
    catalyst_counts = Counter(all_catalysts)
    for i, (catalyst, count) in enumerate(catalyst_counts.most_common(10), 1):
        if i == 1:
            num_color = Fore.YELLOW + Style.BRIGHT
        elif i == 2:
            num_color = Fore.WHITE + Style.BRIGHT
        elif i == 3:
            num_color = Fore.LIGHTYELLOW_EX
        else:
            num_color = Fore.LIGHTBLACK_EX
        print(f"  {num_color}{i:2d}.{Style.RESET_ALL} {Fore.BLUE}{catalyst:30s}{Style.RESET_ALL}: {Fore.GREEN}{count:3d}{Style.RESET_ALL} reactions")
    
    # Alias statistics
    print_subsection("Alias Statistics")
    with_aliases = sum(1 for r in reaction_types if r.get('aliases'))
    total_aliases = sum(len(r.get('aliases', [])) for r in reaction_types)
    avg_aliases = total_aliases / total if total > 0 else 0
    print(f"  {Fore.CYAN}Reaction types with aliases:{Style.RESET_ALL} {Fore.GREEN}{with_aliases}{Style.RESET_ALL}")
    print(f"  {Fore.CYAN}Total aliases:{Style.RESET_ALL}              {Fore.GREEN}{total_aliases}{Style.RESET_ALL}")
    print(f"  {Fore.CYAN}Average aliases per type:{Style.RESET_ALL}   {Fore.YELLOW}{avg_aliases:.2f}{Style.RESET_ALL}")


def show_reactions_by_category(reaction_types: List[Dict[str, Any]], filter_category: Optional[str] = None):
    """Show reaction types organized by category."""
    print_section("REACTION TYPES BY CATEGORY", Fore.MAGENTA + Style.BRIGHT)
    
    # Group by category
    by_category: Dict[str, List[Dict]] = {}
    for r in reaction_types:
        cat = r.get('category', 'uncategorized')
        if cat not in by_category:
            by_category[cat] = []
        by_category[cat].append(r)
    
    # Filter if specified
    if filter_category:
        matching = {k: v for k, v in by_category.items() if filter_category.lower() in k.lower()}
        if not matching:
            print(f"{Fore.RED}No categories matching '{filter_category}'{Style.RESET_ALL}")
            print(f"\n{Fore.CYAN}Available categories:{Style.RESET_ALL}")
            for cat in sorted(by_category.keys()):
                print(f"  {Fore.LIGHTBLACK_EX}-{Style.RESET_ALL} {cat}")
            return
        by_category = matching
    
    for category in sorted(by_category.keys()):
        reactions = by_category[category]
        category_display = category.replace('_', ' ').title()
        print_subsection(f"{category_display} ({len(reactions)} types)", Fore.LIGHTMAGENTA_EX)
        
        for r in sorted(reactions, key=lambda x: x.get('name') or x.get('id', '')):
            rid = r.get('id', 'no-id')
            if r.get('name'):
                name = r.get('name')
            elif _reaction_catalog is not None:
                name = _reaction_catalog.format_reaction_id_display(rid)
            else:
                name = rid
            has_smarts = "Y" if r.get('smarts') else "N"
            smarts_color = Fore.GREEN if r.get('smarts') else Fore.RED
            aliases = r.get('aliases', [])
            alias_str = f" {Fore.LIGHTBLACK_EX}({aliases[0]}){Style.RESET_ALL}" if aliases else ""
            
            print(f"  {smarts_color}[{has_smarts}]{Style.RESET_ALL} {Fore.YELLOW}{rid:30s}{Style.RESET_ALL} {Fore.WHITE}{name}{Style.RESET_ALL}{alias_str}")


def show_smarts_coverage(reaction_types: List[Dict[str, Any]]):
    """Show detailed SMARTS pattern coverage."""
    print_section("SMARTS PATTERN COVERAGE", Fore.GREEN + Style.BRIGHT)
    
    with_smarts = [r for r in reaction_types if r.get('smarts')]
    without_smarts = [r for r in reaction_types if not r.get('smarts')]
    
    print_subsection(f"Reaction Types WITH SMARTS ({len(with_smarts)})", Fore.GREEN)
    for r in sorted(with_smarts, key=lambda x: x.get('category', '')):
        rid = r.get('id', 'no-id')
        if r.get('name'):
            name = r.get('name')
        elif _reaction_catalog is not None:
            name = _reaction_catalog.format_reaction_id_display(rid)
        else:
            name = rid
        category = r.get('category', 'uncategorized')
        smarts = r.get('smarts', '')
        # Truncate long SMARTS
        smarts_display = smarts[:50] + "..." if len(smarts) > 50 else smarts
        print(f"  {Fore.YELLOW}{rid:30s}{Style.RESET_ALL} {Fore.LIGHTBLACK_EX}{smarts_display}{Style.RESET_ALL}")
    
    print_subsection(f"Reaction Types WITHOUT SMARTS ({len(without_smarts)})", Fore.RED)
    # Group by category
    by_cat: Dict[str, List[str]] = {}
    for r in without_smarts:
        cat = r.get('category', 'uncategorized')
        if cat not in by_cat:
            by_cat[cat] = []
        if r.get('name'):
            by_cat[cat].append(r.get('name'))
        elif _reaction_catalog is not None:
            by_cat[cat].append(_reaction_catalog.format_reaction_id_display(r.get('id', 'Unknown')))
        else:
            by_cat[cat].append(r.get('id', 'Unknown'))
    
    for cat in sorted(by_cat.keys()):
        names = by_cat[cat]
        cat_display = cat.replace('_', ' ').title()
        print(f"\n  {Fore.MAGENTA}{cat_display}:{Style.RESET_ALL}")
        for name in sorted(names):
            print(f"    {Fore.RED}-{Style.RESET_ALL} {Fore.WHITE}{name}{Style.RESET_ALL}")


def search_reaction_types(reaction_types: List[Dict[str, Any]], query: str):
    """Search for reaction types by name, alias, or ID."""
    print_section(f"SEARCHING FOR: {Fore.WHITE + Style.BRIGHT}{query}{Style.RESET_ALL}", Fore.CYAN + Style.BRIGHT)
    
    query_lower = query.lower()
    matches = []
    
    for r in reaction_types:
        # Search in ID
        if query_lower in r.get('id', '').lower():
            matches.append((r, 'id'))
            continue
        # Search in name
        name_value = r.get('name')
        if not name_value and _reaction_catalog is not None:
            name_value = _reaction_catalog.format_reaction_id_display(r.get('id', ''))
        if query_lower in (name_value or r.get('id', '')).lower():
            matches.append((r, 'name'))
            continue
        # Search in aliases
        for alias in r.get('aliases', []):
            if query_lower in alias.lower():
                matches.append((r, 'alias'))
                break
        # Search in description
        if query_lower in r.get('description', '').lower():
            matches.append((r, 'description'))
    
    if matches:
        print(f"{Fore.GREEN + Style.BRIGHT}Found {len(matches)} matching reaction type(s):{Style.RESET_ALL}\n")
        for r, match_type in matches:
            print(f"{Fore.LIGHTBLACK_EX}(matched in {match_type}){Style.RESET_ALL}")
            print(format_reaction_type(r))
            print()
    else:
        print(f"{Fore.RED + Style.BRIGHT}No reaction types found matching '{query}'{Style.RESET_ALL}")
        
        # Suggest similar
        print(f"\n{Fore.CYAN}Available reaction type IDs:{Style.RESET_ALL}")
        for r in reaction_types[:10]:
            print(f"  {Fore.LIGHTBLACK_EX}-{Style.RESET_ALL} {r.get('id', 'no-id')}")
        if len(reaction_types) > 10:
            print(f"  {Fore.LIGHTBLACK_EX}... and {len(reaction_types) - 10} more{Style.RESET_ALL}")


def show_reaction_details(reaction_types: List[Dict[str, Any]], reaction_id: str):
    """Show detailed information for a specific reaction type."""
    print_section(f"REACTION TYPE DETAILS: {Fore.WHITE + Style.BRIGHT}{reaction_id}{Style.RESET_ALL}", Fore.CYAN + Style.BRIGHT)
    
    for r in reaction_types:
        if r.get('id', '').lower() == reaction_id.lower():
            print(format_reaction_type(r))
            return
    
    print(f"{Fore.RED + Style.BRIGHT}Reaction type '{reaction_id}' not found{Style.RESET_ALL}")
    
    # Suggest similar IDs
    similar = [r.get('id') for r in reaction_types if reaction_id.lower() in r.get('id', '').lower()]
    if similar:
        print(f"\n{Fore.CYAN}Did you mean:{Style.RESET_ALL}")
        for sid in similar[:5]:
            print(f"  {Fore.LIGHTBLACK_EX}-{Style.RESET_ALL} {sid}")


def export_summary(reaction_types: List[Dict[str, Any]], output_file: str):
    """Export reaction types summary to JSON file."""
    print_section("EXPORTING REACTION TYPES SUMMARY", Fore.YELLOW + Style.BRIGHT)
    
    # Build summary
    summary = {
        "total_types": len(reaction_types),
        "with_smarts": sum(1 for r in reaction_types if r.get('smarts')),
        "without_smarts": sum(1 for r in reaction_types if not r.get('smarts')),
        "by_category": {},
        "reaction_types": []
    }
    
    # Group by category
    for r in reaction_types:
        cat = r.get('category', 'uncategorized')
        if cat not in summary["by_category"]:
            summary["by_category"][cat] = []
        summary["by_category"][cat].append({
            "id": r.get('id'),
            "name": (
                r.get('name')
                or (_reaction_catalog.format_reaction_id_display(r.get('id', '')) if _reaction_catalog is not None else r.get('id'))
            ),
            "has_smarts": bool(r.get('smarts')),
            "aliases": r.get('aliases', [])
        })
        summary["reaction_types"].append(r)
    
    output_path = Path(output_file)
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    
    print(f"{Fore.GREEN + Style.BRIGHT}[OK] Exported summary{Style.RESET_ALL} to: {Fore.CYAN}{output_path}{Style.RESET_ALL}")
    print(f"\n{Fore.CYAN}Summary:{Style.RESET_ALL}")
    print(f"  {Fore.WHITE}Total types:{Style.RESET_ALL} {summary['total_types']}")
    print(f"  {Fore.WHITE}With SMARTS:{Style.RESET_ALL} {summary['with_smarts']}")
    print(f"  {Fore.WHITE}Categories:{Style.RESET_ALL} {len(summary['by_category'])}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Reaction Types Database Analytics Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python app/reaction_types_analytics.py                     # Full analytics
  python app/reaction_types_analytics.py --search suzuki     # Search by name
  python app/reaction_types_analytics.py --category c_c      # Filter by category
  python app/reaction_types_analytics.py --smarts            # Show SMARTS coverage
  python app/reaction_types_analytics.py --id C_N_Coupling  # Show details
  python app/reaction_types_analytics.py --export summary.json
        """
    )
    parser.add_argument("--search", help="Search for reaction types by name/alias/ID")
    parser.add_argument("--category", help="Filter by category (partial match)")
    parser.add_argument("--smarts", action="store_true", help="Show SMARTS pattern coverage")
    parser.add_argument("--id", dest="reaction_id", help="Show details for specific reaction ID")
    parser.add_argument("--export", help="Export summary to JSON file")
    parser.add_argument("--full", action="store_true", help="Run full analytics")
    
    args = parser.parse_args()
    
    try:
        reaction_types = load_reaction_types()
    except FileNotFoundError as e:
        print(f"{Fore.RED + Style.BRIGHT}[X] Error:{Style.RESET_ALL} {e}")
        return 1
    
    if args.search:
        search_reaction_types(reaction_types, args.search)
        return 0
    
    if args.reaction_id:
        show_reaction_details(reaction_types, args.reaction_id)
        return 0
    
    if args.smarts:
        show_smarts_coverage(reaction_types)
        return 0
    
    if args.category:
        show_reactions_by_category(reaction_types, args.category)
        return 0
    
    if args.export:
        export_summary(reaction_types, args.export)
        return 0
    
    # Default: run full analytics
    print_section("REACTION TYPES DATABASE ANALYTICS", Fore.CYAN + Style.BRIGHT)
    print(f"\n{Fore.WHITE}This tool provides analytics for the {Fore.CYAN + Style.BRIGHT}{REACTION_TYPES_PATH.name}{Style.RESET_ALL + Fore.WHITE} database:{Style.RESET_ALL}")
    print(f"  {Fore.GREEN}1.{Style.RESET_ALL} Database statistics and coverage")
    print(f"  {Fore.GREEN}2.{Style.RESET_ALL} Reaction types by category")
    print(f"  {Fore.GREEN}3.{Style.RESET_ALL} SMARTS pattern analysis")
    print(f"  {Fore.GREEN}4.{Style.RESET_ALL} Search and filtering")
    
    try:
        show_database_statistics(reaction_types)
        show_reactions_by_category(reaction_types)
        
        print_section("ANALYTICS COMPLETE", Fore.GREEN + Style.BRIGHT)
        print(f"\n{Fore.YELLOW + Style.BRIGHT}For more options:{Style.RESET_ALL}")
        print(f"  {Fore.CYAN}python app/reaction_types_analytics.py --search <query>{Style.RESET_ALL}")
        print(f"  {Fore.CYAN}python app/reaction_types_analytics.py --category <category>{Style.RESET_ALL}")
        print(f"  {Fore.CYAN}python app/reaction_types_analytics.py --smarts{Style.RESET_ALL}")
        print(f"  {Fore.CYAN}python app/reaction_types_analytics.py --id <reaction_id>{Style.RESET_ALL}")
        print(f"  {Fore.CYAN}python app/reaction_types_analytics.py --export summary.json{Style.RESET_ALL}")
        
    except Exception as e:
        print(f"\n{Fore.RED + Style.BRIGHT}[X] Error during analysis:{Style.RESET_ALL} {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
