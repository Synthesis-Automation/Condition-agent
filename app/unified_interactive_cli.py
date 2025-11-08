#!/usr/bin/env python
"""
Interactive CLI for UnifiedRecommender Testing
==============================================

Interactive REPL for testing the UnifiedRecommender system with real-time
feedback, adjustable parameters, and detailed result display.

Usage:
    # Interactive mode (default)
    python app/unified_interactive_cli.py
    
    # Direct query
    python app/unified_interactive_cli.py --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    # Split mode - show top protocol and top rule with detailed conditions
    python app/unified_interactive_cli.py --rxn "..." --split
    
    # With custom parameters
    python app/unified_interactive_cli.py --rxn "..." -k 5 --min-sim 0.5 --type rule

Features:
    - Interactive REPL with command history
    - Real-time parameter adjustment
    - Split mode: Shows top protocol + top rule with full detailed conditions
    - Detailed or compact output modes
    - Full source detail loading
    - Statistics display
    - Color-coded output (if available)
"""

import sys
import json
from pathlib import Path
from typing import Optional, List

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools.recommend import UnifiedRecommender, RecommendationResult

# Try to import colorama for colored output
try:
    from colorama import init, Fore, Style
    init(autoreset=True)
    HAS_COLOR = True
except ImportError:
    HAS_COLOR = False
    # Dummy placeholders
    class Fore:
        GREEN = RED = YELLOW = CYAN = MAGENTA = BLUE = ""
    class Style:
        BRIGHT = RESET_ALL = ""


def print_separator(char="=", length=80, color=""):
    """Print a separator line with optional color."""
    if HAS_COLOR and color:
        print(color + char * length + Style.RESET_ALL)
    else:
        print(char * length)


def format_result(result: RecommendationResult, rank: int, show_details: bool = False) -> str:
    """Format a single recommendation result for display."""
    icon = "📋" if result.source_type == "protocol" else "📖"
    
    # Color-code by similarity
    if HAS_COLOR:
        if result.similarity >= 0.9:
            sim_color = Fore.GREEN + Style.BRIGHT
        elif result.similarity >= 0.7:
            sim_color = Fore.GREEN
        elif result.similarity >= 0.5:
            sim_color = Fore.YELLOW
        else:
            sim_color = Fore.RED
    else:
        sim_color = ""
    
    lines = [
        f"{icon} [{rank}] {result.name}",
        f"    Family: {result.family}",
        f"    Similarity: {sim_color}{result.similarity:.3f}{Style.RESET_ALL if HAS_COLOR else ''}",
        f"    Type: {result.source_type}",
        f"    Version: {result.version}",
    ]
    
    if result.tags:
        tags_str = ", ".join(result.tags[:5])
        if len(result.tags) > 5:
            tags_str += f" (+{len(result.tags) - 5} more)"
        lines.append(f"    Tags: {tags_str}")
    
    if show_details:
        lines.append(f"    ID: {result.id}")
        lines.append(f"    Source: {result.source_file}")
    
    return "\n".join(lines)


def display_detailed_conditions(recommender: UnifiedRecommender, result: RecommendationResult):
    """Display detailed recommended conditions for a source."""
    details = recommender.get_source_details(result.id)
    
    if not details:
        print(f"  {Fore.YELLOW if HAS_COLOR else ''}(No detailed conditions available){Style.RESET_ALL if HAS_COLOR else ''}")
        return
    
    # Extract key condition fields based on source type
    if result.source_type == "protocol":
        # Protocol has reaction_setup with chemicals and conditions
        print(f"  {Fore.CYAN if HAS_COLOR else ''}Protocol Conditions:{Style.RESET_ALL if HAS_COLOR else ''}")
        
        if 'reaction_setup' in details and details['reaction_setup']:
            setup = details['reaction_setup'][0]  # First setup step
            
            # Extract chemicals by role
            chemicals = setup.get('chemicals', [])
            catalyst_chems = [c for c in chemicals if 'catalyst' in c.get('role', '').lower()]
            solvent_chems = [c for c in chemicals if 'solvent' in c.get('role', '').lower()]
            base_chems = [c for c in chemicals if 'base' in c.get('role', '').lower()]
            ligand_chems = [c for c in chemicals if 'ligand' in c.get('role', '').lower()]
            
            if catalyst_chems:
                cat_names = ', '.join([c['name'] for c in catalyst_chems])
                print(f"    • Catalyst: {cat_names}")
            
            if ligand_chems:
                lig_names = ', '.join([c['name'] for c in ligand_chems])
                print(f"    • Ligand: {lig_names}")
            
            if solvent_chems:
                solv_names = ', '.join([c['name'] for c in solvent_chems])
                print(f"    • Solvent: {solv_names}")
            
            if base_chems:
                base_names = ', '.join([c['name'] for c in base_chems])
                print(f"    • Base: {base_names}")
            
            # Extract conditions (temperature, time, atmosphere)
            conditions = setup.get('conditions', [])
            if conditions:
                # Show main reaction conditions (usually last step)
                main_cond = conditions[-1] if conditions else {}
                
                if 'temperature_C' in main_cond:
                    print(f"    • Temperature: {main_cond['temperature_C']} °C")
                if 'time_h' in main_cond:
                    print(f"    • Time: {main_cond['time_h']} h")
                if 'atmosphere' in main_cond:
                    atm = main_cond['atmosphere']
                    # Simplify atmosphere display
                    if 'N2' in atm or 'nitrogen' in atm.lower():
                        atm_simple = 'N₂ (inert atmosphere)'
                    elif 'Ar' in atm or 'argon' in atm.lower():
                        atm_simple = 'Ar (inert atmosphere)'
                    else:
                        atm_simple = atm
                    print(f"    • Atmosphere: {atm_simple}")
        
        # Show yield if available in metadata
        if 'metadata' in details and 'yield_percent' in details['metadata']:
            print(f"    • Reported Yield: {details['metadata']['yield_percent']}%")
    
    elif result.source_type == "rule":
        # Rule has default_rule and base_rules with conditions
        # Show default rule conditions
        if 'default_rule' in details and 'conditions' in details['default_rule']:
            cond = details['default_rule']['conditions']
            rule_desc = details['default_rule'].get('description', 'Default conditions')
            print(f"  {Fore.CYAN if HAS_COLOR else ''}Default Conditions:{Style.RESET_ALL if HAS_COLOR else ''}")
            print(f"    {Fore.YELLOW if HAS_COLOR else ''}({rule_desc}){Style.RESET_ALL if HAS_COLOR else ''}")
            
            if 'pd_source' in cond and cond['pd_source']:
                print(f"    • Pd Source: {cond['pd_source']}")
            if 'precatalyst' in cond and cond['precatalyst']:
                print(f"    • Precatalyst: {cond['precatalyst']}")
            if 'catalyst_loading_molpct' in cond and cond['catalyst_loading_molpct']:
                print(f"    • Catalyst Loading: {cond['catalyst_loading_molpct']} mol%")
            if 'ligand' in cond and cond['ligand']:
                print(f"    • Ligand: {cond['ligand']}")
            if 'solvent' in cond and cond['solvent']:
                print(f"    • Solvent: {cond['solvent']}")
            if 'base' in cond and cond['base']:
                print(f"    • Base: {cond['base']}")
            if 'base_equiv' in cond and cond['base_equiv']:
                print(f"    • Base Equivalents: {cond['base_equiv']}")
            if 'temperature_C' in cond and cond['temperature_C']:
                print(f"    • Temperature: {cond['temperature_C']} °C")
            if 'time_h' in cond and cond['time_h']:
                print(f"    • Time: {cond['time_h']} h")
            if 'atmosphere' in cond and cond['atmosphere']:
                print(f"    • Atmosphere: {cond['atmosphere']}")
        
        # Show alternative base rules (top 3)
        if 'base_rules' in details and details['base_rules']:
            print(f"\n  {Fore.MAGENTA if HAS_COLOR else ''}Alternative Conditions:{Style.RESET_ALL if HAS_COLOR else ''}")
            for i, rule in enumerate(details['base_rules'][:3], 1):
                rule_name = rule.get('name', f'Alternative {i}')
                rule_desc = rule.get('description', '')
                print(f"    {i}. {Fore.YELLOW if HAS_COLOR else ''}{rule_name}{Style.RESET_ALL if HAS_COLOR else ''}")
                if rule_desc:
                    print(f"       {rule_desc}")
                
                if 'conditions' in rule:
                    cond = rule['conditions']
                    # Show key differences from default
                    key_conditions = []
                    if 'pd_source' in cond:
                        key_conditions.append(f"Pd: {cond['pd_source']}")
                    if 'precatalyst' in cond:
                        key_conditions.append(f"Catalyst: {cond['precatalyst']}")
                    if 'ligand' in cond:
                        key_conditions.append(f"Ligand: {cond['ligand']}")
                    if 'base' in cond:
                        key_conditions.append(f"Base: {cond['base']}")
                    if 'temperature_C' in cond:
                        key_conditions.append(f"Temp: {cond['temperature_C']} °C")
                    
                    if key_conditions:
                        print(f"       {' | '.join(key_conditions)}")
                print()


def display_results(
    results: List[RecommendationResult],
    reaction: str,
    show_details: bool = False,
    compact: bool = False,
    show_top_split: bool = False,
    recommender: Optional[UnifiedRecommender] = None
):
    """Display recommendation results."""
    if compact:
        # Compact mode - one line per result
        for r in results:
            icon = "📋" if r.source_type == "protocol" else "📖"
            print(f"{icon} [{r.rank}] {r.name} (sim: {r.similarity:.3f}, family: {r.family})")
    
    elif show_top_split and recommender:
        # Show top protocol and top rule separately with detailed conditions
        print()
        print_separator(color=Fore.CYAN if HAS_COLOR else "")
        print(f"Query: {reaction}")
        print_separator(color=Fore.CYAN if HAS_COLOR else "")
        print()
        
        # Find top protocol and top rule
        top_protocol = None
        top_rule = None
        
        for r in results:
            if r.source_type == "protocol" and top_protocol is None:
                top_protocol = r
            elif r.source_type == "rule" and top_rule is None:
                top_rule = r
            
            if top_protocol and top_rule:
                break
        
        # Display top protocol
        if top_protocol:
            print(f"{Fore.GREEN if HAS_COLOR else ''}━━━ TOP PROTOCOL ━━━{Style.RESET_ALL if HAS_COLOR else ''}")
            print(f"📋 {Fore.YELLOW if HAS_COLOR else ''}{top_protocol.name}{Style.RESET_ALL if HAS_COLOR else ''}")
            print(f"   Similarity: {Fore.GREEN if HAS_COLOR else ''}{top_protocol.similarity:.3f}{Style.RESET_ALL if HAS_COLOR else ''}")
            print(f"   Family: {top_protocol.family}")
            print(f"   ID: {top_protocol.id}")
            print()
            display_detailed_conditions(recommender, top_protocol)
            print()
        else:
            print(f"{Fore.YELLOW if HAS_COLOR else ''}No protocol recommendations found.{Style.RESET_ALL if HAS_COLOR else ''}")
            print()
        
        # Display top rule
        if top_rule:
            print(f"{Fore.GREEN if HAS_COLOR else ''}━━━ TOP RULE ━━━{Style.RESET_ALL if HAS_COLOR else ''}")
            print(f"📖 {Fore.YELLOW if HAS_COLOR else ''}{top_rule.name}{Style.RESET_ALL if HAS_COLOR else ''}")
            print(f"   Similarity: {Fore.GREEN if HAS_COLOR else ''}{top_rule.similarity:.3f}{Style.RESET_ALL if HAS_COLOR else ''}")
            print(f"   Family: {top_rule.family}")
            print(f"   ID: {top_rule.id}")
            print()
            display_detailed_conditions(recommender, top_rule)
            print()
        else:
            print(f"{Fore.YELLOW if HAS_COLOR else ''}No rule recommendations found.{Style.RESET_ALL if HAS_COLOR else ''}")
            print()
        
        print_separator(color=Fore.CYAN if HAS_COLOR else "")
    
    else:
        # Full mode
        print()
        print_separator(color=Fore.CYAN if HAS_COLOR else "")
        print(f"Query: {reaction}")
        print(f"Found {len(results)} recommendation(s)")
        print_separator(color=Fore.CYAN if HAS_COLOR else "")
        print()
        
        for r in results:
            print(format_result(r, r.rank, show_details=show_details))
            print()


def display_statistics(recommender: UnifiedRecommender):
    """Display index statistics."""
    stats = recommender.get_statistics()
    
    print_separator(color=Fore.MAGENTA if HAS_COLOR else "")
    print("Index Statistics")
    print_separator(color=Fore.MAGENTA if HAS_COLOR else "")
    print(f"Version: {stats['build_info']['version']}")
    print(f"Build Date: {stats['build_info']['build_date']}")
    print()
    print(f"Protocols: {stats['protocols']['count']}")
    print(f"Rules: {stats['rules']['count']}")
    print(f"Total Fingerprints: {stats['drfp']['computed']}")
    print(f"Failed Fingerprints: {stats['drfp']['failed']}")
    print()
    
    # Show protocol families
    if 'families' in stats['protocols']:
        families = stats['protocols']['families']
        print(f"Protocol Families ({len(families)}):")
        for family, count in sorted(families.items()):
            print(f"  - {family}: {count}")
        print()
    
    # Show rule families
    if 'families' in stats['rules']:
        families = stats['rules']['families']
        print(f"Rule Families ({len(families)}):")
        for family, count in sorted(families.items()):
            print(f"  - {family}: {count}")
    
    print_separator(color=Fore.MAGENTA if HAS_COLOR else "")


def load_and_display_source(recommender: UnifiedRecommender, source_id: str):
    """Load and display full source details."""
    details = recommender.get_source_details(source_id)
    
    if not details:
        print(f"{Fore.RED if HAS_COLOR else ''}Error: Source '{source_id}' not found{Style.RESET_ALL if HAS_COLOR else ''}")
        return
    
    print_separator(color=Fore.CYAN if HAS_COLOR else "")
    print(f"Source Details: {source_id}")
    print_separator(color=Fore.CYAN if HAS_COLOR else "")
    
    # Pretty print JSON
    print(json.dumps(details, indent=2))
    print_separator(color=Fore.CYAN if HAS_COLOR else "")


def interactive_mode(
    recommender: UnifiedRecommender,
    initial_k: int = 5,
    initial_min_sim: float = 0.0,
    initial_type: Optional[str] = None
):
    """Interactive REPL mode."""
    print_separator("=", color=Fore.CYAN if HAS_COLOR else "")
    print(f"{Fore.CYAN if HAS_COLOR else ''}UnifiedRecommender Interactive Mode{Style.RESET_ALL if HAS_COLOR else ''}")
    print_separator("=", color=Fore.CYAN if HAS_COLOR else "")
    print()
    print("Enter a reaction SMILES (reactants>>products) to get recommendations.")
    print()
    print(f"{Fore.YELLOW if HAS_COLOR else ''}Commands:{Style.RESET_ALL if HAS_COLOR else ''}")
    print("  /k <int>              -> set top_k (current: %d)" % initial_k)
    print("  /sim <float>          -> set min_similarity 0.0-1.0 (current: %.2f)" % initial_min_sim)
    print("  /type <protocol|rule|all> -> filter by source type (current: %s)" % (initial_type or "all"))
    print("  /details on|off       -> toggle detailed output (current: off)")
    print("  /compact on|off       -> toggle compact mode (current: off)")
    print("  /split on|off         -> toggle top protocol + top rule split view (current: off)")
    print("  /stats                -> show index statistics")
    print("  /load <id>            -> load full source details by ID")
    print("  /show                 -> display current settings")
    print("  /help                 -> show this message")
    print("  /exit or /quit        -> exit interactive mode")
    print_separator("=", color=Fore.CYAN if HAS_COLOR else "")
    print()
    
    # Current settings
    k = initial_k
    min_sim = initial_min_sim
    source_type = initial_type
    show_details = False
    compact_mode = False
    split_mode = False
    
    while True:
        try:
            raw = input(f"{Fore.GREEN if HAS_COLOR else ''}reaction> {Style.RESET_ALL if HAS_COLOR else ''}").strip()
        except (KeyboardInterrupt, EOFError):
            print("\nExiting.")
            break
        
        if not raw:
            continue
        
        # Handle commands
        if raw.startswith("/"):
            cmd_parts = raw.split(maxsplit=1)
            cmd = cmd_parts[0].lower()
            arg = cmd_parts[1] if len(cmd_parts) > 1 else None
            
            if cmd in ["/exit", "/quit"]:
                print("Exiting.")
                break
            
            elif cmd == "/help":
                print("\nCommands:")
                print("  /k <int>              -> set top_k")
                print("  /sim <float>          -> set min_similarity")
                print("  /type <protocol|rule|all> -> filter by source type")
                print("  /details on|off       -> toggle detailed output")
                print("  /compact on|off       -> toggle compact mode")
                print("  /split on|off         -> toggle top protocol + top rule split view")
                print("  /stats                -> show index statistics")
                print("  /load <id>            -> load full source details")
                print("  /show                 -> display current settings")
                print("  /help                 -> show this message")
                print("  /exit or /quit        -> exit\n")
            
            elif cmd == "/k":
                if arg and arg.isdigit():
                    k = int(arg)
                    print(f"Set top_k = {k}")
                else:
                    print("Usage: /k <integer>")
            
            elif cmd == "/sim":
                if arg:
                    try:
                        min_sim = float(arg)
                        if 0.0 <= min_sim <= 1.0:
                            print(f"Set min_similarity = {min_sim}")
                        else:
                            print("Error: min_similarity must be between 0.0 and 1.0")
                            min_sim = 0.0
                    except ValueError:
                        print("Usage: /sim <float>")
                else:
                    print("Usage: /sim <float>")
            
            elif cmd == "/type":
                if arg:
                    arg_lower = arg.lower()
                    if arg_lower in ["protocol", "rule"]:
                        source_type = arg_lower
                        print(f"Set source_type = {source_type}")
                    elif arg_lower == "all":
                        source_type = None
                        print("Set source_type = all")
                    else:
                        print("Usage: /type <protocol|rule|all>")
                else:
                    print("Usage: /type <protocol|rule|all>")
            
            elif cmd == "/details":
                if arg:
                    arg_lower = arg.lower()
                    if arg_lower == "on":
                        show_details = True
                        print("Detailed output enabled")
                    elif arg_lower == "off":
                        show_details = False
                        print("Detailed output disabled")
                    else:
                        print("Usage: /details on|off")
                else:
                    show_details = not show_details
                    print(f"Detailed output {'enabled' if show_details else 'disabled'}")
            
            elif cmd == "/compact":
                if arg:
                    arg_lower = arg.lower()
                    if arg_lower == "on":
                        compact_mode = True
                        split_mode = False  # Disable split mode
                        print("Compact mode enabled")
                    elif arg_lower == "off":
                        compact_mode = False
                        print("Compact mode disabled")
                    else:
                        print("Usage: /compact on|off")
                else:
                    compact_mode = not compact_mode
                    if compact_mode:
                        split_mode = False
                    print(f"Compact mode {'enabled' if compact_mode else 'disabled'}")
            
            elif cmd == "/split":
                if arg:
                    arg_lower = arg.lower()
                    if arg_lower == "on":
                        split_mode = True
                        compact_mode = False  # Disable compact mode
                        print("Split mode enabled - will show top protocol and top rule with detailed conditions")
                    elif arg_lower == "off":
                        split_mode = False
                        print("Split mode disabled")
                    else:
                        print("Usage: /split on|off")
                else:
                    split_mode = not split_mode
                    if split_mode:
                        compact_mode = False
                    print(f"Split mode {'enabled' if split_mode else 'disabled'}")
            
            elif cmd == "/stats":
                display_statistics(recommender)
            
            elif cmd == "/load":
                if arg:
                    load_and_display_source(recommender, arg)
                else:
                    print("Usage: /load <source_id>")
            
            elif cmd == "/show":
                print(f"\nCurrent Settings:")
                print(f"  top_k: {k}")
                print(f"  min_similarity: {min_sim}")
                print(f"  source_type: {source_type or 'all'}")
                print(f"  show_details: {show_details}")
                print(f"  compact_mode: {compact_mode}")
                print(f"  split_mode: {split_mode}\n")
            
            else:
                print(f"Unknown command: {cmd}. Type /help for available commands.")
            
            continue
        
        # Process as reaction SMILES
        reaction_smiles = raw
        
        # Validate basic format
        if ">>" not in reaction_smiles:
            print(f"{Fore.RED if HAS_COLOR else ''}Error: Invalid reaction SMILES (missing >>){Style.RESET_ALL if HAS_COLOR else ''}")
            continue
        
        # Get recommendations
        try:
            source_types = [source_type] if source_type else None
            results = recommender.recommend(
                reaction_smiles=reaction_smiles,
                top_k=k,
                min_similarity=min_sim,
                source_types=source_types
            )
            
            if not results:
                print(f"\n{Fore.YELLOW if HAS_COLOR else ''}No recommendations found matching criteria.{Style.RESET_ALL if HAS_COLOR else ''}")
                print("Try lowering min_similarity or changing source_type filter.\n")
            else:
                display_results(
                    results, 
                    reaction_smiles, 
                    show_details=show_details, 
                    compact=compact_mode,
                    show_top_split=split_mode,
                    recommender=recommender
                )
        
        except Exception as e:
            print(f"{Fore.RED if HAS_COLOR else ''}Error: {str(e)}{Style.RESET_ALL if HAS_COLOR else ''}")


def main():
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Interactive CLI for UnifiedRecommender testing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive mode
  %(prog)s
  
  # Direct query
  %(prog)s --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
  
  # With filters
  %(prog)s --rxn "..." -k 3 --min-sim 0.5 --type protocol
  
  # Show statistics and exit
  %(prog)s --stats
        """
    )
    
    parser.add_argument(
        "--rxn", "--reaction",
        help="Reaction SMILES (reactants>>products). If not provided, enters interactive mode."
    )
    parser.add_argument(
        "-k", "--top-k",
        type=int,
        default=5,
        help="Number of recommendations to return (default: 5)"
    )
    parser.add_argument(
        "--min-sim",
        type=float,
        default=0.0,
        help="Minimum similarity threshold 0.0-1.0 (default: 0.0)"
    )
    parser.add_argument(
        "--type",
        choices=["protocol", "rule", "all"],
        default="all",
        help="Filter by source type (default: all)"
    )
    parser.add_argument(
        "--index-dir",
        help="Path to custom unified index directory (optional)"
    )
    parser.add_argument(
        "--details",
        action="store_true",
        help="Show detailed output (IDs, file paths)"
    )
    parser.add_argument(
        "--compact",
        action="store_true",
        help="Use compact output mode"
    )
    parser.add_argument(
        "--split",
        action="store_true",
        help="Show top protocol and top rule separately with detailed conditions"
    )
    parser.add_argument(
        "--stats",
        action="store_true",
        help="Show index statistics and exit"
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output results as JSON (non-interactive only)"
    )
    
    args = parser.parse_args()
    
    # Initialize recommender
    try:
        if args.index_dir:
            recommender = UnifiedRecommender(index_dir=args.index_dir)
        else:
            recommender = UnifiedRecommender()
    except Exception as e:
        print(f"Error loading recommender: {e}", file=sys.stderr)
        return 1
    
    # Stats mode
    if args.stats:
        display_statistics(recommender)
        return 0
    
    # Direct query mode
    if args.rxn:
        source_types = [args.type] if args.type != "all" else None
        
        try:
            results = recommender.recommend(
                reaction_smiles=args.rxn,
                top_k=args.top_k,
                min_similarity=args.min_sim,
                source_types=source_types
            )
            
            if args.json:
                # JSON output
                output = {
                    "query": args.rxn,
                    "count": len(results),
                    "filters": {
                        "top_k": args.top_k,
                        "min_similarity": args.min_sim,
                        "source_type": args.type
                    },
                    "recommendations": [
                        {
                            "rank": r.rank,
                            "id": r.id,
                            "name": r.name,
                            "source_type": r.source_type,
                            "family": r.family,
                            "similarity": round(r.similarity, 3),
                            "tags": r.tags,
                            "version": r.version,
                        }
                        for r in results
                    ]
                }
                print(json.dumps(output, indent=2))
            else:
                # Human-readable output
                if not results:
                    print(f"{Fore.YELLOW if HAS_COLOR else ''}No recommendations found matching criteria.{Style.RESET_ALL if HAS_COLOR else ''}")
                else:
                    display_results(
                        results, 
                        args.rxn, 
                        show_details=args.details, 
                        compact=args.compact,
                        show_top_split=args.split,
                        recommender=recommender
                    )
        
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1
        
        return 0
    
    # Interactive mode
    source_type = args.type if args.type != "all" else None
    interactive_mode(
        recommender,
        initial_k=args.top_k,
        initial_min_sim=args.min_sim,
        initial_type=source_type
    )
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
