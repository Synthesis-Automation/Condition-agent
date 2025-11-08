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
    
    # With custom parameters
    python app/unified_interactive_cli.py --rxn "..." -k 5 --min-sim 0.5 --type rule

Features:
    - Interactive REPL with command history
    - Real-time parameter adjustment
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


def display_results(
    results: List[RecommendationResult],
    reaction: str,
    show_details: bool = False,
    compact: bool = False
):
    """Display recommendation results."""
    if compact:
        # Compact mode - one line per result
        for r in results:
            icon = "📋" if r.source_type == "protocol" else "📖"
            print(f"{icon} [{r.rank}] {r.name} (sim: {r.similarity:.3f}, family: {r.family})")
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
                        print("Compact mode enabled")
                    elif arg_lower == "off":
                        compact_mode = False
                        print("Compact mode disabled")
                    else:
                        print("Usage: /compact on|off")
                else:
                    compact_mode = not compact_mode
                    print(f"Compact mode {'enabled' if compact_mode else 'disabled'}")
            
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
                print(f"  compact_mode: {compact_mode}\n")
            
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
                display_results(results, reaction_smiles, show_details=show_details, compact=compact_mode)
        
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
                    display_results(results, args.rxn, show_details=args.details, compact=args.compact)
        
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
