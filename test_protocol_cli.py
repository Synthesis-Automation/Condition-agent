#!/usr/bin/env python3
"""
Interactive Protocol Recommendation CLI

This script provides an interactive command-line interface for testing the
protocol recommendation system. Users can input reaction SMILES and get
recommendations with detailed information.

Usage:
    python test_protocol_cli.py
    
    # Or with a reaction SMILES directly
    python test_protocol_cli.py "CCBr.c1ccccc1B(O)O>>CCc1ccccc1"
"""

import sys
import argparse
from typing import Optional

from chemtools.protocol import ProtocolRecommender


def print_separator(char='=', length=70):
    """Print a separator line"""
    print(char * length)


def print_match_details(match: dict, index: int, show_conditions: bool = True):
    """Print detailed information about a protocol match"""
    print(f"\n{index}. Match Score: {match['similarity']:.4f} ({match['similarity']*100:.1f}%)")
    print_separator('-', 70)
    
    # Basic info
    print(f"📄 File: {match['filename']}")
    print(f"🔬 Family: {match['reaction_family']}")
    print(f"📖 Title: {match['source_title']}")
    
    # Source info
    if match['source_journal']:
        print(f"📰 Journal: {match['source_journal']} ({match['source_year']})")
    if match['source_doi']:
        print(f"🔗 DOI: {match['source_doi']}")
    if match['source_url']:
        print(f"🌐 URL: {match['source_url']}")
    
    # Reaction
    print(f"⚗️  Reaction: {match['reaction_smiles']}")
    
    # Tags
    if match['tags']:
        print(f"🏷️  Tags: {', '.join(match['tags'][:10])}")
        if len(match['tags']) > 10:
            print(f"        ... and {len(match['tags']) - 10} more")
    
    # Notes
    if match['notes']:
        notes = match['notes']
        if len(notes) > 150:
            notes = notes[:150] + "..."
        print(f"📝 Notes: {notes}")
    
    # Conditions (if available)
    if show_conditions and 'conditions' in match:
        print_conditions(match['conditions'])


def print_conditions(conditions: dict):
    """Print extracted reaction conditions"""
    print(f"\n🧪 Extracted Conditions:")
    
    if conditions['catalyst']:
        print(f"   Catalyst: {conditions['catalyst']}")
    
    if conditions['ligand']:
        print(f"   Ligand: {conditions['ligand']}")
    
    if conditions['base']:
        print(f"   Base: {conditions['base']}")
    
    if conditions['solvent']:
        print(f"   Solvent: {conditions['solvent']}")
    
    if conditions['additives']:
        print(f"   Additives: {', '.join(conditions['additives'])}")
    
    if conditions['temperature_C'] is not None:
        print(f"   Temperature: {conditions['temperature_C']} °C")
    
    if conditions['time_h'] is not None:
        print(f"   Time: {conditions['time_h']} h")
    
    if conditions['atmosphere']:
        print(f"   Atmosphere: {conditions['atmosphere']}")


def run_interactive_mode(recommender: ProtocolRecommender):
    """Run interactive mode where user can enter multiple reactions"""
    print_separator()
    print("🧪 Interactive Protocol Recommendation CLI")
    print_separator()
    print()
    print("Enter a reaction SMILES to find matching protocols.")
    print("Type 'help' for options, 'quit' to exit.")
    print()
    
    # Default settings
    k = 5
    show_conditions = True
    filter_family = None
    filter_tags = []
    
    while True:
        print_separator('-', 70)
        print()
        
        # Get user input
        try:
            user_input = input("🔬 Reaction SMILES (or command): ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\n\nGoodbye! 👋")
            break
        
        if not user_input:
            continue
        
        # Handle commands
        if user_input.lower() in ('quit', 'exit', 'q'):
            print("\nGoodbye! 👋")
            break
        
        elif user_input.lower() == 'help':
            print_help()
            continue
        
        elif user_input.lower().startswith('set '):
            handle_set_command(user_input, locals())
            continue
        
        elif user_input.lower() == 'settings':
            print_settings(k, show_conditions, filter_family, filter_tags)
            continue
        
        elif user_input.lower() == 'clear':
            filter_family = None
            filter_tags = []
            print("✅ Filters cleared")
            continue
        
        # Treat as reaction SMILES
        reaction_smiles = user_input
        
        # Validate SMILES format
        if '>>' not in reaction_smiles:
            print("❌ Invalid reaction SMILES. Must contain '>>' separator.")
            print("   Example: CCBr.c1ccccc1B(O)O>>CCc1ccccc1")
            continue
        
        # Get recommendations
        print(f"\n🔍 Searching for similar protocols...")
        print(f"   Query: {reaction_smiles}")
        if filter_family:
            print(f"   Family filter: {filter_family}")
        if filter_tags:
            print(f"   Tag filter: {', '.join(filter_tags)}")
        print()
        
        try:
            results = recommender.recommend_with_details(
                reaction_smiles=reaction_smiles,
                k=k,
                reaction_family=filter_family,
                tags=filter_tags if filter_tags else None
            )
            
            # Check if any matches found
            if not results['matches']:
                print("❌ No matching protocols found.")
                print(f"   Searched {results['metadata']['num_candidates']} candidates")
                
                if filter_family or filter_tags:
                    print("\n💡 Tip: Try removing filters with 'clear' command")
                
                continue
            
            # Print summary
            print(f"✅ Found {len(results['matches'])} matches")
            print(f"   (Searched {results['metadata']['num_candidates']} candidates)")
            
            # Print matches
            for i, match in enumerate(results['matches'], 1):
                print_match_details(match, i, show_conditions)
            
            # Ask if user wants to see full protocol
            print()
            print_separator('-', 70)
            print("💡 Tip: Type 'help' for more options")
        
        except Exception as e:
            print(f"❌ Error: {e}")
            if '--debug' in sys.argv:
                import traceback
                traceback.print_exc()


def print_help():
    """Print help message"""
    print()
    print_separator()
    print("📚 Available Commands")
    print_separator()
    print()
    print("BASIC USAGE:")
    print("  <reaction_smiles>     Find protocols for this reaction")
    print("                        Example: CCBr.c1ccccc1B(O)O>>CCc1ccccc1")
    print()
    print("SETTINGS:")
    print("  set k <number>        Set number of results (default: 5)")
    print("  set family <name>     Filter by reaction family")
    print("  set tags <tag1,tag2>  Filter by tags (comma-separated)")
    print("  set conditions on/off Show/hide detailed conditions")
    print("  settings              Show current settings")
    print("  clear                 Clear all filters")
    print()
    print("NAVIGATION:")
    print("  help                  Show this help message")
    print("  quit / exit / q       Exit the program")
    print()
    print("EXAMPLES:")
    print("  set k 3")
    print("  set family Suzuki_Cu_alkyl_halide+aryl_boron")
    print("  set tags suzuki,palladium")
    print("  set conditions off")
    print()


def handle_set_command(command: str, context: dict):
    """Handle 'set' commands"""
    parts = command.split(maxsplit=2)
    
    if len(parts) < 3:
        print("❌ Invalid set command. Usage: set <option> <value>")
        print("   Type 'help' for available options")
        return
    
    option = parts[1].lower()
    value = parts[2]
    
    if option == 'k':
        try:
            new_k = int(value)
            if new_k < 1 or new_k > 50:
                print("❌ k must be between 1 and 50")
                return
            context['k'] = new_k
            print(f"✅ Set k = {new_k}")
        except ValueError:
            print("❌ k must be a number")
    
    elif option == 'family':
        context['filter_family'] = value if value.lower() != 'none' else None
        print(f"✅ Set family filter: {value}")
    
    elif option == 'tags':
        if value.lower() == 'none':
            context['filter_tags'] = []
            print("✅ Cleared tag filter")
        else:
            tags = [t.strip() for t in value.split(',')]
            context['filter_tags'] = tags
            print(f"✅ Set tag filter: {', '.join(tags)}")
    
    elif option == 'conditions':
        if value.lower() in ('on', 'true', 'yes', '1'):
            context['show_conditions'] = True
            print("✅ Condition display: ON")
        elif value.lower() in ('off', 'false', 'no', '0'):
            context['show_conditions'] = False
            print("✅ Condition display: OFF")
        else:
            print("❌ Invalid value. Use: on/off")
    
    else:
        print(f"❌ Unknown option: {option}")
        print("   Type 'help' for available options")


def print_settings(k: int, show_conditions: bool, filter_family: Optional[str], filter_tags: list):
    """Print current settings"""
    print()
    print_separator('-', 70)
    print("⚙️  Current Settings")
    print_separator('-', 70)
    print(f"  Number of results (k): {k}")
    print(f"  Show conditions: {'ON' if show_conditions else 'OFF'}")
    print(f"  Family filter: {filter_family or 'None'}")
    print(f"  Tag filter: {', '.join(filter_tags) if filter_tags else 'None'}")
    print()


def run_single_query(
    reaction_smiles: str,
    k: int = 5,
    family: Optional[str] = None,
    tags: Optional[str] = None,
    show_conditions: bool = True
):
    """Run a single query and print results"""
    print_separator()
    print("🧪 Protocol Recommendation")
    print_separator()
    print()
    
    # Initialize recommender
    print("Loading protocol index...")
    recommender = ProtocolRecommender()
    print(f"✅ Loaded {len(recommender.indexer.records)} protocols")
    print()
    
    # Parse tags
    tag_list = None
    if tags:
        tag_list = [t.strip() for t in tags.split(',')]
    
    # Get recommendations
    print(f"🔍 Searching for: {reaction_smiles}")
    if family:
        print(f"   Family filter: {family}")
    if tag_list:
        print(f"   Tag filter: {', '.join(tag_list)}")
    print()
    
    results = recommender.recommend_with_details(
        reaction_smiles=reaction_smiles,
        k=k,
        reaction_family=family,
        tags=tag_list
    )
    
    # Print results
    if not results['matches']:
        print("❌ No matching protocols found")
        print(f"   Searched {results['metadata']['num_candidates']} candidates")
        return
    
    print(f"✅ Found {len(results['matches'])} matches")
    print(f"   (Searched {results['metadata']['num_candidates']} candidates)")
    
    for i, match in enumerate(results['matches'], 1):
        print_match_details(match, i, show_conditions)
    
    print()
    print_separator()


def main():
    """Main CLI entry point"""
    parser = argparse.ArgumentParser(
        description='Interactive Protocol Recommendation CLI',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive mode
  python test_protocol_cli.py
  
  # Single query
  python test_protocol_cli.py "CCBr.c1ccccc1B(O)O>>CCc1ccccc1"
  
  # With options
  python test_protocol_cli.py "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" -k 3
  python test_protocol_cli.py "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --family Suzuki
  python test_protocol_cli.py "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --tags suzuki,palladium
        """
    )
    
    parser.add_argument(
        'reaction',
        nargs='?',
        help='Reaction SMILES (e.g., "CCBr.c1ccccc1B(O)O>>CCc1ccccc1")'
    )
    parser.add_argument(
        '-k', '--top-k',
        type=int,
        default=5,
        help='Number of top matches to return (default: 5)'
    )
    parser.add_argument(
        '--family',
        help='Filter by reaction family'
    )
    parser.add_argument(
        '--tags',
        help='Filter by tags (comma-separated, e.g., "suzuki,palladium")'
    )
    parser.add_argument(
        '--no-conditions',
        action='store_true',
        help='Hide detailed conditions'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Show debug information on errors'
    )
    
    args = parser.parse_args()
    
    try:
        # Initialize recommender
        recommender = ProtocolRecommender()
        
        # Single query mode
        if args.reaction:
            run_single_query(
                reaction_smiles=args.reaction,
                k=args.top_k,
                family=args.family,
                tags=args.tags,
                show_conditions=not args.no_conditions
            )
        
        # Interactive mode
        else:
            run_interactive_mode(recommender)
    
    except FileNotFoundError as e:
        print()
        print("❌ Protocol index not found!")
        print()
        print("Please build the index first:")
        print("  python -m chemtools.protocol.cli build")
        print()
        if args.debug:
            import traceback
            traceback.print_exc()
        return 1
    
    except KeyboardInterrupt:
        print("\n\nInterrupted by user. Goodbye! 👋")
        return 0
    
    except Exception as e:
        print()
        print(f"❌ Error: {e}")
        print()
        if args.debug:
            import traceback
            traceback.print_exc()
        return 1
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
