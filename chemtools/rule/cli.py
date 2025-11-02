#!/usr/bin/env python
"""
Rule-Based Recommendation CLI
==============================

Command-line interface for testing the rule-based recommendation engine.

Usage:
    python -m chemtools.rule.cli <reaction_smiles> [options]
    python -m chemtools.rule.cli --file <reactions.txt> [options]
    python -m chemtools.rule.cli --interactive [options]

Examples:
    # Single reaction
    python -m chemtools.rule.cli "Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1"
    
    # With symptoms
    python -m chemtools.rule.cli "Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1" --symptoms low_yield
    
    # Batch from file
    python -m chemtools.rule.cli --file reactions.txt --database data/suzuki.json
    
    # Interactive mode
    python -m chemtools.rule.cli --interactive
"""

import argparse
import sys
from pathlib import Path
from typing import List, Optional
import logging

from .engine import RuleEngine
from .models import ConditionRecommendation

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s'
)
logger = logging.getLogger(__name__)


def find_database(name: str = "suzuki") -> Path:
    """
    Find a database file by name.
    
    Args:
        name: Database name (without .json extension)
    
    Returns:
        Path to database file
    """
    # Try common locations
    locations = [
        Path(f"data/rule_db/{name}.json"),
        Path(f"data/{name}.json"),
        Path(f"data/protocol_db/{name}.json"),
        Path(f"{name}.json"),
        Path(__file__).parent.parent.parent / "data" / "rule_db" / f"{name}.json",
        Path(__file__).parent.parent.parent / "data" / f"{name}.json",
        Path(__file__).parent.parent.parent / "data" / "protocol_db" / f"{name}.json"
    ]
    
    for path in locations:
        if path.exists():
            return path
    
    raise FileNotFoundError(
        f"Database '{name}' not found. Searched:\n" +
        "\n".join(f"  - {p}" for p in locations)
    )


def process_reaction(
    engine: RuleEngine,
    reaction_smiles: str,
    symptoms: Optional[List[str]] = None,
    verbose: bool = False
) -> None:
    """
    Process a single reaction and print recommendation.
    
    Args:
        engine: RuleEngine instance
        reaction_smiles: Reaction SMILES string
        symptoms: Optional list of symptoms
        verbose: If True, print detailed output
    """
    try:
        # Generate recommendation
        rec = engine.recommend(
            reaction_smiles,
            symptoms=symptoms,
            include_reasoning=verbose
        )
        
        # Print formatted output
        if verbose:
            print("\n" + "=" * 70)
            print(rec.format_summary())
            print("=" * 70)
            
            # Also print JSON for debugging
            print("\nJSON Output:")
            import json
            print(json.dumps(rec.to_dict(), indent=2))
        else:
            print(rec.format_summary())
    
    except Exception as e:
        logger.error(f"Error processing reaction: {e}")
        if verbose:
            import traceback
            traceback.print_exc()


def interactive_mode(engine: RuleEngine, verbose: bool = False) -> None:
    """
    Run in interactive mode.
    
    Args:
        engine: RuleEngine instance
        verbose: If True, print detailed output
    """
    print("\n" + "=" * 70)
    print("Rule-Based Recommendation Interactive Mode")
    print("=" * 70)
    print(f"\nDatabase: {engine.get_database_summary()}\n")
    print("Commands:")
    print("  <reaction_smiles>          - Get recommendation")
    print("  symptoms: <symptom_list>   - Add symptoms (comma-separated)")
    print("  verbose                    - Toggle verbose mode")
    print("  help                       - Show this help")
    print("  quit                       - Exit")
    print("=" * 70 + "\n")
    
    current_symptoms = []
    
    while True:
        try:
            user_input = input(">>> ").strip()
            
            if not user_input:
                continue
            
            # Handle commands
            if user_input.lower() in ["quit", "exit", "q"]:
                print("Goodbye!")
                break
            
            elif user_input.lower() == "help":
                print("\nCommands:")
                print("  <reaction_smiles>          - Get recommendation")
                print("  symptoms: <symptom_list>   - Add symptoms (comma-separated)")
                print("  verbose                    - Toggle verbose mode")
                print("  help                       - Show this help")
                print("  quit                       - Exit\n")
            
            elif user_input.lower() == "verbose":
                verbose = not verbose
                print(f"Verbose mode: {'ON' if verbose else 'OFF'}")
            
            elif user_input.lower().startswith("symptoms:"):
                symptom_str = user_input[9:].strip()
                current_symptoms = [s.strip() for s in symptom_str.split(",") if s.strip()]
                print(f"Current symptoms: {current_symptoms}")
            
            else:
                # Treat as reaction SMILES
                process_reaction(engine, user_input, current_symptoms, verbose)
        
        except KeyboardInterrupt:
            print("\nGoodbye!")
            break
        except EOFError:
            print("\nGoodbye!")
            break
        except Exception as e:
            logger.error(f"Error: {e}")
            if verbose:
                import traceback
                traceback.print_exc()


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Rule-based condition recommendation engine",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    # Input options
    parser.add_argument(
        "reaction",
        nargs="?",
        help="Reaction SMILES string"
    )
    parser.add_argument(
        "--file", "-f",
        type=Path,
        help="File with reaction SMILES (one per line)"
    )
    parser.add_argument(
        "--interactive", "-i",
        action="store_true",
        help="Run in interactive mode"
    )
    
    # Database options
    parser.add_argument(
        "--database", "-d",
        type=str,
        default="suzuki",
        help="Database name or path (default: suzuki)"
    )
    
    # Processing options
    parser.add_argument(
        "--symptoms", "-s",
        nargs="+",
        help="Observed symptoms (e.g., low_yield side_products)"
    )
    parser.add_argument(
        "--combine",
        choices=["union", "all", "first", "separate"],
        default="union",
        help="How to combine features from multiple reactants (default: union)"
    )
    
    # Output options
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Print detailed output with features and reasoning"
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output in JSON format"
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Validate database and exit"
    )
    
    args = parser.parse_args()
    
    # Load database
    try:
        if args.database.endswith(".json"):
            db_path = Path(args.database)
        else:
            db_path = find_database(args.database)
        
        logger.info(f"Loading database from {db_path}")
        engine = RuleEngine.from_file(db_path)
        
        # Validate if requested
        if args.validate:
            issues = engine.validate_database()
            if issues:
                print("Validation issues:")
                for issue in issues:
                    print(f"  - {issue}")
                sys.exit(1)
            else:
                print("Database is valid")
                print(engine.get_database_summary())
                sys.exit(0)
    
    except Exception as e:
        logger.error(f"Failed to load database: {e}")
        sys.exit(1)
    
    # Interactive mode
    if args.interactive:
        interactive_mode(engine, args.verbose)
        return
    
    # Batch file mode
    if args.file:
        try:
            with open(args.file, 'r') as f:
                reactions = [line.strip() for line in f if line.strip()]
            
            logger.info(f"Processing {len(reactions)} reactions from {args.file}")
            
            for i, reaction in enumerate(reactions, 1):
                print(f"\n--- Reaction {i}/{len(reactions)} ---")
                process_reaction(engine, reaction, args.symptoms, args.verbose)
        
        except Exception as e:
            logger.error(f"Error reading file: {e}")
            sys.exit(1)
        
        return
    
    # Single reaction mode
    if args.reaction:
        process_reaction(engine, args.reaction, args.symptoms, args.verbose)
        return
    
    # No input provided
    parser.print_help()
    sys.exit(1)


if __name__ == "__main__":
    main()
