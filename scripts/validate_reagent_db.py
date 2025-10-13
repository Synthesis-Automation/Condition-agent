"""
Demo script for reagent database validation.

Usage:
    python scripts/validate_reagent_db.py
    python scripts/validate_reagent_db.py --verbose
    python scripts/validate_reagent_db.py --role ligand
"""

import sys
from pathlib import Path

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools.reagent import (
    validate_database,
    validate_entry,
    print_validation_summary,
    DEFAULT_REGISTRY_DIR,
)


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Validate reagent database")
    parser.add_argument(
        "--registry-dir",
        default=DEFAULT_REGISTRY_DIR,
        help="Path to reagent registry directory",
    )
    parser.add_argument(
        "--role",
        help="Validate specific role only (e.g., ligand, base)",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Show detailed issues for each entry",
    )
    parser.add_argument(
        "--strict",
        "-s",
        action="store_true",
        help="Strict type checking",
    )
    
    args = parser.parse_args()
    
    # Run validation
    roles = [args.role] if args.role else None
    
    print(f"Validating reagent database at: {args.registry_dir}")
    if roles:
        print(f"Role filter: {roles}")
    
    results = validate_database(
        args.registry_dir,
        strict=args.strict,
        roles=roles,
    )
    
    # Print results
    print_validation_summary(results, verbose=args.verbose)
    
    # Return exit code
    if results.get("error") or results.get("invalid_entries", 0) > 0:
        return 1
    return 0


def demo_single_entry_validation():
    """Demo: validate a single entry"""
    from chemtools.reagent import validate_entry
    
    # Valid entry
    entry = {
        "id": "ZEMZPXWZVTUONV-UHFFFAOYSA-N",
        "name": "DavePhos",
        "abbreviation": ["DavePhos"],
        "aliases": ["2-Dicyclohexylphosphino-2'-(dimethylamino)biphenyl"],
        "cas": "213697-53-1",
        "inchi_key": "ZEMZPXWZVTUONV-UHFFFAOYSA-N",
        "smiles": None,
        "roles": {
            "ligand": {
                "families": ["phosphine"],
                "donors": ["P"],
                "denticity": 1,
            }
        },
        "embedding_text": "type: LIGAND | name: DavePhos",
    }
    
    print("Validating example entry...")
    issues = validate_entry(entry, role="ligand")
    
    if not issues:
        print("✅ Entry is valid!")
    else:
        print("❌ Entry has issues:")
        for issue in issues:
            print(f"  {issue['severity'].upper()}: {issue['field']} - {issue['message']}")
    
    # Invalid entry (missing fields)
    print("\n\nValidating invalid entry...")
    bad_entry = {
        "name": "Test",
        "roles": {"ligand": {}},  # Missing required fields
    }
    
    issues = validate_entry(bad_entry, role="ligand")
    print(f"Found {len(issues)} issues:")
    for issue in issues:
        print(f"  {issue['severity'].upper()}: {issue['field']} - {issue['message']}")


if __name__ == "__main__":
    # Uncomment to see single entry validation demo:
    # demo_single_entry_validation()
    
    sys.exit(main())
