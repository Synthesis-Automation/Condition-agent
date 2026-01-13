"""
Demo script for reagent database validation.

Usage:
    python scripts/validate_reagent_db.py
    python scripts/validate_reagent_db.py --verbose
    python scripts/validate_reagent_db.py --role ligand
    python scripts/validate_reagent_db.py --registry-dir path/to/reagent_db
"""

import sys
from pathlib import Path

# Add parent to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from chemtools.reagent import (
    validate_database,
    validate_entry,
    print_validation_summary,
    print_critical_errors_summary,
)

# Default to data/reagent_db in project root
DEFAULT_REGISTRY_PATH = project_root / "data" / "reagent_db"


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Validate reagent database")
    parser.add_argument(
        "--registry-dir",
        default=str(DEFAULT_REGISTRY_PATH),
        help=f"Path to reagent registry directory (default: {DEFAULT_REGISTRY_PATH})",
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
    
    # Print critical errors summary (always show, even without --verbose)
    if not args.verbose:
        print_critical_errors_summary(results)
    
    # Return exit code
    if results.get("error") or results.get("invalid_entries", 0) > 0:
        return 1
    return 0


def demo_single_entry_validation():
    """Demo: validate a single entry"""
    from chemtools.reagent import validate_entry

    # Valid entry
    entry = {
        "name": "DavePhos",
        "abbreviation": "DavePhos",
        "cas": "213697-53-1",
        "smiles": "",
        "role_1": "ligand",
        "family_1": "trialkyl_triaryl_phosphines",
        "tag_1": "",
    }

    print("Validating example entry...")
    issues = validate_entry(entry, role="ligand")

    if not issues:
        print("OK: entry is valid.")
    else:
        print("Entry has issues:")
        for issue in issues:
            print(f"  {issue['severity'].upper()}: {issue['field']} - {issue['message']}")

    # Invalid entry (missing fields)
    print("\n\nValidating invalid entry...")
    bad_entry = {
        "name": "Test",
    }

    issues = validate_entry(bad_entry, role="ligand")
    print(f"Found {len(issues)} issues:")
    for issue in issues:
        print(f"  {issue['severity'].upper()}: {issue['field']} - {issue['message']}")


if __name__ == "__main__":
    # Uncomment to see single entry validation demo:
    # demo_single_entry_validation()
    
    sys.exit(main())
