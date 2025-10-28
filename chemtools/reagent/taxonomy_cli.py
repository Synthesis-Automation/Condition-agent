#!/usr/bin/env python3

"""
Command-line reagent taxonomy generator.

Usage:
    python -m chemtools.reagent.taxonomy_cli --cas "14221-01-3" --name "Pd(PPh3)4"
    python -m chemtools.reagent.taxonomy_cli --list-families
    python -m chemtools.reagent.taxonomy_cli --cas "14221-01-3"  # Auto-resolve name
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .constants import ROLE_FILES
from .registry_addition import (
    ReagentAdditionError,
    add_reagent_entry,
    list_available_families,
)

# Ensure UTF-8 encoding for output
for _stream in ("stdout", "stderr"):
    try:
        getattr(sys, _stream).reconfigure(encoding="utf-8")
    except Exception:
        pass


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Assign reagent role/family and update taxonomy.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Add reagent with auto-detection
  %(prog)s --cas "14221-01-3" --name "Pd(PPh3)4"
  
  # Auto-resolve name from CAS (requires internet)
  %(prog)s --cas "14221-01-3"
  
  # Specify role and family explicitly
  %(prog)s --cas "14221-01-3" --name "Pd(PPh3)4" \\
           --role metal_precursor --family pd_0_complexes
  
  # List all available families
  %(prog)s --list-families
  
  # Dry run (preview without saving)
  %(prog)s --cas "14221-01-3" --name "Pd(PPh3)4" --dry-run
        """
    )
    parser.add_argument("--cas", help="CAS number for the reagent.")
    parser.add_argument("--name", help="Primary name for the reagent.")
    parser.add_argument("--abbr", help="Short abbreviation.")
    parser.add_argument(
        "--synonym",
        action="append",
        default=[],
        help="Additional synonym (repeatable)."
    )
    parser.add_argument(
        "--role",
        choices=sorted(ROLE_FILES.keys()),
        help="Override detected role."
    )
    parser.add_argument("--family", help="Explicit family identifier to use.")
    parser.add_argument("--smiles", help="Optional SMILES string to store with the entry.")
    parser.add_argument(
        "--taxonomy-dir",
        default="data/compound_taxonomy",
        help="Path to taxonomy directory."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Infer role/family but do not write files."
    )
    parser.add_argument(
        "--list-families",
        action="store_true",
        help="Print known families and exit."
    )
    parser.add_argument(
        "--allow-default-family",
        action="store_true",
        help="Allow fallback default family when inference fails."
    )
    parser.add_argument(
        "--no-auto-resolve",
        action="store_true",
        help="Disable CAS-based identity lookup when name is missing."
    )
    parser.add_argument(
        "--resolver-timeout",
        type=float,
        default=DEFAULT_RESOLVER_TIMEOUT,
        help="Timeout (seconds) for resolver HTTP requests."
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Emit extra debug details."
    )
    return parser.parse_args()


def main() -> None:
    """Main CLI entry point."""
    args = parse_args()
    taxonomy_path = Path(args.taxonomy_dir).expanduser()

    if args.list_families:
        families = list_available_families(taxonomy_path)
        print(json.dumps(families, indent=2, ensure_ascii=False))
        return

    if not args.cas:
        raise SystemExit("--cas is required unless --list-families is used.")

    try:
        result = add_reagent_entry(
            cas=args.cas,
            taxonomy_dir=taxonomy_path,
            name=args.name,
            synonyms=args.synonym,
            abbreviation=args.abbr,
            role=args.role,
            family_id=args.family,
            smiles=args.smiles,
            allow_default_family=args.allow_default_family,
            dry_run=args.dry_run,
            auto_resolve=not args.no_auto_resolve,
            resolver_timeout=args.resolver_timeout,
        )
    except ReagentAdditionError as exc:
        raise SystemExit(str(exc))

    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
