"""
Simple validation entrypoint for taxonomy data files.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

from . import load_registry, reset_registry


def validate(root: Optional[Path] = None) -> None:
    """
    Load the taxonomy registry from ``root`` (or default location) and print a summary.
    """
    registry = load_registry(root)
    print("Taxonomy validation succeeded:")
    print(f"  Root: {registry.root}")
    print(f"  Schema version: {registry.manifest.schema_version}")
    print(f"  Taxonomy version: {registry.manifest.taxonomy_version}")
    print(f"  Reaction categories: {len(list(registry.iter_reaction_categories()))}")
    print(f"  Reaction types: {len(list(registry.iter_reaction_types()))}")
    print(f"  Reactant types: {len(list(registry.iter_reactant_types()))}")
    print(f"  Reagent roles: {len(list(registry.iter_reagent_roles()))}")
    print(f"  Reagent families: {len(list(registry.iter_reagent_families()))}")
    print(f"  Aliases: {len(registry.aliases)}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate unified taxonomy data files.")
    parser.add_argument(
        "--root",
        type=Path,
        default=None,
        help="Directory containing taxonomy JSON files (defaults to chemtools/taxonomy/data).",
    )
    parser.add_argument(
        "--reset-cache",
        action="store_true",
        help="Clear cached registry instance before loading.",
    )
    args = parser.parse_args()

    if args.reset_cache:
        reset_registry(args.root)

    validate(args.root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
