"""
Command-line interface for the rule builder wizard.

Usage:
    python -m chem_assistant.rule_builder_cli --family Suzuki_Miyaura --output data/rule_db_v2/Suzuki_new.json
    python -m chem_assistant.rule_builder_cli --load data/rule_db_v2/Suzuki_db.json
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List

from chemtools.rule import RuleBuilder

from .rule_builder_session import RuleBuilderSession


def run_cli(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Interactive wizard for creating or editing rule databases."
    )
    parser.add_argument(
        "--family",
        help="Reaction family name (required when creating a new database).",
    )
    parser.add_argument(
        "--load",
        help="Path to an existing rule DB JSON file to edit.",
    )
    parser.add_argument(
        "--output",
        help="Destination path for saving the rule DB.",
    )
    args = parser.parse_args(argv)

    if args.load:
        source_path = Path(args.load)
        builder = RuleBuilder.from_file(source_path)
        target_path = Path(args.output) if args.output else source_path
        print(f"Loaded existing database from {source_path}")
    else:
        if not args.family:
            parser.error("--family is required when not using --load")
        if not args.output:
            parser.error("--output is required when creating a new database")
        builder = RuleBuilder.new(args.family)
        target_path = Path(args.output)
        print(f"Creating new database for family '{args.family}'")

    session = RuleBuilderSession(builder)
    session.run_wizard()

    issues = builder.validate(strict=False)
    if issues:
        print("\nValidation results:")
        for issue in issues:
            print(f"  [{issue.severity.upper()}] {issue.field}: {issue.message}")
        blocking = [i for i in issues if i.severity == "error"]
        if blocking:
            print("\nErrors detected; please re-run the wizard to fix them.")
            return 1

    diff = builder.diff()
    if diff.strip():
        print("\nChanges preview:\n")
        print(diff)
    else:
        print("\nNo changes detected relative to the source snapshot.")

    save_choice = input(f"\nSave to {target_path}? (y/n) [y]: ").strip().lower()
    if save_choice in {"", "y", "yes"}:
        builder.save(target_path)
        print(f"Saved rule database to {target_path}")
    else:
        print("Aborted save; no files were written.")
    return 0


def main() -> None:
    sys.exit(run_cli())


if __name__ == "__main__":
    main()
