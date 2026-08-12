"""Command-line builder for reusable canonical molecule indexes."""

from __future__ import annotations

import argparse
import json
from typing import Optional, Sequence

from .molecule_index import build_canonical_molecule_index


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a canonical SQLite molecule-identity index"
    )
    parser.add_argument("source_csv")
    parser.add_argument("output_sqlite")
    parser.add_argument("--smiles-column", default="compound_smiles")
    parser.add_argument(
        "--provenance-column",
        action="append",
        default=[],
        help="source column retained with each match; repeat as needed",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Build the requested index and print its deterministic report."""

    arguments = _parser().parse_args(argv)
    report = build_canonical_molecule_index(
        arguments.source_csv,
        arguments.output_sqlite,
        smiles_column=arguments.smiles_column,
        provenance_columns=arguments.provenance_column,
    )
    print(json.dumps(report.to_dict(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":  # pragma: no cover - module entry point
    raise SystemExit(main())

