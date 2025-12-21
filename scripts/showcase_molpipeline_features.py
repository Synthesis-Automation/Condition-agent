#!/usr/bin/env python3
"""Showcase MolPipeline-backed feature generation inside ChemTools.

Usage:
    python scripts/showcase_molpipeline_features.py \
        --electrophile Brc1ccccc1 \
        --nucleophile Nc1ccccc1

The script prints the standard C-N coupling feature dictionary and, when
MolPipeline is installed, the richer fingerprint/descriptor payload that is now
available via ``include_molpipeline=True``.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

# Ensure project root is on sys.path when run directly
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from chemtools.featurizers.reaction_pair import featurize_pair

try:
    from chemtools.integrations import molpipeline as mp
except Exception:  # pragma: no cover - MolPipeline optional
    mp = None  # type: ignore


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Demonstrate ChemTools to MolPipeline integration.",
    )
    parser.add_argument(
        "--electrophile",
        default="Brc1ccccc1",
        help="Electrophile SMILES (default: %(default)s)",
    )
    parser.add_argument(
        "--nucleophile",
        default="Nc1ccccc1",
        help="Nucleophile SMILES (default: %(default)s)",
    )
    parser.add_argument(
        "--indent",
        type=int,
        default=2,
        help="Pretty-print JSON with the given indentation (default: %(default)s)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    include_vectors = bool(mp and mp.is_available())
    if not include_vectors:
        print(
            "MolPipeline not detected - standard ChemTools features will be returned "
            "without the optional fingerprint payload.",
            file=sys.stderr,
        )

    features = featurize_pair(
        args.electrophile,
        args.nucleophile,
        include_molpipeline=include_vectors,
    )

    print(json.dumps(features, indent=args.indent, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
