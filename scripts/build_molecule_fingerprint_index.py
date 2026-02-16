#!/usr/bin/env python3
"""
Build a bulk molecule fingerprint NPZ index from .smi/.sdf(.gz) files.

Examples:
  python scripts/build_molecule_fingerprint_index.py \
    --input data/reagent_db/reagents.smi \
    --output results/reagents_morgan_fp.npz

  python scripts/build_molecule_fingerprint_index.py \
    --input data/large.sdf.gz \
    --output results/large_rdkit_fp.npz \
    --fp-type rdkit --fp-size 4096 --num-threads 8
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from chemtools.util.rdkit_bulk_fingerprints import (  # noqa: E402
    fingerprints_from_molecule_file,
    save_fingerprint_batch,
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Bulk fingerprint molecules from file into compressed NPZ.",
    )
    parser.add_argument("--input", required=True, help="Input .smi/.sdf(.gz) file")
    parser.add_argument("--output", required=True, help="Output .npz path")
    parser.add_argument(
        "--fp-type",
        choices=("morgan", "rdkit", "atom_pair", "topological_torsion"),
        default="morgan",
        help="Fingerprint algorithm",
    )
    parser.add_argument("--fp-size", type=int, default=2048, help="Fingerprint size")
    parser.add_argument("--radius", type=int, default=2, help="Morgan radius")
    parser.add_argument(
        "--include-chirality",
        action="store_true",
        help="Include chirality when supported",
    )
    parser.add_argument("--num-threads", type=int, default=0, help="0 means all cores")
    parser.add_argument(
        "--delimiter",
        default="\t",
        help="Delimiter for SMILES-like files (default: tab)",
    )
    parser.add_argument("--smiles-column", type=int, default=0, help="SMILES column index")
    parser.add_argument("--name-column", type=int, default=1, help="Name column index")
    parser.add_argument(
        "--no-title-line",
        action="store_true",
        help="Input has no header row for SMILES-like files",
    )
    parser.add_argument(
        "--no-sanitize",
        action="store_true",
        help="Disable RDKit sanitize for faster ingest on trusted inputs",
    )
    parser.add_argument(
        "--keep-invalid",
        action="store_true",
        help="Keep invalid rows as all-zero vectors instead of dropping",
    )
    parser.add_argument(
        "--force-fallback",
        action="store_true",
        help="Disable rdMolProcessing fast path and use Python supplier loop",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()

    batch = fingerprints_from_molecule_file(
        args.input,
        fp_type=args.fp_type,
        fp_size=int(args.fp_size),
        radius=int(args.radius),
        include_chirality=bool(args.include_chirality),
        num_threads=int(args.num_threads),
        delimiter=str(args.delimiter),
        smiles_column=int(args.smiles_column),
        name_column=int(args.name_column),
        title_line=not bool(args.no_title_line),
        sanitize=not bool(args.no_sanitize),
        skip_invalid=not bool(args.keep_invalid),
        prefer_bulk=not bool(args.force_fallback),
    )
    save_fingerprint_batch(batch, args.output)

    summary = batch.to_summary()
    print(f"Input: {Path(args.input).resolve()}")
    print(f"Output: {Path(args.output).resolve()}")
    print(f"Engine: {summary['engine']}")
    print(f"FP type: {summary['fp_type']}")
    print(f"FP size: {summary['fp_size']}")
    print(f"Records: total={summary['total_records']}, valid={summary['valid_records']}, invalid={summary['invalid_records']}")
    print(f"Matrix shape: {summary['shape']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
