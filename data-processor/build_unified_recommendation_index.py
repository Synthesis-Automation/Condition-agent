#!/usr/bin/env python3
"""
Build the unified recommendation index (dataset + protocol + HTE).

Usage:
  python data-processor/build_unified_recommendation_index.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from chemtools.recommend.index_builder import build_unified_recommendation_index


def main() -> int:
    parser = argparse.ArgumentParser(description="Build unified recommendation index")
    parser.add_argument("--reaction-dataset", default="data/reaction_dataset")
    parser.add_argument("--protocols", default="data/protocol_db_v2")
    parser.add_argument("--hte", default="data/HTE_db")
    parser.add_argument("--output", default="build/unified_recommendation_index")
    parser.add_argument("--no-hte", action="store_true")
    parser.add_argument("--max-records", type=int, default=None)
    parser.add_argument("--skip-drfp", action="store_true")
    args = parser.parse_args()

    stats = build_unified_recommendation_index(
        reaction_dataset_dir=args.reaction_dataset,
        protocol_dir=args.protocols,
        hte_dir=args.hte,
        output_dir=args.output,
        include_hte=not args.no_hte,
        max_records=args.max_records,
        skip_drfp=args.skip_drfp,
    )

    print("Build complete.")
    print(f"Total entries: {stats.total_entries}")
    print(f"Sources: {stats.by_source}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
