from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.recommend import warm_hte_cache


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prebuild HTE recommender cache/index files for faster startup."
    )
    parser.add_argument(
        "--db-path",
        type=str,
        default="data/HTE_db",
        help="HTE database path (folder or CSV). Default: data/HTE_db",
    )
    parser.add_argument(
        "--source-group",
        type=str,
        default="all",
        help="Optional source group target: all, experiments, literature, protocols, rules",
    )
    parser.add_argument(
        "--clear-memory-cache",
        action="store_true",
        help="Clear in-process LRU cache before warming.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print JSON summary.",
    )
    args = parser.parse_args()

    summary = warm_hte_cache(
        args.db_path,
        source_group=args.source_group,
        clear_memory_cache=args.clear_memory_cache,
    )

    if args.json:
        print(json.dumps(summary, indent=2))
        return 0

    print("HTE cache warm completed")
    print(f"Base path: {summary.get('base_path')}")
    print(f"Source group: {summary.get('source_group')}")
    for item in summary.get("targets", []):
        print("-" * 72)
        print(f"Target: {item.get('target')}")
        print(f"Cache dir: {item.get('cache_dir')}")
        print(f"Rows: {item.get('num_rows')}")
        print(f"Reactant index keys: {item.get('reactant_index_keys')}")
        print(f"Transformation keys: {item.get('transformation_index_keys')}")
        print(f"Elapsed: {item.get('elapsed_s')} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

