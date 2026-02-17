from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.taxonomy import loader
from chemtools.taxonomy.motif_scope_index import build_and_save_motif_scope_index


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prebuild motif scope index for scope-aware recommendation matching."
    )
    parser.add_argument(
        "--output",
        type=str,
        default=str(loader.MOTIF_SCOPE_INDEX_FILE),
        help=f"Output JSON path. Default: {loader.MOTIF_SCOPE_INDEX_FILE}",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print JSON summary.",
    )
    args = parser.parse_args()

    output_path = Path(args.output)
    payload = build_and_save_motif_scope_index(output_path)
    scope_map = payload.get("scope_map", {}) if isinstance(payload, dict) else {}
    scaffold_parent_map = payload.get("scaffold_parent_map", {}) if isinstance(payload, dict) else {}
    summary = {
        "output_path": str(output_path),
        "version": payload.get("version"),
        "rdkit_available": bool(payload.get("rdkit_available")),
        "scope_parents": len(scope_map),
        "scope_edges": sum(len(v or []) for v in (scope_map or {}).values()),
        "scaffold_children": len(scaffold_parent_map),
    }

    if args.json:
        print(json.dumps(summary, indent=2))
        return 0

    print("Motif scope index build completed")
    print(f"Output: {summary['output_path']}")
    print(f"Version: {summary['version']}")
    print(f"RDKit available: {summary['rdkit_available']}")
    print(f"Scope parents: {summary['scope_parents']}")
    print(f"Scope edges: {summary['scope_edges']}")
    print(f"Scaffold children: {summary['scaffold_children']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

