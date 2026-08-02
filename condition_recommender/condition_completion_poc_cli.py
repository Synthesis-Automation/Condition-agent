"""CLI for exact-core and learned-group condition context completion."""

from __future__ import annotations

import argparse
import json

from .condition_completion_poc import complete_condition_core


def main() -> None:
    """Complete protocol context for one role-aware decisive core."""

    parser = argparse.ArgumentParser(
        description="Find observed or group-prior protocol context for a condition core"
    )
    parser.add_argument("artifact_dir")
    parser.add_argument("condition_core")
    parser.add_argument("--top-k", type=int, default=5)
    args = parser.parse_args()
    result = complete_condition_core(
        args.condition_core,
        args.artifact_dir,
        top_k=args.top_k,
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
