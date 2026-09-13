"""Build and audit experimental shared-core artifacts from canonical indices."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .generic_api import GenericConditionRecommender
from .generic_indexing import load_generic_index
from .shared_core_index import build_shared_core_index


def main(argv: list[str] | None = None) -> int:
    """Build derived keys or emit an auditable dual-channel query result."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    build = commands.add_parser(
        "build", help="Backfill stored graph projections atomically"
    )
    build.add_argument("index", type=Path)
    build.add_argument("output", type=Path)
    audit = commands.add_parser(
        "query", help="Run experimental shared-core recommendation"
    )
    audit.add_argument("index", type=Path)
    audit.add_argument("projections", type=Path)
    audit.add_argument("reaction")
    audit.add_argument(
        "--scope", choices=("same_handle", "automatic", "broad"), default="automatic"
    )
    audit.add_argument("--retro-precedent", action="append", default=[])
    audit.add_argument("--top-k", type=int, default=5)
    args = parser.parse_args(argv)
    if args.command == "build":
        if args.index.resolve() == args.output.resolve():
            parser.error("output must differ from the canonical source index")
        result = build_shared_core_index(load_generic_index(args.index), args.output)
    else:
        recommender = GenericConditionRecommender.from_path(
            args.index, shared_core_path=args.projections
        )
        result = recommender.recommend(
            args.reaction,
            search_scope=args.scope,
            preferred_reaction_ids=tuple(args.retro_precedent),
            top_k=args.top_k,
        ).to_dict()
    print(json.dumps(result, indent=2, ensure_ascii=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
