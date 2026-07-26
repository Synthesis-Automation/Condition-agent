"""CLI for restartable full-corpus generic conversion."""

from __future__ import annotations

import argparse
import json

from .conversion.sharded import convert_datasets_sharded


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert reaction datasets into restartable canonical shards"
    )
    parser.add_argument("dataset_path")
    parser.add_argument("output_dir")
    parser.add_argument("--shard-size", type=int, default=1_000)
    parser.add_argument("--max-shards", type=int)
    parser.add_argument("--mode", default="full")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--checkpoint-interval", type=int, default=10)
    args = parser.parse_args()
    report = convert_datasets_sharded(
        args.dataset_path,
        args.output_dir,
        shard_size=args.shard_size,
        max_shards=args.max_shards,
        mode=args.mode,
        workers=args.workers,
        checkpoint_interval=args.checkpoint_interval,
        progress=True,
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
