"""Temporary guarded runner for migrating the local Full recommendation library."""

from __future__ import annotations

import json
from pathlib import Path

from condition_recommender.conversion.artifacts import (
    combine_saved_recommendation_batches,
    resume_saved_conversion_batch,
)


def main() -> None:
    """Resume the retained legacy batch and rebuild all Full artifacts."""
    full = Path("datasets/literature/full").resolve()
    manifest = full / "shard_manifest.json"
    last = {"phase": "", "row_bucket": -1, "shards": -1}

    def progress(value: object) -> None:
        phase = str(getattr(value, "phase"))
        row_count = int(getattr(value, "row_count"))
        shard_count = int(getattr(value, "shard_count"))
        row_bucket = row_count // 25_000
        important = (
            phase != last["phase"]
            or row_bucket != last["row_bucket"]
            or (
                shard_count
                and shard_count % 25 == 0
                and shard_count != last["shards"]
            )
        )
        if important:
            print(
                f"[{phase}] sources={getattr(value, 'source_file_count')} "
                f"shards={shard_count} rows={row_count}: "
                f"{getattr(value, 'message')}",
                flush=True,
            )
            last.update(
                phase=phase,
                row_bucket=row_bucket,
                shards=shard_count,
            )

    print(f"Resuming {manifest}", flush=True)
    resumed = resume_saved_conversion_batch(
        manifest,
        workers=6,
        progress_callback=progress,
    )
    print(
        "RESUME_COMPLETE "
        + json.dumps(
            {
                "record_count": resumed.get("record_count"),
                "reused_shard_count": resumed.get("reused_shard_count"),
                "shard_count": resumed.get("shard_count"),
                "output_dir": resumed.get("output_dir"),
            },
            sort_keys=True,
        ),
        flush=True,
    )
    last.update(phase="", row_bucket=-1, shards=-1)
    print(f"Combining saved Full batches in {full}", flush=True)
    combined = combine_saved_recommendation_batches(
        full,
        progress_callback=progress,
    )
    print(
        "COMBINE_COMPLETE "
        + json.dumps(
            {
                "record_count": combined.get("record_count"),
                "input_record_count": combined.get("input_record_count"),
                "duplicate_record_count": combined.get(
                    "duplicate_record_count"
                ),
                "trusted_precedent_count": combined.get(
                    "trusted_precedent_count"
                ),
                "review_core_precedent_count": combined.get(
                    "review_core_precedent_count"
                ),
                "batch_count": combined.get("batch_count"),
                "output_dir": combined.get("output_dir"),
            },
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
