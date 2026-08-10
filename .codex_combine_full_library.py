"""Temporary runner for rebuilding the migrated Full recommendation library."""

from __future__ import annotations

import json
from pathlib import Path

from condition_recommender.conversion.artifacts import (
    combine_saved_recommendation_batches,
)


def main() -> None:
    """Combine completed Full batches and print throttled progress."""
    full = Path("datasets/literature/full").resolve()
    last = {"phase": "", "row_bucket": -1}

    def progress(value: object) -> None:
        phase = str(getattr(value, "phase"))
        row_count = int(getattr(value, "row_count"))
        row_bucket = row_count // 25_000
        if phase != last["phase"] or row_bucket != last["row_bucket"]:
            print(
                f"[{phase}] rows={row_count}: {getattr(value, 'message')}",
                flush=True,
            )
            last.update(phase=phase, row_bucket=row_bucket)

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
