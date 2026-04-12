from __future__ import annotations

import argparse
import tempfile
import time
from pathlib import Path

from app.A_convert_to_hte_format import process_reaction_dataset


def _parse_workers(text: str) -> list[int]:
    values: list[int] = []
    for part in str(text or "").split(","):
        token = part.strip()
        if not token:
            continue
        values.append(int(token))
    return values or [1]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Benchmark HTE CSV conversion with different worker counts."
    )
    parser.add_argument("--input", "-i", required=True, help="Input CSV file.")
    parser.add_argument(
        "--workers",
        default="1,0",
        help="Comma-separated worker counts. Use 0 for auto. Default: 1,0",
    )
    parser.add_argument(
        "--repeat",
        type=int,
        default=1,
        help="Runs per worker count. Default: 1",
    )
    parser.add_argument(
        "--keep-no-catalyst",
        action="store_true",
        help="Keep reactions without catalyst during benchmark.",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    worker_values = _parse_workers(args.workers)
    print(f"Benchmarking {input_path}")
    print(f"Worker settings: {worker_values}")
    print(f"Repeats: {args.repeat}")

    for workers in worker_values:
        elapsed_runs: list[float] = []
        for run_idx in range(1, args.repeat + 1):
            with tempfile.TemporaryDirectory(prefix="hte_bench_") as tmpdir:
                output_path = Path(tmpdir) / f"{input_path.stem}_canonical.csv"
                started = time.perf_counter()
                process_reaction_dataset(
                    str(input_path),
                    str(output_path),
                    drop_no_catalyst=not args.keep_no_catalyst,
                    num_workers=workers,
                )
                elapsed = time.perf_counter() - started
                elapsed_runs.append(elapsed)
                print(
                    f"workers={workers} run={run_idx}/{args.repeat} elapsed={elapsed:.2f}s"
                )

        avg = sum(elapsed_runs) / len(elapsed_runs)
        best = min(elapsed_runs)
        worst = max(elapsed_runs)
        print(
            f"workers={workers} summary avg={avg:.2f}s best={best:.2f}s worst={worst:.2f}s"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
