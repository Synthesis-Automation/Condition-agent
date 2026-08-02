"""CLI for the experimental weak-label condition grouping POC."""

from __future__ import annotations

import argparse
import json

from .condition_grouping_poc import run_condition_grouping_poc


def main() -> None:
    """Run the role-aware grouping POC and print its artifact report."""

    parser = argparse.ArgumentParser(
        description=(
            "Group structured weak-label condition recipes without using "
            "reaction labels or outcomes as clustering features"
        )
    )
    parser.add_argument("dataset_path")
    parser.add_argument("output_dir")
    parser.add_argument("--clusters", type=int, default=256)
    parser.add_argument("--latent-dimensions", type=int, default=48)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--silhouette-sample-size", type=int, default=3000)
    args = parser.parse_args()
    result = run_condition_grouping_poc(
        args.dataset_path,
        args.output_dir,
        cluster_count=args.clusters,
        latent_dimensions=args.latent_dimensions,
        seed=args.seed,
        silhouette_sample_size=args.silhouette_sample_size,
    )
    print(json.dumps(result.report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
