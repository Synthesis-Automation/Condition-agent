"""CLI for composing final generic-recommender release gates."""

from __future__ import annotations

import argparse
import json

from .release_validation import build_release_readiness_report


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compose machine and chemist-review release gates"
    )
    parser.add_argument("conversion_report_path")
    parser.add_argument("index_integrity_path")
    parser.add_argument("grouped_evaluation_path")
    parser.add_argument("scaffold_evaluation_path")
    parser.add_argument("untouched_evaluation_path")
    parser.add_argument("calibration_report_path")
    parser.add_argument("output_path")
    parser.add_argument("--chemist-review-summary-path")
    args = parser.parse_args()
    report = build_release_readiness_report(
        args.conversion_report_path,
        args.index_integrity_path,
        args.grouped_evaluation_path,
        args.scaffold_evaluation_path,
        args.untouched_evaluation_path,
        args.calibration_report_path,
        args.output_path,
        chemist_review_summary_path=args.chemist_review_summary_path,
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if not report["machine_release_ready"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
