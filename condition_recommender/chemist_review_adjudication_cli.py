"""CLI for adjudicating a completed blind chemist-review form."""

from __future__ import annotations

import argparse
import json

from .chemist_review_adjudication import adjudicate_chemist_review


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Unblind and summarize a completed chemist-review packet"
    )
    parser.add_argument("review_dir")
    parser.add_argument("output_path")
    parser.add_argument("--reviewer", required=True)
    parser.add_argument("--independent-reviewer", action="store_true")
    parser.add_argument("--sign-off", action="store_true")
    parser.add_argument(
        "--unresolved-defect",
        action="append",
        default=[],
        help="Describe an unresolved systematic defect; repeat as needed",
    )
    args = parser.parse_args()
    summary = adjudicate_chemist_review(
        args.review_dir,
        args.output_path,
        reviewer_name=args.reviewer,
        independent_reviewer=args.independent_reviewer,
        reviewer_signoff=args.sign_off,
        unresolved_systematic_defects=args.unresolved_defect,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
