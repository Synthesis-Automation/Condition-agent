"""Validate and adjudicate a completed blind chemist-review packet."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Sequence

from .chemist_review import CHEMIST_REVIEW_PACKET_SCHEMA_VERSION

CHEMIST_REVIEW_SUMMARY_SCHEMA_VERSION = "1.0"
ALLOWED_CHEMIST_DECISIONS = frozenset(
    {
        "compatible",
        "compatible_with_caution",
        "incompatible",
        "uncertain",
    }
)
ACCEPTABLE_CHEMIST_DECISIONS = frozenset(
    {"compatible", "compatible_with_caution"}
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_jsonl(path: Path) -> list[Dict[str, Any]]:
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                value = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(
                    f"Invalid review JSONL at {path}:{line_number}: {exc.msg}"
                ) from exc
            if not isinstance(value, dict):
                raise ValueError(
                    f"Review JSONL row is not an object at {path}:{line_number}"
                )
            rows.append(value)
    return rows


def _read_completed_form(path: Path) -> list[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = [
            {str(key): str(value or "").strip() for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]
    if not rows:
        raise ValueError("Chemist review form is empty")
    missing = []
    invalid = []
    seen = set()
    for row_number, row in enumerate(rows, start=2):
        pair = (row.get("case_id", ""), row.get("candidate_id", ""))
        if not all(pair):
            missing.append(f"row {row_number}:missing_case_or_candidate_id")
        if pair in seen:
            invalid.append(f"row {row_number}:duplicate_case_candidate")
        seen.add(pair)
        decision = row.get("chemist_decision", "").casefold()
        if not decision:
            missing.append(f"row {row_number}:missing_decision")
        elif decision not in ALLOWED_CHEMIST_DECISIONS:
            invalid.append(f"row {row_number}:invalid_decision:{decision}")
        row["chemist_decision"] = decision
    if missing:
        raise ValueError(
            "Review must be complete before unblinding: " + "; ".join(missing[:10])
        )
    if invalid:
        raise ValueError("Invalid review form: " + "; ".join(invalid[:10]))
    return rows


def _rate(numerator: int, denominator: int) -> float | None:
    return round(numerator / denominator, 6) if denominator else None


def _normalize_defects(values: Iterable[str]) -> tuple[str, ...]:
    return tuple(
        sorted(
            {
                str(value).strip()
                for value in values
                if str(value).strip()
            },
            key=str.casefold,
        )
    )


def adjudicate_chemist_review(
    review_dir: str | Path,
    output_path: str | Path,
    *,
    reviewer_name: str,
    independent_reviewer: bool,
    reviewer_signoff: bool,
    unresolved_systematic_defects: Sequence[str] = (),
) -> Dict[str, Any]:
    """Unblind a complete form and create a provenance-bound review summary.

    The answer key is intentionally read only after every candidate has a valid
    decision. ``reviewer_signoff`` is an explicit human assertion that the
    reported findings were considered and no undeclared systematic defect
    remains.
    """
    reviewer = reviewer_name.strip()
    if not reviewer:
        raise ValueError("reviewer_name is required")
    source = Path(review_dir)
    form_path = source / "review_form.csv"
    answer_path = source / "answer_key.jsonl"
    packet_path = source / "review_packet.jsonl"
    html_packet_path = source / "review_packet.html"
    report_path = source / "review_report.json"
    for path in (
        form_path,
        answer_path,
        packet_path,
        html_packet_path,
        report_path,
    ):
        if not path.is_file():
            raise ValueError(f"Missing chemist-review artifact: {path}")

    # Preserve blinding by refusing to open the key until the form is complete.
    form_rows = _read_completed_form(form_path)
    answer_rows = _read_jsonl(answer_path)
    answer_by_pair = {}
    for row in answer_rows:
        pair = (str(row.get("case_id") or ""), str(row.get("candidate_id") or ""))
        if not all(pair) or pair in answer_by_pair:
            raise ValueError("Answer key has missing or duplicate candidate identity")
        answer_by_pair[pair] = row
    form_by_pair = {
        (row["case_id"], row["candidate_id"]): row for row in form_rows
    }
    if set(form_by_pair) != set(answer_by_pair):
        missing_form = sorted(set(answer_by_pair) - set(form_by_pair))
        missing_key = sorted(set(form_by_pair) - set(answer_by_pair))
        raise ValueError(
            "Review form and answer key candidate sets differ: "
            f"missing_form={len(missing_form)}, missing_key={len(missing_key)}"
        )

    packet_report = json.loads(report_path.read_text(encoding="utf-8"))
    if (
        packet_report.get("schema_version")
        != CHEMIST_REVIEW_PACKET_SCHEMA_VERSION
        or packet_report.get("artifact_type")
        != "generic_chemist_review_packet"
    ):
        raise ValueError("Unsupported chemist-review packet report")
    if int(packet_report.get("candidate_count") or -1) != len(form_rows):
        raise ValueError("Review form count does not match packet report")

    decision_counts = Counter()
    source_decision_counts: Dict[str, Counter[str]] = defaultdict(Counter)
    cases_with_acceptable_recommendation = set()
    recommendation_cases = set()
    accepted_controls = []
    uncertain_pairs = []
    for pair, form_row in form_by_pair.items():
        answer = answer_by_pair[pair]
        decision = form_row["chemist_decision"]
        source_kind = str(answer.get("source_kind") or "")
        if source_kind not in {"recommendation", "negative_control"}:
            raise ValueError(f"Unknown answer-key source kind: {source_kind}")
        decision_counts[decision] += 1
        source_decision_counts[source_kind][decision] += 1
        if source_kind == "recommendation":
            recommendation_cases.add(pair[0])
            if decision in ACCEPTABLE_CHEMIST_DECISIONS:
                cases_with_acceptable_recommendation.add(pair[0])
        elif decision in ACCEPTABLE_CHEMIST_DECISIONS:
            accepted_controls.append(
                {
                    "case_id": pair[0],
                    "candidate_id": pair[1],
                    "decision": decision,
                    "comments": form_row.get("chemist_comments", ""),
                }
            )
        if decision == "uncertain":
            uncertain_pairs.append(
                {"case_id": pair[0], "candidate_id": pair[1]}
            )

    recommendation_count = sum(source_decision_counts["recommendation"].values())
    accepted_recommendation_count = sum(
        source_decision_counts["recommendation"][decision]
        for decision in ACCEPTABLE_CHEMIST_DECISIONS
    )
    control_count = sum(source_decision_counts["negative_control"].values())
    rejected_control_count = source_decision_counts["negative_control"][
        "incompatible"
    ]
    defects = _normalize_defects(unresolved_systematic_defects)
    review_complete = bool(
        independent_reviewer
        and reviewer_signoff
        and not defects
    )
    summary: Dict[str, Any] = {
        "schema_version": CHEMIST_REVIEW_SUMMARY_SCHEMA_VERSION,
        "artifact_type": "generic_chemist_review_summary",
        "completed": review_complete,
        "reviewer": {
            "name": reviewer,
            "independent": independent_reviewer,
            "signed_off": reviewer_signoff,
        },
        "unresolved_systematic_defect_count": len(defects),
        "unresolved_systematic_defects": list(defects),
        "metrics": {
            "case_count": int(packet_report.get("case_count") or 0),
            "candidate_count": len(form_rows),
            "decision_counts": dict(sorted(decision_counts.items())),
            "recommendation_candidate_count": recommendation_count,
            "acceptable_recommendation_count": accepted_recommendation_count,
            "acceptable_recommendation_rate": _rate(
                accepted_recommendation_count,
                recommendation_count,
            ),
            "cases_with_recommendation_count": len(recommendation_cases),
            "cases_with_acceptable_recommendation_count": len(
                cases_with_acceptable_recommendation
            ),
            "case_acceptable_recommendation_rate": _rate(
                len(cases_with_acceptable_recommendation),
                len(recommendation_cases),
            ),
            "negative_control_count": control_count,
            "rejected_negative_control_count": rejected_control_count,
            "negative_control_rejection_rate": _rate(
                rejected_control_count,
                control_count,
            ),
            "accepted_negative_control_count": len(accepted_controls),
            "uncertain_candidate_count": len(uncertain_pairs),
        },
        "review_findings": {
            "accepted_negative_controls": accepted_controls,
            "uncertain_candidates": uncertain_pairs,
        },
        "provenance": {
            "packet_report_sha256": _sha256(report_path),
            "review_packet_sha256": _sha256(packet_path),
            "html_packet_sha256": _sha256(html_packet_path),
            "review_form_sha256": _sha256(form_path),
            "answer_key_sha256": _sha256(answer_path),
        },
    }
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def valid_completed_chemist_review(summary: Mapping[str, Any]) -> bool:
    """Return whether a review summary satisfies the independent human gate."""
    reviewer = summary.get("reviewer")
    provenance = summary.get("provenance")
    return bool(
        summary.get("schema_version") == CHEMIST_REVIEW_SUMMARY_SCHEMA_VERSION
        and summary.get("artifact_type") == "generic_chemist_review_summary"
        and summary.get("completed")
        and isinstance(reviewer, Mapping)
        and str(reviewer.get("name") or "").strip()
        and reviewer.get("independent") is True
        and reviewer.get("signed_off") is True
        and int(summary.get("unresolved_systematic_defect_count") or 0) == 0
        and isinstance(provenance, Mapping)
        and all(
            str(provenance.get(key) or "").strip()
            for key in (
                "packet_report_sha256",
                "review_packet_sha256",
                "html_packet_sha256",
                "review_form_sha256",
                "answer_key_sha256",
            )
        )
    )


__all__ = [
    "ACCEPTABLE_CHEMIST_DECISIONS",
    "ALLOWED_CHEMIST_DECISIONS",
    "CHEMIST_REVIEW_SUMMARY_SCHEMA_VERSION",
    "adjudicate_chemist_review",
    "valid_completed_chemist_review",
]
