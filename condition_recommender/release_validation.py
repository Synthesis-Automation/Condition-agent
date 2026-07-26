"""Compose machine and human release gates for the generic recommender."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional

from .chemist_review_adjudication import valid_completed_chemist_review


def _read(path: str | Path) -> Dict[str, Any]:
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Release artifact is not a JSON object: {path}")
    return value


def build_release_readiness_report(
    conversion_report_path: str | Path,
    index_integrity_path: str | Path,
    grouped_evaluation_path: str | Path,
    scaffold_evaluation_path: str | Path,
    untouched_evaluation_path: str | Path,
    calibration_report_path: str | Path,
    output_path: str | Path,
    *,
    chemist_review_summary_path: Optional[str | Path] = None,
) -> Dict[str, Any]:
    """Validate required artifacts without converting review into model evidence."""
    conversion = _read(conversion_report_path)
    index = _read(index_integrity_path)
    grouped = _read(grouped_evaluation_path)
    scaffold = _read(scaffold_evaluation_path)
    untouched = _read(untouched_evaluation_path)
    calibration = _read(calibration_report_path)
    review = (
        _read(chemist_review_summary_path)
        if chemist_review_summary_path is not None
        else None
    )
    gates = {
        "conversion_complete": (
            int(conversion.get("failed_shard_count") or 0) == 0
            and int(conversion.get("input_row_count") or 0)
            == int(conversion.get("output_row_count") or -1)
            and bool((conversion.get("integrity") or {}).get("valid"))
        ),
        "index_integrity": bool(index.get("valid")),
        "reference_reaction_leakage_zero": all(
            int(report["split"]["reference_overlap_count"]) == 0
            and int(report["split"]["canonical_reaction_overlap_count"]) == 0
            for report in (grouped, scaffold, untouched)
        ),
        "strict_scaffold_leakage_zero": (
            int(scaffold["split"]["scaffold_overlap_count"]) == 0
        ),
        "hard_incompatible_recommendations_zero": all(
            int(
                report["metrics"][
                    "hard_incompatible_recommendation_count"
                ]
            )
            == 0
            for report in (grouped, scaffold, untouched)
        ),
        "sample_calibration_promoted": bool(calibration.get("promoted")),
        "untouched_test_executed": (
            int(untouched["metrics"]["query_count"]) > 0
        ),
    }
    machine_ready = all(gates.values())
    human_review_complete = bool(
        review and valid_completed_chemist_review(review)
    )
    report: Dict[str, Any] = {
        "schema_version": "1.0",
        "artifact_type": "generic_recommender_release_readiness",
        "machine_release_ready": machine_ready,
        "production_release_ready": machine_ready and human_review_complete,
        "gates": gates,
        "human_review": {
            "summary_supplied": review is not None,
            "complete_without_systematic_defect": human_review_complete,
        },
        "supported_scope": {
            "indexed_transformation_class_counts": index.get(
                "transformation_class_counts", {}
            ),
            "indexed_named_family_counts": index.get(
                "named_family_counts", {}
            ),
            "indexed_reaction_scope_counts": index.get(
                "reaction_scope_counts", {}
            ),
        },
        "abstaining_scope": {
            "admission_reason_counts": conversion.get(
                "admission_reason_counts", {}
            ),
            "ineligible_record_count": (
                conversion.get("index_eligibility_counts", {}).get(
                    "ineligible", 0
                )
            ),
            "review_only_record_count": (
                conversion.get("index_eligibility_counts", {}).get(
                    "review_only", 0
                )
            ),
        },
        "untouched_test_metrics": untouched["metrics"],
        "definition_versions": untouched.get("definition_versions", {}),
        "artifact_paths": {
            "conversion_report": str(Path(conversion_report_path)),
            "index_integrity": str(Path(index_integrity_path)),
            "grouped_evaluation": str(Path(grouped_evaluation_path)),
            "scaffold_evaluation": str(Path(scaffold_evaluation_path)),
            "untouched_evaluation": str(Path(untouched_evaluation_path)),
            "calibration_report": str(Path(calibration_report_path)),
            "chemist_review_summary": (
                str(Path(chemist_review_summary_path))
                if chemist_review_summary_path is not None
                else None
            ),
        },
    }
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return report


__all__ = ["build_release_readiness_report"]
