"""Deterministic per-reaction audit of generic source round trips."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Literal, Sequence

from .core_admission import (
    CoreAdmissionPolicyName,
    load_generic_core_admission_policy,
)
from .generic_compiler import compile_generic_templates
from .sources import iter_library_rows, source_shard_files


ROUND_TRIP_AUDIT_SCHEMA_VERSION = "1.0"
ROUND_TRIP_AUDIT_ALGORITHM_VERSION = "generic_round_trip_audit.v1"
RoundTripLevel = Literal["L0", "L1", "L2"]


def _sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def audit_generic_round_trips(
    source: str | Path,
    destination: str | Path,
    *,
    reaction_ids: Sequence[str],
    core_admission_policy: CoreAdmissionPolicyName = "pass_only",
    levels: Sequence[RoundTripLevel] = ("L0", "L1", "L2"),
) -> dict[str, Any]:
    """Compile selected rows and write transparent per-level round-trip results."""

    requested = tuple(dict.fromkeys(str(item) for item in reaction_ids if item))
    if not requested:
        raise ValueError("round-trip audit requires reaction IDs")
    if len(requested) != len(reaction_ids):
        raise ValueError("round-trip audit reaction IDs must be unique and nonempty")
    normalized_levels = tuple(dict.fromkeys(levels))
    if not normalized_levels or any(
        item not in {"L0", "L1", "L2"} for item in normalized_levels
    ):
        raise ValueError("round-trip audit levels must be L0, L1, and/or L2")
    policy = load_generic_core_admission_policy(core_admission_policy)
    requested_set = set(requested)
    found: dict[str, dict[str, Any]] = {}
    for row in iter_library_rows(source):
        reaction_id = str(row.get("reaction_id") or "")
        if reaction_id not in requested_set:
            continue
        if reaction_id in found:
            raise ValueError(f"duplicate audited reaction ID: {reaction_id}")
        result = compile_generic_templates(
            row,
            levels=normalized_levels,
            admission_mode="data_driven",
            core_admission_policy=core_admission_policy,
        )
        diagnostics = dict(result.diagnostics or {})
        templates = tuple(
            {
                "abstraction_level": template.abstraction_level,
                "operator_id": template.operator_id,
                "reaction_smarts": template.reaction_smarts,
                "realization_id": template.realization_id,
                "stereo_policy": template.stereo_policy,
                "template_id": template.template_id,
            }
            for template in sorted(
                result.templates,
                key=lambda item: (item.abstraction_level, item.template_id),
            )
        )
        core = row.get("reaction_core") or {}
        quality = core.get("quality") or {}
        found[reaction_id] = {
            "reaction_id": reaction_id,
            "status": "accepted" if templates else "rejected",
            "rejection_reason": result.rejection_reason,
            "rejection_stage": result.rejection_stage,
            "stored_core_status": str(quality.get("status") or "missing"),
            "stored_review_reasons": list(quality.get("review_reasons") or ()),
            "core_admission_reason": diagnostics.get("core_admission_reason"),
            "core_admission_issues": list(
                diagnostics.get("core_admission_issues") or ()
            ),
            "departure_edit_descriptors": list(
                diagnostics.get("core_admission_departure_edit_descriptors")
                or ()
            ),
            "round_trip_levels": list(
                diagnostics.get("round_trip_levels") or ()
            ),
            "templates": list(templates),
        }
    missing = sorted(requested_set.difference(found))
    if missing:
        raise ValueError(f"round-trip audit reactions not found: {missing}")
    source_files = tuple(source_shard_files(source))
    results = [found[reaction_id] for reaction_id in requested]
    report = {
        "artifact_type": "generic_operator_round_trip_audit",
        "schema_version": ROUND_TRIP_AUDIT_SCHEMA_VERSION,
        "algorithm_version": ROUND_TRIP_AUDIT_ALGORITHM_VERSION,
        "source": str(Path(source).resolve()),
        "source_files": [
            {
                "path": str(path.resolve()),
                "sha256": _sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for path in source_files
        ],
        "core_admission_policy": core_admission_policy,
        "core_admission_policy_definition_id": policy.definition_id,
        "requested_reaction_ids": list(requested),
        "accepted_reaction_count": sum(
            item["status"] == "accepted" for item in results
        ),
        "rejected_reaction_count": sum(
            item["status"] == "rejected" for item in results
        ),
        "results": results,
    }
    output = Path(destination)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return {**report, "output_path": str(output.resolve())}


__all__ = [
    "ROUND_TRIP_AUDIT_ALGORITHM_VERSION",
    "ROUND_TRIP_AUDIT_SCHEMA_VERSION",
    "audit_generic_round_trips",
]
