"""Find or complete protocol context from a decisive condition core.

Completion follows an evidence ladder.  An observed exact core returns its own
protocol variants.  An unseen core may fall back to the nearest experimental
condition group, but group-level context is labeled as a broader prior and is
never represented as an observed exact protocol.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Tuple

from .condition_grouping_poc import (
    CONDITION_GROUPING_POC_DEFINITION_VERSION,
    CONDITION_GROUPING_POC_SCHEMA_VERSION,
    condition_group_assignment_status,
    parse_condition_display,
)


CONDITION_COMPLETION_POC_SCHEMA_VERSION = "1.0"


def _read_jsonl(path: Path) -> Iterable[Dict[str, Any]]:
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            value = line.strip()
            if not value:
                continue
            payload = json.loads(value)
            if not isinstance(payload, dict):
                raise ValueError(f"{path}:{line_number}: expected JSON object")
            yield payload


def _find_record(
    path: Path, field: str, value: str
) -> Optional[Dict[str, Any]]:
    for record in _read_jsonl(path):
        if str(record.get(field) or "") == value:
            return record
    return None


def _top_values(values: Any, top_k: int) -> Tuple[Dict[str, Any], ...]:
    if not isinstance(values, (list, tuple)):
        return ()
    return tuple(dict(value) for value in values[:top_k] if isinstance(value, dict))


def _validate_artifact_record(record: Mapping[str, Any]) -> None:
    if str(record.get("schema_version") or "") != CONDITION_GROUPING_POC_SCHEMA_VERSION:
        raise ValueError("condition grouping artifact schema version mismatch")
    if (
        str(record.get("definition_version") or "")
        != CONDITION_GROUPING_POC_DEFINITION_VERSION
    ):
        raise ValueError("condition grouping artifact definition version mismatch")


def _group_context(
    group: Mapping[str, Any], top_k: int
) -> Dict[str, Any]:
    context = group.get("context_cross_references")
    context = context if isinstance(context, Mapping) else {}
    return {
        "condition_group_id": str(group.get("condition_group_id") or ""),
        "prevalence_default_core_id": str(
            group.get("prevalence_default_core_id") or ""
        ),
        "prevalence_default_core_display": str(
            group.get("prevalence_default_core_display") or ""
        ),
        "statistical_prototype_core_id": str(
            group.get("statistical_prototype_core_id") or ""
        ),
        "statistical_prototype_core_display": str(
            group.get("statistical_prototype_core_display") or ""
        ),
        "core_role_consensus": dict(group.get("core_role_consensus") or {}),
        "most_common_cores": _top_values(group.get("most_common_cores"), top_k),
        "solvent_options": _top_values(context.get("solvent_systems"), top_k),
        "additive_options": _top_values(context.get("additive_systems"), top_k),
        "temperature_c": dict(context.get("temperature_c") or {}),
        "time_h": dict(context.get("time_h") or {}),
    }


def _exact_completion(
    query_display: str,
    core: Mapping[str, Any],
    group: Optional[Mapping[str, Any]],
    top_k: int,
) -> Dict[str, Any]:
    _validate_artifact_record(core)
    variants = _top_values(core.get("protocol_variants"), top_k)
    return {
        "artifact_type": "condition_core_completion_poc",
        "schema_version": CONDITION_COMPLETION_POC_SCHEMA_VERSION,
        "valid": True,
        "query_condition_display": query_display,
        "normalized_core_display": str(core.get("core_display") or ""),
        "condition_core_id": str(core.get("condition_core_id") or ""),
        "condition_group_id": str(core.get("condition_group_id") or ""),
        "completion_level": "exact_core",
        "assignment_status": "observed_exact_core",
        "centroid_similarity": core.get("centroid_similarity"),
        "assignment_margin": core.get("assignment_margin"),
        "core_observation_count": int(core.get("observation_count") or 0),
        "material_variant_count": int(core.get("material_variant_count") or 0),
        "suggested_protocol_variants": variants,
        "solvent_options": _top_values(core.get("solvent_systems"), top_k),
        "additive_options": _top_values(core.get("additive_systems"), top_k),
        "temperature_c": dict(core.get("temperature_c") or {}),
        "time_h": dict(core.get("time_h") or {}),
        "group_prior": _group_context(group, top_k) if group is not None else None,
        "warnings": (
            "OBSERVED_PROTOCOL_FREQUENCY_IS_NOT_OPTIMALITY",
            "RUN_COMPATIBILITY_FILTERS_BEFORE_RECOMMENDATION",
        ),
    }


def _learned_group_assignment(
    feature_tokens: Tuple[str, ...], model_path: Path
) -> Tuple[str, float, float, str]:
    try:
        import joblib
        import numpy as np
        from sklearn.preprocessing import normalize
    except ImportError as exc:  # pragma: no cover - optional POC environment
        raise RuntimeError(
            "condition completion POC requires scikit-learn, numpy, and joblib"
        ) from exc
    artifact = joblib.load(model_path)
    if artifact.get("schema_version") != CONDITION_GROUPING_POC_SCHEMA_VERSION:
        raise ValueError("condition grouping model schema version mismatch")
    if artifact.get("definition_version") != CONDITION_GROUPING_POC_DEFINITION_VERSION:
        raise ValueError("condition grouping model definition version mismatch")
    vectorizer = artifact["vectorizer"]
    svd = artifact["svd"]
    cluster_model = artifact["cluster_model"]
    matrix = vectorizer.transform([feature_tokens])
    if matrix.nnz == 0:
        raise ValueError("query core has no recognized grouping features")
    latent = normalize(svd.transform(matrix), norm="l2")
    centers = normalize(cluster_model.cluster_centers_, norm="l2")
    scores = np.asarray(latent @ centers.T)[0]
    order = np.argsort(scores)[::-1]
    label = int(order[0])
    similarity = float(scores[label])
    margin = similarity - float(scores[int(order[1])])
    group_id = str(artifact["group_id_by_label"][label])
    status = condition_group_assignment_status(similarity, margin)
    return group_id, similarity, margin, status


def complete_condition_core(
    condition_display: str,
    artifact_dir: str | Path,
    *,
    top_k: int = 5,
) -> Dict[str, Any]:
    """Complete solvent and protocol context from an exact or nearby core."""

    if top_k < 1:
        raise ValueError("top_k must be positive")
    destination = Path(artifact_dir)
    cores_path = destination / "condition_cores.jsonl"
    groups_path = destination / "condition_groups.jsonl"
    model_path = destination / "condition_grouping_model.joblib"
    for path in (cores_path, groups_path, model_path):
        if not path.is_file():
            raise FileNotFoundError(f"missing condition grouping artifact: {path}")

    parsed = parse_condition_display(condition_display)
    core_id = parsed.core_id
    if core_id is None:
        return {
            "artifact_type": "condition_core_completion_poc",
            "schema_version": CONDITION_COMPLETION_POC_SCHEMA_VERSION,
            "valid": False,
            "query_condition_display": condition_display,
            "normalized_core_display": "",
            "condition_core_id": None,
            "condition_group_id": None,
            "completion_level": "none",
            "error": "NO_DECISIVE_CORE_COMPONENT",
            "warnings": tuple(parsed.warnings),
        }

    core_record = _find_record(cores_path, "condition_core_id", core_id)
    if core_record is not None:
        group_id = str(core_record.get("condition_group_id") or "")
        group_record = _find_record(groups_path, "condition_group_id", group_id)
        if group_record is not None:
            _validate_artifact_record(group_record)
        return _exact_completion(
            condition_display,
            core_record,
            group_record,
            top_k,
        )

    group_id, similarity, margin, status = _learned_group_assignment(
        parsed.core_feature_tokens,
        model_path,
    )
    group_record = _find_record(groups_path, "condition_group_id", group_id)
    if group_record is None:
        raise ValueError("learned condition group is missing from group catalog")
    _validate_artifact_record(group_record)
    group_context = _group_context(group_record, top_k)
    return {
        "artifact_type": "condition_core_completion_poc",
        "schema_version": CONDITION_COMPLETION_POC_SCHEMA_VERSION,
        "valid": True,
        "query_condition_display": condition_display,
        "normalized_core_display": parsed.core_display,
        "condition_core_id": core_id,
        "condition_group_id": group_id,
        "completion_level": "learned_group",
        "assignment_status": status,
        "centroid_similarity": round(similarity, 6),
        "assignment_margin": round(margin, 6),
        "core_observation_count": 0,
        "material_variant_count": 0,
        "suggested_protocol_variants": (),
        "solvent_options": group_context["solvent_options"],
        "additive_options": group_context["additive_options"],
        "temperature_c": group_context["temperature_c"],
        "time_h": group_context["time_h"],
        "group_prior": group_context,
        "warnings": (
            "UNSEEN_CORE_GROUP_LEVEL_CONTEXT_ONLY",
            "DO_NOT_PRESENT_GROUP_PRIOR_AS_OBSERVED_EXACT_PROTOCOL",
            "RUN_COMPATIBILITY_FILTERS_BEFORE_RECOMMENDATION",
        ),
    }


__all__ = [
    "CONDITION_COMPLETION_POC_SCHEMA_VERSION",
    "complete_condition_core",
]
