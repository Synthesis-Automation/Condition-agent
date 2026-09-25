"""Application records for attributed investigator proposals, never chemistry rules."""

from __future__ import annotations

from dataclasses import asdict
from typing import Any, TYPE_CHECKING

from condition_registry import CONDITION_RECIPE_COMPONENT_BUCKETS

from .store import canonical_bytes

if TYPE_CHECKING:
    from .operations import ScientificOperations


def record_adaptation(
    operations: ScientificOperations, source_ref: str, observation_id: str,
    components: list[dict[str, Any]], operating_conditions: dict[str, Any],
    change_reasons: dict[str, str], evidence_refs: list[str],
    assumptions: list[str], risks: list[str],
) -> dict[str, Any]:
    """Compare complete canonical recipes and preserve proposal/support separately."""
    source = operations.store.read_artifact(source_ref)
    if source.get("operation") != "inspect_condition_precedents" or source.get("execution_status") != "completed":
        raise ValueError("source_ref must identify completed precedent inspection")
    selected = [row for row in source["result"]["precedents"] if row["observation"]["observation_id"] == observation_id]
    if len(selected) != 1:
        raise ValueError("Select exactly one inspected observation_id")
    for values in (evidence_refs, assumptions, risks):
        if not isinstance(values, list) or not all(isinstance(value, str) and value.strip() for value in values):
            raise ValueError("Evidence, assumptions and risks must be lists of nonempty strings")
    if not evidence_refs or not assumptions or not risks:
        raise ValueError("Adaptations require explicit evidence, assumptions and unresolved risks")
    kinds = {event.artifact_ref: event.kind for event in operations.store.events()}
    for reference in evidence_refs:
        operations.store.read_artifact(reference)
        if kinds.get(reference) not in {"call", "derived_file", "replay", "custom_execution"}:
            raise ValueError("Adaptation evidence must be recorded scientific evidence")
    allowed = ("temperature_c", "time_h", "concentration_m", "atmosphere")
    if set(operating_conditions) - set(allowed):
        raise ValueError("Unsupported operating condition field")
    proposed = asdict(operations.resolve_recipe(components, **operating_conditions))
    original = selected[0]["observation"]["resolved_recipe"]
    changes = []
    for key in (*CONDITION_RECIPE_COMPONENT_BUCKETS, *allowed):
        empty = [] if key in CONDITION_RECIPE_COMPONENT_BUCKETS else None
        before, after = original.get(key, empty), proposed.get(key, empty)
        if canonical_bytes(before) != canonical_bytes(after):
            reason = change_reasons.get(key)
            if not isinstance(reason, str) or not reason.strip():
                raise ValueError(f"Missing change reason for {key}")
            changes.append({"field": key, "before": before, "after": after, "reason": reason,
                            "basis": "proposed", "evidence_refs": evidence_refs})
    if not changes or set(change_reasons) != {change["field"] for change in changes}:
        raise ValueError("Reasons must describe exactly the changed recipe fields")
    # Unsupported stage editing must not silently discard a reported protocol.
    if original.get("stages") or original.get("declared_absences"):
        raise ValueError("Stage/declared-absence adaptation is unsupported; preserve the source protocol")
    assessment = operations.assess_recipe(source["result"]["query_reaction_smiles"], proposed)
    return {
        "schema_version": "condition_adaptation.v1", "origin": "agent_proposal",
        "review_status": "unreviewed", "transfer_status": "not_established",
        "source_ref": source_ref, "observation_id": observation_id,
        "reaction_smiles": source["result"]["query_reaction_smiles"],
        "original_recipe": original, "proposed_recipe": proposed, "changes": changes,
        "assumptions": assumptions, "risks": risks, "evidence_refs": evidence_refs,
        "compatibility": asdict(assessment),
        "limitations": ["Citation existence does not establish support for a proposed change.",
                        "No predicted yield or experimentally validated transfer is implied."],
    }
