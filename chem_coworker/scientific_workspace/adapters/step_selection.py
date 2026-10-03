"""Resolve exact saved step selections without new scientific retrieval."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from ..core.store import InvestigationStore


SCHEMA_VERSION = "route_step_precedents.v1"

SOURCE_OPERATIONS = {"disconnect_target", "assess_route_step", "assess_route_proposal", "revise_route_branch",
                     "assess_retro_validity"}

def _call(store: InvestigationStore, reference: str, operations: set[str]) -> dict[str, Any]:
    if not any(event.artifact_ref == reference and event.kind == "call" for event in store.events()):
        raise ValueError("Precedent source must be a recorded scientific call")
    payload = store.read_artifact(reference)
    if payload.get("operation") not in operations or payload.get("execution_status") != "completed":
        raise ValueError("Precedent source must be a completed supported scientific call")
    return payload

def _selection(payload: dict[str, Any], step_id: str | None, realization_id: str | None) -> tuple[dict, dict]:
    result = payload["result"]
    if payload["operation"] == "disconnect_target":
        if step_id is not None or not realization_id:
            raise ValueError("Select a saved disconnect_target realization_id")
        candidates = [candidate for strategy in result.get("strategies", [])
                      for candidate in [strategy.get("representative"), *strategy.get("alternate_realizations", [])]
                      if candidate and candidate.get("realization_id") == realization_id]
        if not candidates:
            raise ValueError("realization_id is not present in this saved disconnection")
        return candidates[0], {}
    if realization_id is not None:
        raise ValueError("realization_id is only valid for a saved disconnection")
    if payload["operation"] in {"assess_route_step", "assess_retro_validity"}:
        if step_id is not None:
            raise ValueError("A single-step assessment does not require step_id")
        assessment = (result["validity"]["structural_assessment"]
                      if payload["operation"] == "assess_retro_validity" else result["assessment"])
        return result["proposal"], assessment
    proposals = [step for step in result["proposal"]["steps"] if step["external_step_id"] == step_id]
    if not proposals:
        raise ValueError("step_id is not present in this saved route")
    assessment = next(item["assessment"] for item in result["assessment"]["step_assessments"]
                      if item["external_step_id"] == step_id)
    return proposals[0], assessment
