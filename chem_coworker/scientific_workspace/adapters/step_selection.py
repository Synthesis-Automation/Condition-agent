"""Resolve exact saved step selections without new scientific retrieval."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from reactive_taxonomy.structure_audit import audit_structure

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

def _selection(
    payload: dict[str, Any], step_id: str | None, realization_id: str | None,
    strategy_id: str | None = None, precursor_smiles: str | None = None,
) -> tuple[dict, dict]:
    """Resolve a concrete candidate; template realization IDs alone may be shared."""
    result = payload["result"]
    if payload["operation"] == "disconnect_target":
        if step_id is not None or not realization_id:
            raise ValueError("Select a saved disconnect_target realization_id")
        expected = audit_structure(precursor_smiles) if precursor_smiles is not None else None
        if expected is not None and not expected.valid:
            raise ValueError("precursor_smiles must identify valid complete precursor graphs")
        candidates = [candidate for strategy in result.get("strategies", [])
                      if strategy_id is None or strategy.get("strategy_id") == strategy_id
                      for candidate in [strategy.get("representative"), *strategy.get("alternate_realizations", [])]
                      if candidate and candidate.get("realization_id") == realization_id
                      and (expected is None or audit_structure(candidate["precursor_smiles"]).canonical_smiles
                           == expected.canonical_smiles)]
        if not candidates:
            raise ValueError("No saved candidate matches realization_id, strategy_id and precursor_smiles")
        identities = {(audit_structure(item["precursor_smiles"]).canonical_smiles,
                       audit_structure(item["target_smiles"]).canonical_smiles,
                       item.get("operator_id"), item.get("template_id")) for item in candidates}
        if len(identities) != 1:
            raise ValueError("Ambiguous realization_id: supply strategy_id and/or precursor_smiles from the saved candidate")
        return candidates[0], {}
    if strategy_id is not None or precursor_smiles is not None:
        raise ValueError("strategy_id and precursor_smiles selectors are only valid for saved disconnections")
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
