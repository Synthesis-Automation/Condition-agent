"""Resolve exact saved step selections without new scientific retrieval."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from pydantic import BaseModel, ConfigDict, Field

from reactive_taxonomy.structure_audit import audit_structure

if TYPE_CHECKING:
    from ..core.store import InvestigationStore


SCHEMA_VERSION = "route_step_precedents.v1"

SOURCE_OPERATIONS = {"disconnect_target", "assess_route_step", "assess_route_proposal", "revise_route_branch",
                     "assess_retro_validity"}


class SavedCandidateSelection(BaseModel):
    """Select recorded operator evidence without manually copying an atom mapping."""

    model_config = ConfigDict(extra="forbid", strict=True)
    source_ref: str = Field(pattern=r"^sha256:[0-9a-f]{64}$")
    realization_id: str = Field(min_length=1)
    strategy_id: str | None = None


def resolve_saved_candidate(
    store: InvestigationStore, proposal: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    """Copy only an explicitly selected mapping; canonical assessment still validates it.

    This is proposed operator reconstruction, never observed source atom mapping.
    Explicit mappings and candidate mappings cannot silently override one another.
    """
    if not isinstance(proposal, dict) or "saved_candidate" not in proposal:
        return proposal, None
    value = dict(proposal)
    selector = SavedCandidateSelection.model_validate(value.pop("saved_candidate"))
    if value.get("mapped_reaction_smiles") is not None:
        raise ValueError("Supply saved_candidate or mapped_reaction_smiles, not both; conflicting mapping sources require separate assessments")
    payload = _call(store, selector.source_ref, {"disconnect_target"})
    selected, _ = _selection(payload, None, selector.realization_id, selector.strategy_id,
                             value.get("precursor_smiles"))
    for name in ("precursor_smiles", "target_smiles"):
        actual, expected = audit_structure(value.get(name, "")), audit_structure(selected[name])
        if not actual.valid or not expected.valid or actual.canonical_smiles != expected.canonical_smiles:
            raise ValueError(f"saved_candidate does not match {name}, including stereochemistry and chemical form")
    mapping = selected.get("condition_query_reaction_smiles")
    if not mapping:
        raise ValueError("Saved candidate has no mapped reconstruction; assess the explicit structures without saved_candidate")
    value["mapped_reaction_smiles"] = mapping
    return value, {**selector.model_dump(), "origin": "saved_operator_reconstruction",
                   "validation": "revalidated_by_proposal_assessor", "observed_reaction": False}

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
                       item.get("operator_id"), item.get("template_id"),
                       item.get("condition_query_reaction_smiles")) for item in candidates}
        if len(identities) != 1:
            raise ValueError("Ambiguous realization_id or conflicting saved mappings: supply strategy_id and/or precursor_smiles from the saved candidate")
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
