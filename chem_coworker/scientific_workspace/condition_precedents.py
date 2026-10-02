"""Bind saved condition inspections to explicit answer reactions for presentation."""

from __future__ import annotations

from typing import Any, TYPE_CHECKING

from core_retrosynthesis.chemistry import canonical_smiles

if TYPE_CHECKING:
    from .store import InvestigationStore


def load_condition_precedent_evidence(
    store: InvestigationStore, reference: str, precursor_smiles: str, target_smiles: str,
) -> dict[str, Any]:
    """Reject unrelated, incomplete or unrecorded condition inspections.

    Compare both complete sides, preserving specified stereo and chemical form.
    This validates attribution only; it does not establish experimental transfer.
    """
    if not any(event.artifact_ref == reference and event.kind == "call" for event in store.events()):
        raise ValueError("Condition precedent source must be a recorded scientific call")
    payload = store.read_artifact(reference)
    if (payload.get("operation") != "inspect_condition_precedents"
            or payload.get("execution_status") != "completed"):
        raise ValueError("Condition precedents require a completed inspect_condition_precedents call")
    result = payload["result"]
    if result.get("schema_version") != "condition_evidence_comparison.v1":
        raise ValueError("Unsupported condition precedent evidence schema")
    query = result.get("query_reaction_smiles", "")
    sides = query.split(">")
    if len(sides) != 3 or query != payload.get("arguments", {}).get("reaction_smiles"):
        raise ValueError("Condition inspection query does not match its recorded input")
    for actual, expected in ((precursor_smiles, sides[0]), (target_smiles, sides[2])):
        canonical = canonical_smiles(actual)
        if not canonical or canonical != canonical_smiles(expected):
            raise ValueError("Condition inspection does not match this answer step's structures and stereochemistry")
    return result


def condition_precedent_view(result: dict[str, Any], reference: str) -> dict[str, Any]:
    """Adapt saved observations to the shared precedent display without new chemistry."""
    records = []
    for entry in result.get("precedents", []):
        observation = entry["observation"]
        records.append({
            "match_id": observation["observation_id"],
            "reaction_id": observation["reaction_id"],
            "reference_id": observation.get("reference_id"),
            "reference_record": entry.get("reference_record"),
            "reaction_smiles": observation["reaction_smiles"],
            "support_kind": "condition_observation",
            "observations": [observation],
            "procedures": entry.get("procedure_observations", []),
            "reaction_level_procedures": entry.get("reaction_level_procedures", []),
            "experimental_link_scope": "exact_observation_id",
            "structural_comparison": entry.get("structural_comparison"),
            "compatibility": entry.get("compatibility"),
            "missing_operating_fields": entry.get("missing_operating_fields", []),
            "limitations": ["Conditions and yield belong to this source observation; transfer remains unverified."],
        })
    return {
        "artifact_ref": reference, "evidence_origin": "condition_inspection",
        "status": "precedents_available" if records else "no_precedents_retrieved",
        "scope": "selected_condition_observations", "saved_match_count": len(records),
        "distinct_references_on_page": result.get("distinct_reference_count", 0),
        "page": result.get("page", {}), "precedents": records,
        "limitations": result.get("limitations", []),
    }
