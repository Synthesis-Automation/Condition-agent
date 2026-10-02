"""Display saved condition observations without transferring their reported outcomes."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    pass


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
