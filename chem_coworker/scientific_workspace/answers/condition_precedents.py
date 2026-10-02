"""Validate condition-inspection attribution against exact answer structures."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from core_retrosynthesis.chemistry import canonical_smiles

if TYPE_CHECKING:
    from ..core.store import InvestigationStore


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
