"""Validate exact step attribution and require available supporting evidence."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING, Any

from core_retrosynthesis.chemistry import canonical_smiles

if TYPE_CHECKING:
    from ..core.store import InvestigationStore

from ..adapters.step_selection import (
    SCHEMA_VERSION,
    SOURCE_OPERATIONS,
    _call,
    _selection,
)


def load_step_precedent_evidence(
    store: InvestigationStore, reference: str, precursor_smiles: str, target_smiles: str,
) -> dict[str, Any]:
    """Resolve a recorded inspection and reject citations belonging to another step."""
    record = _call(store, reference, {"inspect_step_precedents"})["result"]
    if record.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Unsupported step precedent evidence schema")
    selection = record["selection"]
    for actual, expected in ((precursor_smiles, selection["precursor_smiles"]),
                             (target_smiles, selection["target_smiles"])):
        canonical = canonical_smiles(actual)
        if not canonical or canonical != expected:
            raise ValueError("Precedent inspection does not match this answer step's structures and stereochemistry")
    return record

def _structure_key(precursors: str, target: str) -> tuple[str, str] | None:
    left, right = canonical_smiles(precursors), canonical_smiles(target)
    return (left, right) if left and right else None

def _saved_support(store: InvestigationStore) -> dict[tuple[str, str], list[dict[str, Any]]]:
    """Index verified saved selections by both sides, retaining stereo and provenance.

    This reads call artifacts only. It does not retrieve from a live corpus, rerun
    chemistry, or infer that an agent reviewed a retrieved match.
    """
    support: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for event in reversed(store.events()):
        if event.kind != "call":
            continue
        try:
            payload = store.read_artifact(event.artifact_ref)
            operation = payload.get("operation")
            if payload.get("execution_status") != "completed":
                continue
            result = payload.get("result", {})
            selections = []
            if operation == "inspect_step_precedents" and result.get("schema_version") == SCHEMA_VERSION:
                selections.append((result["selection"], {}, {}, result))
            elif operation == "disconnect_target":
                for strategy in result.get("strategies", []):
                    for selected in [strategy.get("representative"), *strategy.get("alternate_realizations", [])]:
                        if selected:
                            selections.append((selected, {}, {"realization_id": selected["realization_id"]}, None))
            elif operation in SOURCE_OPERATIONS:
                selectors = [{}] if operation == "assess_route_step" else [
                    {"step_id": item["external_step_id"]} for item in result.get("proposal", {}).get("steps", [])
                ]
                for selector in selectors:
                    selected, assessment = _selection(payload, selector.get("step_id"), None)
                    selections.append((selected, assessment, selector, None))
            for selected, assessment, selector, inspection in selections:
                key = _structure_key(selected["precursor_smiles"], selected["target_smiles"])
                if key:
                    support.setdefault(key, []).append({
                        "source_ref": event.artifact_ref, "selector": selector,
                        "assessment": assessment, "inspection": inspection,
                        "available": bool((inspection or {}).get("precedents") or assessment.get("precedent_matches")
                                          or selected.get("precedent_reaction_ids")),
                    })
        except (OSError, ValueError, KeyError, TypeError):
            # Unrelated malformed/corrupt calls cannot become inferred support.
            # Explicitly cited artifacts are validated separately and fail closed.
            continue
    return support

def require_available_step_precedents(store: InvestigationStore, answer: dict[str, Any]) -> None:
    """Reject publication that silently omits available evidence for a final step."""
    if not answer.get("steps"):
        return
    support = _saved_support(store)
    molecules = {item["id"]: item["smiles"] for item in answer["molecules"]}
    missing = []
    for step in answer["steps"]:
        precursors = ".".join(molecules[key] for key in step["reactant_ids"])
        target = ".".join(molecules[key] for key in step["product_ids"])
        candidates = [item for item in support.get(_structure_key(precursors, target), []) if item["available"]]
        if not candidates:
            continue
        if any(load_step_precedent_evidence(store, reference, precursors, target).get("precedents")
               for reference in step.get("precedent_refs", [])):
            continue
        existing = next((item for item in candidates if item["inspection"]), None)
        if existing:
            action = f"attach existing inspection {existing['source_ref']} in precedent_refs"
        else:
            candidate = candidates[0]
            arguments = {"source_ref": candidate["source_ref"], **candidate["selector"], "limit": 3}
            action = f"run inspect_step_precedents with {json.dumps(arguments)}, review it, and attach its artifact_ref in precedent_refs"
        missing.append(f"{step['id']}: {action}")
    if missing:
        raise ValueError("Available supporting reactions must accompany the final presented steps. " + "; ".join(missing))
