"""Present saved precedent evidence and explicitly labelled historical recovery."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from core_retrosynthesis.chemistry import canonical_smiles

from ..answers.condition_precedents import load_condition_precedent_evidence
from .condition_precedents import condition_precedent_view

if TYPE_CHECKING:
    from ..core.store import InvestigationStore

from ..answers.step_precedents import (
    _saved_support,
    _structure_key,
    load_step_precedent_evidence,
)


def _assessment_preview(candidate: dict[str, Any], key: tuple[str, str]) -> dict[str, Any]:
    """Expose saved source structures without claiming a detailed inspection occurred."""
    assessment = candidate["assessment"]
    unique = {}
    for item in assessment.get("precedent_matches", []):
        unique.setdefault((item["reaction_id"], item["reference_id"], item["mapped_reaction_smiles"]), item)
    matches = list(unique.values())
    records = [{
        **item, "reaction_smiles": f"{item['precursor_smiles']}>>{item['product_smiles']}",
        "same_recorded_precursors": canonical_smiles(item["precursor_smiles"]) == key[0],
        "same_recorded_product": canonical_smiles(item["product_smiles"]) == key[1],
        "support_kind": "template_precedent", "observations": [], "procedures": [],
        "limitations": ["Detailed inspection and experimental conditions were not attached to this answer."],
    } for item in matches[:3]]
    return {
        "artifact_ref": candidate["source_ref"], "source_ref": candidate["source_ref"],
        "evidence_origin": "saved_assessment", "inspection_status": "not_attached",
        "selection": {"precursor_smiles": key[0], "target_smiles": key[1], **candidate["selector"]},
        "status": "precedents_available" if matches else "no_precedents_retrieved",
        "scope": "saved_assessment_matches", "saved_match_count": len(matches),
        "assessment_status": assessment.get("status"), "assessment_warnings": assessment.get("warnings", []),
        "page": {"offset": 0, "returned": len(records), "next_offset": 3 if len(matches) > 3 else None},
        "distinct_references_on_page": len({item["reference_id"] for item in records if item["reference_id"]}),
        "precedents": records, "experimental_feasibility": "not_established",
    }

def answer_step_precedents(store: InvestigationStore, answer: dict[str, Any]) -> dict[str, list[dict]]:
    """Resolve saved answer links for display; preserve unavailable evidence explicitly."""
    molecules = {item["id"]: item["smiles"] for item in answer.get("molecules", [])}
    result = {}
    support, support_error = {}, None
    try:
        if any(not step.get("precedent_refs") for step in answer.get("steps", [])):
            support = _saved_support(store)
    except (OSError, ValueError, KeyError) as exc:
        support_error = str(exc)
    for step in answer.get("steps", []):
        records = []
        for reference in dict.fromkeys(step.get("precedent_refs", [])):
            try:
                record = load_step_precedent_evidence(
                    store, reference, ".".join(molecules[key] for key in step["reactant_ids"]),
                    ".".join(molecules[key] for key in step["product_ids"]),
                )
                records.append({"artifact_ref": reference, **record})
            except (OSError, ValueError, KeyError) as exc:
                records.append({"artifact_ref": reference, "status": "evidence_unavailable",
                                "error": str(exc), "precedents": []})
        if not step.get("precedent_refs"):
            if support_error:
                records.append({"status": "evidence_unavailable", "error": support_error, "precedents": []})
            key = _structure_key(".".join(molecules[item] for item in step["reactant_ids"]),
                                 ".".join(molecules[item] for item in step["product_ids"]))
            candidates = support.get(key, [])
            existing = next((item for item in candidates if item["inspection"] and item["available"]), None)
            if existing:
                records.append({"artifact_ref": existing["source_ref"], **existing["inspection"],
                                "evidence_origin": "recovered_inspection"})
            else:
                assessed = next((item for item in candidates if item["assessment"].get("precedent_matches")), None)
                if assessed:
                    records.append(_assessment_preview(assessed, key))
        for reference in dict.fromkeys(step.get("condition_precedent_refs", [])):
            try:
                inspected = load_condition_precedent_evidence(
                    store, reference, ".".join(molecules[key] for key in step["reactant_ids"]),
                    ".".join(molecules[key] for key in step["product_ids"]),
                )
                records.append(condition_precedent_view(inspected, reference))
            except (OSError, ValueError, KeyError, TypeError) as exc:
                records.append({"artifact_ref": reference, "status": "evidence_unavailable",
                                "error": str(exc), "precedents": []})
        result[step["id"]] = records
    return result
