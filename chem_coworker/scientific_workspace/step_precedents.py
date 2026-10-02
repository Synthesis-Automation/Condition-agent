"""Recorded inspection of template precedents, linked to an exact proposed step.

Scientific matching stays in the canonical domain packages. This adapter joins
saved selections and baseline-pinned source records without upgrading admission.
"""

from __future__ import annotations

import json
from typing import Any, TYPE_CHECKING

from core_retrosynthesis.chemistry import canonical_smiles
from core_retrosynthesis.step_precedents import lookup_reaction_precedents

from .source_catalogs import reference_records
from .condition_precedents import condition_precedent_view, load_condition_precedent_evidence

if TYPE_CHECKING:
    from .operations import ScientificOperations
    from .store import InvestigationStore


SCHEMA_VERSION = "route_step_precedents.v1"
SOURCE_OPERATIONS = {"disconnect_target", "assess_route_step", "assess_route_proposal", "revise_route_branch"}


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
    if payload["operation"] == "assess_route_step":
        if step_id is not None:
            raise ValueError("A single-step assessment does not require step_id")
        return result["proposal"], result["assessment"]
    proposals = [step for step in result["proposal"]["steps"] if step["external_step_id"] == step_id]
    if not proposals:
        raise ValueError("step_id is not present in this saved route")
    assessment = next(item["assessment"] for item in result["assessment"]["step_assessments"]
                      if item["external_step_id"] == step_id)
    return proposals[0], assessment


def inspect_step_precedents(
    operations: ScientificOperations, source_ref: str, step_id: str | None = None,
    realization_id: str | None = None, offset: int = 0, limit: int = 3,
) -> dict[str, Any]:
    """Inspect at most five reactions on a page from one saved step's support.

    Disconnections use the canonical top-20 template lookup; assessments retain
    their saved bounded match set. Counts never imply whole-corpus coverage.
    """
    if type(offset) is not int or offset < 0 or type(limit) is not int or not 1 <= limit <= 5:
        raise ValueError("offset must be nonnegative and limit must be between one and five")
    payload = _call(operations.store, source_ref, SOURCE_OPERATIONS)
    selected, assessment = _selection(payload, step_id, realization_id)
    target = canonical_smiles(selected["target_smiles"])
    precursors = canonical_smiles(selected["precursor_smiles"])
    if not target or not precursors:
        raise ValueError("Selected step requires valid reactant and product structures")
    library = operations._external_route_library()
    if realization_id:
        lookup = lookup_reaction_precedents(
            step_id=realization_id, template_id=selected["template_id"], operator_id=selected["operator_id"],
            product_smiles=target, precursor_smiles=precursors, library=library, limit=20,
        )
        matches = [item.to_dict() for item in lookup.matches]
        available = lookup.available_precedent_count
        scope = "selected_template_top_20"
    else:
        matches = assessment.get("precedent_matches", [])
        available = None
        scope = "saved_assessment_matches"
    # A repeated record through several operators is one supporting reaction here.
    unique = {}
    for match in matches:
        key = (match["reaction_id"], match["reference_id"], match["mapped_reaction_smiles"])
        unique.setdefault(key, match)
    matches = list(unique.values())
    if offset > len(matches):
        raise ValueError("offset exceeds the saved precedent selection")
    page = matches[offset:offset + limit]
    references, reference_status = reference_records(operations, {item["reference_id"] for item in page})
    ids = list(dict.fromkeys(item["reaction_id"] for item in page))
    observations, procedures = {}, {}
    if ids:
        try:
            observations = operations.get_precedents(ids, limit=20)
        except FileNotFoundError:
            observations = {"records": [], "availability": "index_unavailable"}
        try:
            procedures = operations.get_procedures(ids)
        except FileNotFoundError:
            procedures = {"records": [], "availability": "catalog_unavailable"}
    templates = {item.template_id: item for item in library.templates}
    records = []
    for match in page:
        template = templates.get(match["template_id"])
        try:
            comparison = operations.compare_molecules(target, match["product_smiles"], timeout_seconds=1).to_dict()
        except ValueError as exc:
            comparison = {"status": "unavailable", "warnings": [str(exc)]}
        rows = [{key: row.get(key) for key in (
            "reaction_id", "observation_id", "reference_id", "reaction_smiles", "resolved_recipe", "yield_pct",
            "condition_uncertain", "chemistry_status", "condition_status",
        )} for row in observations.get("records", []) if row["reaction_id"] == match["reaction_id"]]
        source_procedures = [row for row in procedures.get("records", [])
                             if row["reaction_id"] == match["reaction_id"]]
        records.append({
            **match, "reaction_smiles": f"{match['precursor_smiles']}>>{match['product_smiles']}",
            "support_kind": "template_precedent", "reference_record": references.get(match["reference_id"]),
            "same_recorded_product": canonical_smiles(match["product_smiles"]) == target,
            "same_recorded_precursors": canonical_smiles(match["precursor_smiles"]) == precursors,
            "product_comparison": comparison,
            "template_context": {"edit_tokens": list(template.edit_tokens), "handle_signature": template.handle_signature,
                                 "stereo_policy": template.stereo_policy} if template else None,
            "observations": rows, "procedures": source_procedures,
            "experimental_link_scope": "reaction_id_only_template_has_no_observation_id",
            "limitations": ["Template membership and similarity do not establish transfer or experimental feasibility.",
                            "Associated conditions and yields belong to their source observations, not the proposed step.",
                            "Substrate, functional-group and stereochemical transfer require inspection."],
        })
    return {
        "schema_version": SCHEMA_VERSION, "source_ref": source_ref,
        "selection": {"step_id": step_id, "realization_id": realization_id,
                      "target_smiles": target, "precursor_smiles": precursors},
        "status": "precedents_available" if matches else "no_precedents_retrieved",
        "assessment_status": assessment.get("status"), "assessment_warnings": assessment.get("warnings", []),
        "scope": scope, "available_template_records": available, "saved_match_count": len(matches),
        "retrieval_truncated": available > 20 if available is not None else None,
        "page": {"offset": offset, "limit": limit, "returned": len(page),
                 "next_offset": offset + limit if offset + limit < len(matches) else None},
        "distinct_references_on_page": len({item["reference_id"] for item in page if item["reference_id"]}),
        "reference_catalog_status": reference_status,
        "observation_page": {key: value for key, value in observations.items() if key != "records"},
        "procedure_catalog_status": procedures.get("availability", "catalog_available" if ids else "not_requested"),
        "precedents": records, "experimental_feasibility": "not_established",
        "limitations": ["Counts describe the selected template or saved assessment, not an exhaustive literature search.",
                        "An empty result does not establish that no experimental precedent exists.",
                        "This inspection does not change route admission or validate proposed conditions."],
    }


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
