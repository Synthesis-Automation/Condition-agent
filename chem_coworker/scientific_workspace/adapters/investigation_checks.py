"""Compose saved source searches and canonical chemistry checks for investigations."""

from __future__ import annotations

from dataclasses import asdict
import heapq
import math
import re
from typing import TYPE_CHECKING, Any, Iterator

from condition_registry import ConditionComponentInput, build_resolved_recipe_from_inputs
from condition_registry.models import ConditionProcessStage
from reactive_taxonomy.structure_audit import audit_structure

from .literature import load_source_passage
from .route_investigation import _record

if TYPE_CHECKING:
    from .operations import ScientificOperations
    from ..core.store import InvestigationStore


def load_proposed_recipe_check(
    store: InvestigationStore, reference: str, precursors: str, target: str,
) -> dict[str, Any]:
    """Require a completed recipe check for the exact graphs, chemical form and stereo."""
    from .step_selection import _call

    result = _call(store, reference, {"assess_proposed_recipe"})["result"]
    sides = result["reaction_smiles"].split(">")
    if len(sides) != 3:
        raise ValueError("Recipe check must contain explicit reaction sides")
    for saved, expected in ((sides[0], precursors), (sides[2], target)):
        left, right = audit_structure(saved), audit_structure(expected)
        if not left.valid or not right.valid or left.canonical_smiles != right.canonical_smiles:
            raise ValueError("Proposed recipe check does not match this step's structures and stereochemistry")
    return result


def _literal_matches(text: str, query: str) -> Iterator[tuple[int, str, int]]:
    """Yield ordered offsets without materializing matches or their excerpts."""
    for match in re.finditer(re.escape(query), text, re.IGNORECASE):
        yield match.start(), query, match.end()


def search_captured_sources(
    operations: ScientificOperations, queries: list[str], source_refs: list[str] | None,
    offset: int, limit: int,
) -> dict[str, Any]:
    """Search literal terms in recorded text; matches are leads, not chemical assignments."""
    if (not isinstance(queries, list) or not 1 <= len(queries) <= 10
            or any(not isinstance(q, str) or not q.strip() or len(q) > 500 for q in queries)):
        raise ValueError("queries must contain 1..10 literal nonempty terms of at most 500 characters")
    if type(offset) is not int or offset < 0 or type(limit) is not int or not 1 <= limit <= 10:
        raise ValueError("offset must be nonnegative; limit must be 1..10")
    refs = source_refs if source_refs is not None else [
        event.artifact_ref for event in operations.store.events() if event.kind == "literature_source"
        and operations.store.read_artifact(event.artifact_ref).get("extraction", {}).get("text")
    ]
    if not isinstance(refs, list) or len(refs) > 100 or any(not isinstance(ref, str) for ref in refs):
        raise ValueError("source_refs must be a list of at most 100 recorded sources")
    matches, roots = [], []
    total = 0
    for reference in dict.fromkeys(refs):
        source, text, root = load_source_passage(operations.store, reference)
        if root in roots:
            continue
        roots.append(root)
        # Searches always use the complete captured root, including upstream sections.
        text = source["extraction"]["text"]
        streams = (_literal_matches(text, query) for query in dict.fromkeys(queries))
        for match_start, query, match_end in heapq.merge(*streams):
            if offset <= total < offset + limit:
                start, end = max(0, match_start - 200), min(len(text), match_end + 400)
                matches.append({"source_ref": root, "source_url": source["source_url"],
                                "query": query, "match_start": match_start, "start": start, "end": end,
                                "text": text[start:end], "claim_support": "not_assessed"})
            total += 1
    if offset > total:
        raise ValueError("offset exceeds the captured match count")
    return {"schema_version": "captured_source_search.v1", "matches": matches,
            "total": total, "offset": offset,
            "next_offset": offset + limit if offset + limit < total else None,
            "evidence_refs": roots, "scope": "captured_text_only",
            "limitations": ["Literal matches do not establish compound identity or a preparative reaction.",
                            "Missing/OCR-corrupted text and uncaptured pages are not searched."]}


def inspect_route_inputs(
    operations: ScientificOperations, source_ref: str, leaf_queries: list[dict[str, Any]],
    source_refs: list[str] | None,
) -> dict[str, Any]:
    """Inspect every recorded route leaf, preserving assumptions and source-search gaps."""
    route = _record(operations, source_ref)
    leaves = route["assessment"]["leaf_smiles"]
    identities = {audit_structure(smiles).canonical_smiles: smiles for smiles in leaves}
    if not isinstance(leaf_queries, list) or len(leaf_queries) > 100:
        raise ValueError("leaf_queries must contain at most 100 {smiles, terms} objects")
    queries = {}
    for item in leaf_queries:
        if not isinstance(item, dict) or set(item) != {"smiles", "terms"}:
            raise ValueError("Each leaf query needs exactly smiles and terms")
        audit = audit_structure(item["smiles"])
        if not audit.valid or audit.canonical_smiles not in identities or audit.canonical_smiles in queries:
            raise ValueError("Each query must identify one distinct actual route leaf, including stereochemistry")
        queries[audit.canonical_smiles] = item["terms"]
    results, refs = [], [source_ref]
    unavailable = route.get("material_constraints", {}).get("unavailable_starting_materials", [])
    for identity, smiles in identities.items():
        previous = None
        for event in reversed(operations.store.events()):
            if event.kind != "call":
                continue
            call = operations.store.read_artifact(event.artifact_ref)
            if (call.get("operation") == "assess_starting_material" and call.get("execution_status") == "completed"
                    and set(call["arguments"]) <= {"smiles", "unavailable_starting_materials"}
                    and call["arguments"].get("unavailable_starting_materials", []) == unavailable
                    and audit_structure(call["arguments"]["smiles"]).canonical_smiles == identity):
                previous = event
                break
        assessment = (operations.store.read_artifact(previous.artifact_ref)["result"] if previous
                      else operations.assess_starting_material(smiles, unavailable_starting_materials=unavailable))
        if previous:
            refs.append(previous.artifact_ref)
        search = search_captured_sources(operations, queries[identity], source_refs, 0, 3) if identity in queries else None
        if search:
            refs.extend(search["evidence_refs"])
        results.append({"smiles": smiles, "starting_material_assessment": assessment,
                        "assessment_ref": previous.artifact_ref if previous else None,
                        "captured_source_search": search,
                        "source_search_status": "searched" if search is not None else "terms_not_supplied"})
    return {"schema_version": "route_input_inspection.v1", "source_ref": source_ref,
            "leaves": results, "scope": "route_from_stated_inputs",
            "evidence_refs": list(dict.fromkeys(refs)),
            "limitations": ["Molecular-weight stopping and registry membership do not establish availability.",
                            "Source matches require graph assignment and step assessment before route extension.",
                            "This operation does not automatically expand, admit or rank the route."]}


def assess_proposed_recipe(
    operations: ScientificOperations, reaction_smiles: str, components: list[dict[str, Any]],
    operating_conditions: dict[str, Any], evidence_refs: list[str] | None,
) -> dict[str, Any]:
    """Resolve a complete supplied recipe and delegate compatibility to the canonical package."""
    allowed = {"temperature_c", "time_h", "concentration_m", "atmosphere", "stages", "declared_absences"}
    if not isinstance(operating_conditions, dict) or set(operating_conditions) - allowed:
        raise ValueError("operating_conditions supports temperature_c, time_h, concentration_m, atmosphere, stages, declared_absences")
    if not isinstance(components, list) or not 1 <= len(components) <= 100:
        raise ValueError("Supply 1..100 explicit complete condition components")
    for index, item in enumerate(components):
        if not isinstance(item, dict) or not isinstance(item.get("raw_identifier"), str) or not item["raw_identifier"].strip():
            raise ValueError("Every condition component needs a nonempty raw_identifier")
        if not isinstance(item.get("provenance", {}), dict):
            raise ValueError(f"components[{index}].provenance must be a JSON object, e.g. {{'description': 'Source or proposal'}}")
    values = dict(operating_conditions)
    stages = values.get("stages", [])
    if not isinstance(stages, list) or len(stages) > 30:
        raise ValueError("stages must contain at most 30 explicit process stages")
    for index, stage in enumerate(stages):
        if isinstance(stage, dict) and not isinstance(stage.get("provenance", {}), dict):
            raise ValueError(f"stages[{index}].provenance must be a JSON object, e.g. {{'description': 'Proposed stage'}}")
    for item in [values, *stages]:
        if not isinstance(item, dict):
            raise ValueError("Process stages must be objects")
        for key in ("temperature_c", "time_h", "concentration_m"):
            value = item.get(key)
            if value is not None and (type(value) not in (int, float) or not math.isfinite(value)):
                raise ValueError(f"{key} must be a finite number or null")
            if value is not None and key != "temperature_c" and value < 0:
                raise ValueError(f"{key} must be nonnegative")
        if item.get("atmosphere") is not None and not isinstance(item["atmosphere"], str):
            raise ValueError("atmosphere must be a string or null")
    if any(type(stage.get("stage_index")) is not int or stage["stage_index"] < 0 for stage in stages):
        raise ValueError("stage_index must be a nonnegative integer")
    values["stages"] = [ConditionProcessStage(**stage) for stage in stages]
    if len({stage.stage_index for stage in values["stages"]}) != len(stages):
        raise ValueError("Process stage indices must be distinct")
    recipe = asdict(build_resolved_recipe_from_inputs([ConditionComponentInput(**item) for item in components], **values))
    from .route_investigation import _evidence

    provenance = _evidence(operations, evidence_refs)
    return {"schema_version": "proposed_recipe_check.v1", "reaction_smiles": reaction_smiles,
            "proposed_recipe": recipe, "compatibility": asdict(operations.assess_recipe(reaction_smiles, recipe)),
            "origin": "agent_proposal", "experimental_feasibility": "not_established", **provenance,
            "process_coverage": "stages_recorded_not_evaluated" if stages else "no_stages_supplied",
            "limitations": ["Checks apply to this saved recipe, not arbitrary prose conditions or changed stages.",
                            "Normalization and absence of known conflicts do not establish conversion or yield.",
                            "Stage-specific mixtures, addition order and workup chemistry are not evaluated by these rules."]}
