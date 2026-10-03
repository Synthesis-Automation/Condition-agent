"""Deterministic fragment-witness projection and single-step transfer evaluation.

Discovery evidence proposes sites; complete source admission and forward graph
validation independently determine which proposals can be returned.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from rdkit import Chem

from .generic_library import GenericTemplateLibrary

DEFAULT_POLICY = Path(__file__).with_name("definitions") / "fragment_guided_retro_policy.v1.json"


def load_transfer_policy(path: str | Path | None = None) -> dict[str, Any]:
    """Validate bounded alignment-projection work for explicit query transfer."""
    path = Path(path) if path is not None else Path(__file__).with_name("definitions") / "fragment_transfer.v1.json"
    value = json.loads(path.read_text("utf-8"))
    if (not isinstance(value, dict)
            or set(value) != {"schema_version", "definition_version", "max_target_alignments",
                      "max_projection_witnesses"}
            or value["schema_version"] != "fragment_transfer_policy.v1"
            or value["definition_version"] != "fragment_transfer.v1@1.0"):
        raise ValueError("Unsupported fragment transfer policy")
    for key, maximum in (("max_target_alignments", 64), ("max_projection_witnesses", 2048)):
        if type(value[key]) is not int or not 1 <= value[key] <= maximum:
            raise ValueError(f"Invalid transfer policy field: {key}")
    return value


def restore_source_smiles(value: Any) -> str:
    """Restore complete index text chunks, checking offsets, length and source hash."""
    if isinstance(value, str) and value:
        return value
    if not isinstance(value, dict) or value.get("truncated") is not False:
        raise ValueError("Source reaction text is missing or truncated")
    parts, offset = [], 0
    for chunk in value.get("chunks", []):
        text = chunk["text"]
        if chunk["start"] != offset or chunk["end"] != offset + len(text):
            raise ValueError("Source reaction text chunk offsets conflict")
        parts.append(text)
        offset += len(text)
    text = "".join(parts)
    if (not text or len(text) != value.get("total_characters")
            or hashlib.sha256(text.encode()).hexdigest() != value.get("source_sha256")):
        raise ValueError("Source reaction text length or hash conflicts")
    return text


def load_policy(path: str | Path) -> dict[str, Any]:
    """Validate the bounded experiment policy; no executable names are loaded."""
    value = json.loads(Path(path).read_text("utf-8"))
    required = {
        "schema_version", "definition_version", "query_selection", "query_limit",
        "hit_limit", "search_timeout_seconds", "require_complete_search",
        "require_resolved_embedding", "max_focus_bonds", "max_source_observations",
        "levels", "top_k", "max_templates_to_apply", "max_candidates_to_validate",
    }
    if not isinstance(value, dict) or set(value) != required:
        raise ValueError("Unsupported POC policy fields")
    if (value["schema_version"] != "fragment_guided_retro_poc_policy.v1"
            or value["definition_version"] != "fragment_guided_retro_poc.v1@1.0"
            or value["query_selection"] != "suggestion_order"):
        raise ValueError("Unsupported POC policy identity or selection")
    for key, maximum in {
        "query_limit": 5, "hit_limit": 10, "search_timeout_seconds": 30,
        "max_focus_bonds": 5, "max_source_observations": 10, "top_k": 10,
        "max_templates_to_apply": 100, "max_candidates_to_validate": 25,
    }.items():
        if type(value[key]) is not int or not 1 <= value[key] <= maximum:
            raise ValueError(f"Invalid bounded policy field: {key}")
    if (value["require_complete_search"] is not True
            or value["require_resolved_embedding"] is not True
            or value["levels"] != ["L2", "L1", "L0"]):
        raise ValueError("This POC requires complete searches and resolved embeddings")
    return value


def project_construction_guidance(
    suggestions: dict[str, Any], searches: list[dict[str, Any]],
    policy: dict[str, Any],
) -> dict[str, Any]:
    """Project observed formed edges through query alignments, not reaction maps.

    Target bonds remain proposed construction sites. Boundary edits, retention,
    unresolved/truncated embeddings and incomplete searches cannot seed this POC.
    The complete source reaction is separately subjected to compiler admission.
    """
    target = Chem.MolFromSmiles(suggestions["target_smiles"])
    candidates = {c["candidate_id"]: c for c in suggestions["candidates"]}
    bonds: dict[tuple[int, int], list[dict[str, Any]]] = {}
    sources: dict[str, dict[str, Any]] = {}
    exclusions: list[dict[str, Any]] = []
    for search in searches:
        candidate = candidates[search["candidate_id"]]
        result = search.get("result") or {}
        validation = result.get("target_validation") or {}
        if (search.get("execution_status") != "completed"
                or result.get("search_status") != "complete"
                or validation.get("matches_target") is not True
                or validation.get("target_smiles") != suggestions["target_smiles"]):
            exclusions.append({"candidate_id": candidate["candidate_id"],
                               "reason": "incomplete_or_unvalidated_search"})
            continue
        for hit in result.get("hits", []):
            for match_index, match in enumerate(hit["matches"]):
                if (match.get("evidence_status") != "validated_supplied_mapping"
                        or "unresolved" in match["relationships"]
                        or match.get("embedding_truncated")):
                    exclusions.append({"observation_id": hit["observation_id"],
                                       "reason": "unresolved_or_truncated_embedding"})
                    continue
                original_ids = match["query_to_original_product_atoms"]
                target_ids = candidate["query_atom_target_ids"]
                if (len(original_ids) != len(target_ids)
                        or len(set(original_ids)) != len(original_ids)):
                    raise ValueError("Fragment alignment is not a one-to-one query embedding")
                alignment = dict(zip(original_ids, target_ids))
                for witness in match["witnesses"]:
                    if (witness["kind"] != "formed"
                            or witness["relationship"] != "constructed"):
                        continue
                    endpoints = witness["product_atoms"]
                    if len(endpoints) != 2 or any(a not in alignment for a in endpoints):
                        raise ValueError("Construction witness lies outside its query embedding")
                    pair = tuple(sorted(alignment[a] for a in endpoints))
                    if target.GetBondBetweenAtoms(*pair) is None:
                        raise ValueError("Projected construction site is not a target bond")
                    raw_reaction = hit["record"].get("reaction_smiles")
                    if not raw_reaction or (
                        isinstance(raw_reaction, dict) and raw_reaction.get("truncated")
                    ):
                        exclusions.append({"observation_id": hit["observation_id"],
                                           "reason": "missing_or_truncated_source_structure"})
                        continue
                    source_smiles = restore_source_smiles(raw_reaction)
                    support = {
                        "observation_id": hit["observation_id"],
                        "reaction_id": hit["reaction_id"],
                        "reference_id": hit["reference_id"],
                        "candidate_id": candidate["candidate_id"],
                        "search_ref": search["artifact_ref"], "match_index": match_index,
                        "witness": witness,
                        "source_reaction_smiles": source_smiles,
                    }
                    bonds.setdefault(pair, []).append(support)
                    sources.setdefault(hit["observation_id"], {
                        **hit["record"], "reaction_smiles": source_smiles,
                    })
    ranked = sorted(bonds, key=lambda pair: (
        -len({s["reference_id"] or s["observation_id"] for s in bonds[pair]}), pair,
    ))
    return {
        "target_smiles": suggestions["target_smiles"],
        "projection_basis": "query_embedding_alignment_not_reaction_mapping",
        "focus_bonds": [{"target_atom_ids": list(pair), "supports": bonds[pair]}
                        for pair in ranked[:policy["max_focus_bonds"]]],
        "eligible_bond_count": len(ranked),
        "source_records": [sources[key] for key in sorted(sources)
                           [:policy["max_source_observations"]]],
        "eligible_source_count": len(sources), "exclusions": exclusions,
        "scope": "returned_hit_sample_not_all_indexed_construction_records",
    }


def scientific_projection(value: dict[str, Any]) -> dict[str, Any]:
    """Remove execution measurements from a completed fragment search for comparison."""
    return {key: item for key, item in value.items() if key != "execution"}


def project_selected_query_guidance(
    target_smiles: str, result: dict[str, Any], selected_observation_ids: list[str],
    policy: dict[str, Any],
) -> dict[str, Any]:
    """Project selected construction witnesses through every bounded target alignment.

    Inputs are canonical search results, not user-provided source reactions.
    Target-query symmetries remain alternative hypotheses; truncated target
    enumeration cannot seed guidance. Query edits never weaken source admission.
    """
    from reactive_taxonomy.fragment_search import (
        compile_fragment_query, fragment_embeddings,
        indexed_product, validate_fragment_target,
    )

    description = result["query"]
    query = compile_fragment_query(
        description["expression"], description["query_format"], description["topology"],
    )
    target, _ = indexed_product(target_smiles)
    validation = validate_fragment_target(query, target)
    supplied = result.get("target_validation") or {}
    if (not validation.matches_target or supplied.get("matches_target") is not True
            or supplied.get("target_smiles") != target
            or supplied.get("query_id") != query.query_id
            or description.get("query_id") != query.query_id):
        raise ValueError("Selected search evidence does not validate this target and query")
    available_ids = {hit["observation_id"] for hit in result["hits"]}
    if not selected_observation_ids or len(set(selected_observation_ids)) != len(selected_observation_ids):
        raise ValueError("Select unique source observation IDs")
    if any(identity not in available_ids for identity in selected_observation_ids):
        raise ValueError("Selected source is absent from the freshly returned search sample; search again")
    transfer_policy = load_transfer_policy()
    alignments, truncated = fragment_embeddings(
        query, Chem.MolFromSmiles(target),
        maximum=transfer_policy["max_target_alignments"],
    )
    candidates = [{
        "candidate_id": f"{query.query_id}:target:{number}",
        "kind": "chosen_query", "query": query.expression,
        "query_format": query.query_format, "topology": query.topology,
        "query_atom_target_ids": list(alignment), "target_atom_ids": sorted(alignment),
        "reasons": ["User-chosen query matched to the complete canonical target."],
        "cautions": ["TARGET_ALIGNMENTS_ARE_ANALOGUE_HYPOTHESES"],
    } for number, alignment in enumerate(alignments)]
    suggestions = {"target_smiles": target, "candidates": candidates}
    selected_result = {**result, "hits": [
        hit for hit in result["hits"] if hit["observation_id"] in selected_observation_ids
    ]}
    witness_count = len(alignments) * sum(
        witness["kind"] == "formed" and witness["relationship"] == "constructed"
        for hit in selected_result["hits"] for match in hit["matches"]
        for witness in match["witnesses"]
    )
    projection_limited = witness_count > transfer_policy["max_projection_witnesses"]
    searches = [{
        "candidate_id": candidate["candidate_id"], "artifact_ref": f"query:{query.query_id}",
        "execution_status": "completed", "result": selected_result,
    } for candidate in candidates] if not truncated and not projection_limited else []
    guidance = project_construction_guidance(suggestions, searches, policy)
    if truncated:
        guidance["exclusions"].append({"reason": "target_alignment_enumeration_truncated"})
    if projection_limited:
        guidance["exclusions"].append({"reason": "projection_witness_budget_exceeded"})
    guidance["target_alignment_count"] = len(alignments)
    guidance["target_alignments_truncated"] = truncated
    guidance["selected_observation_ids"] = list(selected_observation_ids)
    guidance["transfer_policy"] = transfer_policy
    guidance["estimated_projection_witness_count"] = witness_count
    return {"suggestions": suggestions, "guidance": guidance}


def evaluate_transfers(
    parameters: dict[str, Any], *, library: GenericTemplateLibrary | None = None,
    repeat: bool = True, source_library_path: str | Path | None = "selected_source_operators.json",
) -> dict[str, Any]:
    """Compare baseline, admitted source transfers and witness-directed library search."""
    from core_retrosynthesis import build_generic_library, load_generic_library
    from core_retrosynthesis.generic_search import disconnect_operator_ladder_detailed
    from core_retrosynthesis.generic_library import save_generic_library

    guidance, policy = parameters["guidance"], parameters["policy"]
    target = guidance["target_smiles"]
    library = library if library is not None else load_generic_library(parameters["library"])
    admissions: list[dict[str, Any]] = []
    source_library = build_generic_library(
        guidance["source_records"], engine="reaction_core", admission_mode="data_driven",
        levels=tuple(policy["levels"]), admission_callback=admissions.append,
    )
    # This is a bounded experimental library, never a replacement corpus/index.
    if source_library_path is not None:
        save_generic_library(source_library, Path(source_library_path))
    options = {key: policy[key] for key in (
        "top_k", "max_templates_to_apply", "max_candidates_to_validate",
    )}
    options.update(use_hierarchical_ranking=False, diversify=False)

    def search(selected_library: Any, bond: list[int] | None) -> dict[str, Any]:
        if not selected_library.templates:
            return {"candidates": [], "diagnostics": None,
                    "status": "no_admitted_source_operators"}
        focus = ({"required_disconnection_bond": bond, "focus_target_smiles": target}
                 if bond is not None else {})
        candidates, diagnostics = disconnect_operator_ladder_detailed(
            target, selected_library, **options, **focus,
        )
        serialized = [c.to_dict() for c in candidates]
        for candidate in serialized:
            if candidate["forward_validation_status"] != "verified_signature":
                raise ValueError("POC results require a verified reaction signature")
            if bond is not None and (candidate.get("bond_focus_check") or {}).get("status") != "verified":
                raise ValueError("Focused result lacks final observed-formation validation")
        return {"status": "completed", "candidates": serialized,
                "diagnostics": diagnostics.to_dict()}

    def run_searches() -> dict[str, Any]:
        return {
            "baseline": search(library, None) if parameters.get("include_baseline", True) else {
                "status": "not_requested", "candidates": [], "diagnostics": None,
            },
            "guided": [{"target_atom_ids": b["target_atom_ids"],
                        "direct_source_transfer": search(source_library, b["target_atom_ids"]),
                        "witness_directed_library": search(library, b["target_atom_ids"])}
                       for b in guidance["focus_bonds"]],
        }

    first = run_searches()
    second = run_searches() if repeat else None
    return {
        **first, "comparison": summarize_transfers(first), "source_admissions": admissions,
        "compiled_source_template_count": len(source_library.templates),
        "source_operator_definition": source_library.definition,
        "prepared_library_template_count": len(library.templates),
        "repeat_scientific_results_identical": first == second if repeat else None,
        "budget_scope": "per_specificity_level_per_bond_per_arm; descriptive comparison, unequal total work",
        "source_library_file": str(source_library_path) if source_library_path is not None else None,
    }


def summarize_transfers(transfers: dict[str, Any]) -> dict[str, Any]:
    """Count distinct returned precursor sets and actual work across search arms."""
    baseline = transfers["baseline"]
    guided = [branch[arm] for branch in transfers["guided"] for arm in (
        "direct_source_transfer", "witness_directed_library",
    )]
    baseline_sets = {c["precursor_smiles"] for c in baseline["candidates"]}
    guided_sets = {c["precursor_smiles"] for arm in guided for c in arm["candidates"]}

    def work(arms: list[dict[str, Any]]) -> dict[str, int]:
        diagnostics = [arm.get("diagnostics") or {} for arm in arms]
        return {
            "template_applications": sum(
                level.get("applied_template_count", 0)
                for item in diagnostics for level in item.get("level_diagnostics", {}).values()
            ),
            "validation_attempts": sum(item.get("validation_attempt_count", 0) for item in diagnostics),
        }

    return {
        "baseline_unique_precursor_count": len(baseline_sets),
        "guided_unique_precursor_count": len(guided_sets),
        "baseline_requested": baseline.get("status") != "not_requested",
        "additional_guided_precursor_sets": (
            sorted(guided_sets - baseline_sets) if baseline.get("status") != "not_requested" else []
        ),
        "baseline_work": work([baseline]), "guided_work": work(guided),
        "interpretation": "descriptive_coverage_under_unequal_total_work_not_accuracy",
    }


