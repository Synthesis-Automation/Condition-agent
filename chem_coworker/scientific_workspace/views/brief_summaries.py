"""Brief decision views projected directly from complete saved results."""

from __future__ import annotations

from collections.abc import Mapping
from itertools import islice
import json
from typing import Any


def summarize_call_brief(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Project saved results directly; never count or select from an earlier preview."""
    operation = payload.get("operation")
    source = payload.get("result")
    source = source if isinstance(source, Mapping) else {}
    result = {"summary_schema_version": "scientific_call_brief.v2"}
    result.update(_brief_fields(payload, (
        "operation", "execution_status", "error", "duration_seconds", "result_bytes",
    )))
    if payload.get("schema_version") in {"literature_source.v1", "literature_excerpt.v1"}:
        result["result_summary"] = _brief_fields(payload, (
            "source_url", "title", "acquisition", "retrieval_status", "text_characters",
            "text_sha256", "reported_locator", "reported_reference_id", "capture_file",
            "error", "recovery", "claim_support", "review_status", "network_request_skipped", "blocked_by_ref",
        ))
        result["result_summary"]["extraction"] = _brief_fields(payload.get("extraction", {}), (
            "status", "method", "limitations", "error",
        ))
        result["inspection"] = {"projection_only": True, "detail_available": True,
                                "hint": "Use inspect_source for exact bounded source passages."}
        return result
    if not isinstance(payload.get("result"), Mapping):
        result["result_summary"] = {}
        result["inspection"] = {"result_path": "$.result", "projection_only": True,
                                "detail_available": True, "result_present": "result" in payload}
        return result
    overview = _brief_fields(source, (
        "status", "valid", "error", "warnings", "limitations", "review_status",
        "experimental_feasibility", "ranking", "search_status", "stop_reason",
        "canonical_smiles", "canonical_target_smiles", "retrieval_level",
        "candidate_count", "compatible_candidate_count", "returned_count",
        "definition_version", "selection_mode", "generated_count", "rejected_count",
        "query", "target_smiles", "source_ref", "step_id",
        "source_coverage_complete", "group_count_scope", "ranking_scope",
    ))
    if "warnings" in source and source["warnings"] in ([], ()):
        # Keep the existing root-level warning contract for workspace clients.
        overview["warnings"] = []
    if operation == "inspect_step_precedents":
        overview.update(_brief_fields(source, ("scope", "saved_match_count", "available_template_records",
                                              "retrieval_truncated", "page", "distinct_references_on_page",
                                              "observation_page", "reference_catalog_status", "procedure_catalog_status",
                                              "assessment_status", "assessment_warnings")))
        overview["selection"] = _brief_fields(source.get("selection", {}), (
            "step_id", "realization_id", "target_smiles", "precursor_smiles",
        ))
        overview["precedents"] = [_brief_fields(item, (
            "match_id", "reaction_id", "reference_id", "reaction_smiles", "support_kind",
            "product_similarity", "precursor_similarity", "same_recorded_product", "same_recorded_precursors",
            "limitations",
        )) for item in source.get("precedents", [])[:3]]
        for brief, item in zip(overview["precedents"], source.get("precedents", [])):
            brief["product_comparison"] = _brief_fields(item.get("product_comparison", {}), (
                "status", "same_constitution", "stereo_relationship", "left_coverage", "right_coverage",
                "alignment_ambiguous", "search_timed_out", "warnings",
            ))
            brief["observation_count"] = len(item.get("observations", []))
            brief["procedure_count"] = len(item.get("procedures", []))
    elif operation == "analyze_molecule":
        overview.update(_brief_fields(source, ("formal_charge", "stereocenters")))
        for key, fields in (("motifs", ("motif_id", "chemist_label")),
                            ("reactive_sites", ("chemist_label", "availability", "warnings"))):
            overview[f"{key}_count"] = len(source.get(key, []))
            overview[key] = [_brief_fields(item, fields) for item in source.get(key, [])[:3]]
    elif operation == "assess_starting_material":
        overview.update(_brief_fields(source, (
            "stop_expansion", "availability", "molecular_weight", "mw_threshold",
            "unavailable_starting_materials",
        )))
        for key in ("registry", "literature"):
            overview[key] = _brief_fields(source.get(key, {}), (
                "status", "reason", "candidates", "matched_identifier", "index_id",
                "source_scope", "source_coverage_complete", "precedents", "has_more",
                "error_type", "message", "warnings",
            ))
    elif operation == "compare_molecules":
        overview.update(_brief_fields(source, (
            "same_constitution", "stereo_relationship", "core_method", "core_atom_count",
            "left_coverage", "right_coverage", "alignment_count_observed", "alignment_ambiguous",
            "embeddings_truncated", "alignments_truncated", "search_timed_out",
        )))
        for side in ("left", "right"):
            overview[side] = _brief_fields(source.get(side, {}), ("canonical_smiles", "warnings"))
    elif operation == "inspect_reactive_sites":
        overview.update(_brief_fields(source, ("selected_atom_ids", "radius")))
        overview["sites_count"] = len(source.get("sites", []))
        overview["other_sites_count"] = len(source.get("other_sites", []))
        overview["molecule"] = _brief_fields(source.get("molecule", {}), ("canonical_smiles", "warnings"))
        overview["sites"] = [
            _brief_fields(item, ("hypothesis_id", "chemist_label", "atom_indices", "availability", "warnings"))
            for item in source.get("sites", [])[:3] if isinstance(item, Mapping)
        ]
        overview["other_sites"] = [
            _brief_fields(item, ("hypothesis_id", "chemist_label", "atom_indices", "warnings"))
            for item in source.get("other_sites", [])[:3] if isinstance(item, Mapping)
        ]
        overview["site_profiles"] = []
        for item in source.get("environments", [])[:2]:
            profile = item.get("reactivity_profile", {})
            compact = _brief_fields(item, ("hypothesis_id", "center_atom_index"))
            compact.update(_brief_fields(profile, ("status", "context_kind", "flags")))
            for key, field in (("steric", "accessibility_class"), ("electronic", "activation_class")):
                descriptor = profile.get(key, {})
                compact[key] = _brief_fields(descriptor, (field,))
                if isinstance(descriptor.get("evidence"), Mapping):
                    compact[key]["evidence"] = _brief_fields(descriptor["evidence"], (
                        "source", "method", "confidence", "warnings",
                    ))
            overview["site_profiles"].append(compact)
    elif operation == "analyze_reaction":
        overview.update(_brief_fields(source, ("evidence_quality", "transformation_class", "named_family")))
        for key in ("observation", "interpretation", "reaction_completeness"):
            if isinstance(source.get(key), Mapping):
                overview[key] = _brief_fields(source[key], ("status", "valid", "evidence_quality", "warnings"))
        if isinstance(source.get("reaction_signature"), Mapping):
            signature = source["reaction_signature"]
            overview["bond_changes"] = _brief_fields(signature, ("formed_bond_types", "broken_bond_types"))
            overview["bond_changes"].update({
                f"{key}_count": len(signature.get(key, []))
                for key in ("order_changes", "hydrogen_changes", "stereo_changes")
                if key in signature
            })
    elif operation == "disconnect_composite":
        overview.update(_brief_fields(source, ("catalog_id", "diagnostics")))
        overview["actions"] = []
        for item in source.get("actions", [])[:3]:
            if not isinstance(item, Mapping):
                continue
            action = _brief_fields(item, (
                "action_id", "strategy_id", "intermediate_smiles", "terminal_precursor_smiles",
                "physical_step_count", "physical_step_cost", "condition_compatibility_status", "one_pot_status",
            ))
            action["dependency"] = _brief_fields(item.get("dependency", {}), (
                "admitted", "status", "relationship_class", "dependency_class", "lineage_status", "warnings",
            ))
            action["physical_steps"] = [_brief_fields(step, (
                "forward_step_number", "reaction_smiles", "operator_id", "precedent_reaction_ids",
                "forward_validation_status", "condition_status", "evidence_kind",
                "precursor_compatibility_disposition", "reaction_compatibility_disposition",
            )) for step in item.get("physical_steps", [])]
            overview["actions"].append(action)
        overview["dependency_reviews"] = source.get("dependency_reviews", [])[:3]
    elif operation == "disconnect_target":
        overview.update(_brief_fields(source, ("bond_focus",)))
        if source.get("search_diagnostics"):
            diagnostics = source["search_diagnostics"]
            overview["search_diagnostics"] = _brief_fields(diagnostics, (
                "levels_attempted", "proposed_action_count", "validation_attempt_count", "valid_action_count",
            ))
            overview["search_diagnostics"]["level_diagnostics"] = {
                level: _brief_fields(item, (
                    "applied_template_count", "validation_attempt_count", "valid_candidate_count",
                    "focus_matched_count", "focus_rejected_count", "focus_unresolved_count",
                    "focus_validation_rejected_count", "template_budget_excluded_count",
                    "validation_budget_excluded_count",
                )) for level, item in diagnostics.get("level_diagnostics", {}).items()
            }
        strategies = source.get("strategies", [])
        overview["strategy_count"] = len(strategies)
        overview["strategies"] = []
        for strategy in strategies[:3]:
            item = _brief_fields(strategy, ("strategy_id", "strategy_rank", "independent_reference_support"))
            if isinstance(strategy.get("representative"), Mapping):
                item["representative"] = _brief_fields(strategy["representative"], (
                    "realization_id", "precursor_smiles", "forward_validation_status",
                    "precursor_compatibility_disposition", "reaction_compatibility_disposition",
                    "template_id", "operator_id", "precedent_reaction_ids",
                    "bond_focus_check",
                ))
                warnings = strategy["representative"].get("selectivity_warnings", [])
                if warnings:
                    item["representative"]["selectivity_warnings"] = [
                        _brief_fields(warning, ("code", "message", "conditions_evaluated"))
                        for warning in warnings[:2] if isinstance(warning, Mapping)
                    ]
                    item["representative"]["selectivity_warning_count"] = len(warnings)
            overview["strategies"].append(item)
    elif operation == "assess_retro_validity":
        overview["validity"] = _brief_fields(source.get("validity", {}), (
            "status", "structural_status", "precedent_grade", "evidence_rank", "suggested_action",
            "forward_status", "forward_execution_status", "cautions", "unresolved_checks", "warnings",
            "score_semantics", "ranking_influence",
        ))
        for name in ("operator_precedent_support", "corpus_precedent_support"):
            support = source.get("validity", {}).get(name)
            if isinstance(support, Mapping):
                overview["validity"][name] = _brief_fields(support, (
                    "status", "strongest_level", "level_counts", "distinct_reference_count",
                    "candidate_count", "qualified_count", "candidate_truncated", "matches_truncated", "scope",
                ))
        overview["artifact_warnings"] = source.get("artifact_warnings", [])
    elif operation in {"assess_route_step", "inspect_route_step", "assess_route_proposal", "revise_route_branch"}:
        assessment = source.get("assessment")
        if isinstance(assessment, Mapping):
            overview["assessment"] = _brief_fields(assessment, (
                "status", "actionable", "admission_eligible", "strongest_evidence_tier",
                "warnings", "unresolved_step_ids", "disconnected_step_ids",
            ))
            gate_key = "topology_gates" if operation in {"assess_route_proposal", "revise_route_branch"} else "gates"
            gates = assessment.get(gate_key, [])
            non_pass_gates = [
                gate for gate in gates
                if isinstance(gate, Mapping) and gate.get("status") != "pass"
            ]
            non_pass_gates.sort(key=lambda gate: 0 if gate.get("status") in {
                "failed", "invalid", "conflicting", "unresolved", "unknown",
            } else 1)
            overview["assessment"]["non_pass_gates"] = [
                _brief_fields(gate, ("gate_id", "status", "summary", "warnings"))
                for gate in non_pass_gates[:4]
            ]
            overview["assessment"]["non_pass_gate_count"] = len(non_pass_gates)
            if "step_assessments" in assessment:
                steps = assessment["step_assessments"]
                overview["assessment"]["step_count"] = len(steps)
                if len(steps) > 10:
                    steps = sorted(steps, key=lambda step: (
                        (step.get("assessment") or {}).get("actionable") is not False
                        and not (step.get("assessment") or {}).get("warnings")
                    ))
                    overview["assessment"]["steps_preview_order"] = "non_actionable_or_warned_first"
                overview["assessment"]["steps"] = [
                    {"step_id": step.get("external_step_id"),
                     **(_brief_fields(step["assessment"], ("status", "actionable", "warnings"))
                        if isinstance(step.get("assessment"), Mapping) else {})}
                    for step in steps[:10] if isinstance(step, Mapping)
                ]
        if isinstance(source.get("material_constraints"), Mapping):
            overview["material_constraints"] = _brief_fields(source["material_constraints"], (
                "status", "stock_availability", "unavailable_starting_materials",
            ))
    elif operation == "generate_weak_label_screening_array":
        overview.update(_brief_fields(source, (
            "recommendation_mode", "reaction_type_id", "source_reaction_type_candidates",
            "recipe_count", "excluded_candidate_count",
        )))
        recommendations = source.get("recommendations", [])
        overview["recommendation_count"] = len(recommendations)
        overview["recommendations"] = []
        for item in recommendations[:3]:
            if not isinstance(item, Mapping):
                continue
            compact = _brief_fields(item, (
                "rank", "recipe_id", "support", "source_reaction_types", "source_row_numbers",
                "compatibility_score", "cautions",
            ))
            if isinstance(item.get("resolved_recipe"), Mapping):
                compact["recipe"] = _brief_recipe(item["resolved_recipe"])
            overview["recommendations"].append(compact)
    elif operation == "recommend_conditions":
        recommendations = source.get("recommendations", [])
        overview["recommendation_count"] = len(recommendations)
        overview["recommendations"] = []
        for item in recommendations[:3]:
            if not isinstance(item, Mapping):
                continue
            compact = _brief_fields(item, (
                "rank", "recipe_id", "match_label", "retrieval_level", "compatibility_status",
                "reference_support", "observation_support", "precedent_reaction_ids", "cautions",
            ))
            if isinstance(item.get("resolved_recipe"), Mapping):
                compact["recipe"] = _brief_recipe(item["resolved_recipe"])
            overview["recommendations"].append(compact)
    elif operation == "resolve_recipe":
        overview.update(_brief_fields(source, ("recipe_id",)))
        overview["recipe"] = _brief_recipe(source)
    elif operation in {"get_precedents", "get_procedures"}:
        overview.update(_brief_fields(source, (
            "total", "offset", "next_offset", "missing_reaction_ids", "record_scope", "index_scope",
        )))
        records = source.get("records", [])
        overview["saved_record_count"] = len(records)
        overview["records"] = []
        for index, item in enumerate(records[:3]):
            compact = _brief_fields(item, (
                "reaction_id", "observation_id", "reference_id", "source_url", "source_uri",
                "chemistry_status", "condition_status", "condition_uncertain", "warnings",
            ))
            compact["saved_index"] = index
            if operation == "get_procedures":
                text = item.get("procedure_text")
                compact["procedure_characters"] = len(text) if isinstance(text, str) else 0
                compact["procedure_text_present"] = isinstance(text, str) and bool(text.strip())
            else:
                compact.update(_brief_fields(item, ("reaction_smiles", "recipe_id", "yield_pct")))
            overview["records"].append(compact)
    elif operation == "assess_recipe":
        overview.update(_brief_fields(source, ("compatible", "hard_conflicts", "unresolved_requirements")))
        if isinstance(source.get("coverage"), Mapping):
            overview["coverage"] = _brief_fields(source["coverage"], (
                "schema_version", "assessment_status", "capability_status",
                "condition_identity_status", "unresolved_components", "limitations", "score_meaning",
            ))
    elif operation == "assess_route_step_forward":
        if isinstance(source.get("assessment"), Mapping):
            overview["assessment"] = _brief_fields(source["assessment"], (
                "targeted_replay_status", "intended_match", "disposition", "validity",
                "advisory_only", "warnings",
            ))
        if isinstance(source.get("execution"), Mapping):
            execution = source["execution"]
            overview["execution"] = _brief_fields(execution, ("timings",))
            stages = execution.get("stages", [])
            overview["execution"]["stages"] = [
                _brief_fields(item, ("stage", "status", "duration_seconds")) for item in stages[-2:]
            ]
            overview["execution"]["stage_count"] = len(stages)
    elif operation == "compare_route_proposals":
        alternatives = source.get("alternatives", [])
        overview["alternative_count"] = len(alternatives)
        overview["alternatives"] = [
            _brief_fields(item, ("status", "source_ref", "route_id", "step_count",
                                 "unresolved_step_ids", "limitations"))
            for item in alternatives[:3] if isinstance(item, Mapping)
        ]
    elif operation == "propose_condition_adaptation":
        overview.update(_brief_fields(source, ("transfer_status", "assumptions", "risks")))
        if isinstance(source.get("compatibility"), Mapping):
            overview["compatibility"] = _brief_fields(source["compatibility"], (
                "status", "compatible", "hard_conflicts", "unresolved_requirements",
            ))
        overview["changes"] = [
            _brief_fields(item, ("field", "basis", "reason"))
            for item in source.get("changes", [])[:3] if isinstance(item, Mapping)
        ]
    elif operation == "search_fragment_precedents":
        overview.update(_brief_fields(source, ("counts", "refinement_hints", "output_truncated")))
        if isinstance(source.get("target_validation"), Mapping):
            overview["target_validation"] = _brief_fields(source["target_validation"], (
                "matches_target", "target_smiles", "query_id", "schema_version",
                "definition_version", "compiler_version",
            ))
        if isinstance(source.get("query"), Mapping):
            overview["query"] = _brief_fields(source["query"], ("expression", "query_format", "topology", "query_id"))
        hits = source.get("hits", [])
        overview["saved_hit_count"] = len(hits)
        overview["hits"] = []
        for index, item in enumerate(hits[:3]):
            if not isinstance(item, Mapping):
                continue
            compact = _brief_fields(item, ("hit_id", "reaction_id", "observation_id", "reference_id",
                                          "product_smiles", "relationship_summary",
                                          "procedure_availability", "warnings"))
            record = item.get("record", {})
            compact.update(_brief_fields(record, ("reaction_smiles", "yield_pct", "temperature_c", "time_h")))
            compact["publication"] = _publication(record)
            compact["record_path"] = ["result", "hits", index, "record"]
            overview["hits"].append(compact)
        overview["hit_page"] = {"shown": len(overview["hits"]), "total": len(hits),
                                "next_offset": 3 if len(hits) > 3 else None,
                                "path": ["result", "hits"], "count_scope": "saved_result"}
    elif operation == "prepare_literature_reaction":
        overview.update(_brief_fields(source, ("structure_origin", "participant_count", "source_provenance",
                                               "source_conflicts", "structure_checks", "limitations")))
        overview["hint"] = "Use w.prepared_literature_reaction(ref) in the answer; do not print the full block."
    elif operation == "suggest_search_fragments":
        overview.update(_brief_fields(source, ("generation_truncated", "output_truncated")))
        overview["candidate_count"] = len(source.get("candidates", []))
        overview["candidates"] = [
            _brief_fields(item, ("candidate_id", "kind", "query", "reasons", "cautions"))
            for item in source.get("candidates", [])[:3] if isinstance(item, Mapping)
        ]
    elif operation == "inspect_condition_precedents":
        overview.update(_brief_fields(source, (
            "observation_count", "distinct_reference_count", "page", "missing_reference_observation_ids",
        )))
        if isinstance(source.get("procedure_catalog"), Mapping):
            overview["procedure_catalog"] = _brief_fields(source["procedure_catalog"], ("availability",))
        overview["precedents"] = []
        for item in source.get("precedents", [])[:2]:
            if not isinstance(item, Mapping):
                continue
            precedent = _brief_fields(item, ("missing_operating_fields", "procedure_link_scope"))
            for key, fields in (
                ("observation", ("observation_id", "reaction_id", "reference_id", "chemistry_status", "condition_status")),
                ("compatibility", ("status", "compatible", "hard_conflicts")),
            ):
                if isinstance(item.get(key), Mapping):
                    precedent[key] = _brief_fields(item[key], fields)
            observation = item.get("observation", {})
            if isinstance(observation.get("resolved_recipe"), Mapping):
                precedent["recipe"] = _brief_recipe(observation["resolved_recipe"])
            overview["precedents"].append(precedent)
    else:
        overview.update(_brief_fields(source, (
            "compatible", "hard_conflicts", "unresolved_requirements", "total", "next_offset",
            "observation_count", "distinct_reference_count", "transfer_status", "question",
        )))
    result["result_summary"] = overview
    result["inspection"] = {
        "result_path": "$.result", "projection_only": True,
        "detail_available": True,
        "count_scope": "saved_result",
    }
    if not overview and source:
        result["inspection"]["result_fields"] = [str(key)[:80] for key in islice(source, 20)]
    return result


def _publication(record: Mapping[str, Any]) -> dict[str, Any]:
    """Project explicit source bibliography only; never infer it from reaction names."""
    identity = record.get("reference_identity", {})
    identity = identity if isinstance(identity, Mapping) else {}
    result = _brief_fields(identity, ("reference_id", "doi", "url", "title", "authors", "citation", "year"))
    raw = record
    for key in ("source", "raw_fields", "source_observation", "raw_fields"):
        raw = raw.get(key, {}) if isinstance(raw, Mapping) else {}
    if isinstance(raw, Mapping):
        for key in ("authors", "citation", "doi", "title", "url"):
            if key not in result and raw.get(key):
                result.update(_brief_fields(raw, (key,)))
    return result


def bound_call_summary(summary: dict[str, Any], payload: Mapping[str, Any]) -> dict[str, Any]:
    """Cap console JSON at 16 KiB and explicitly require detail when a view is too large."""
    def size(value: dict[str, Any]) -> int:
        return len(json.dumps(value, ensure_ascii=False, separators=(",", ":")).encode("utf-8"))

    if size(summary) <= 16000:
        return summary
    source = payload.get("result", {})
    source = source if isinstance(source, Mapping) else {}
    result = {key: value for key, value in summary.items() if key not in {"result_summary", "inspection"}}
    event = dict(result.get("event", {}))
    refs = event.get("evidence_refs", [])
    event["evidence_refs"] = refs[:5]
    event["evidence_ref_count"] = len(refs)
    event["evidence_refs_omitted"] = max(0, len(refs) - 5)
    result["event"] = event
    result["result_summary"] = _brief_fields(source, (
        "status", "valid", "error", "compatible", "hard_conflicts", "warnings", "limitations",
        "execution_status", "experimental_feasibility", "coverage",
    ))
    result["inspection"] = {
        "projection_only": True, "detail_available": True, "summary_omitted": True,
        "serialized_byte_limit": 16000, "reason": "console_summary_budget",
        "hint": "Large result omitted. Inspect decision-critical gates, warnings and saved structures before relying on it.",
    }
    if size(result) > 16000:
        result["result_summary"] = {
            key: value for key in ("status", "valid", "compatible", "execution_status")
            if (value := source.get(key)) is None or isinstance(value, bool)
            or isinstance(value, str) and len(value) <= 80
        }
        result["result_summary"].update({
            f"{key}_count": len(source[key]) for key in ("warnings", "hard_conflicts", "limitations")
            if isinstance(source.get(key), (list, tuple))
        })
        if isinstance(payload.get("error"), Mapping):
            result["error"] = {key: str(payload["error"].get(key, ""))[:500] for key in ("type", "message")}
    return result

def _brief_fields(
    source: Mapping[str, Any], keys: tuple[str, ...], *, _depth: int = 0,
) -> dict[str, Any]:
    """Bound nested values, keeping counts from the saved source and literal IDs."""
    selected: dict[str, Any] = {}
    for key in keys:
        if key not in source:
            continue
        value = source[key]
        if key == "error" and value is None:
            continue
        if key in {"warnings", "cautions", "limitations"} and value in ([], (), None):
            continue
        if isinstance(value, str):
            # Structure strings and evidence locators must stay usable. If even
            # this larger bound is exceeded, explicitly label the partial text.
            literal = key.endswith(("smiles", "_id", "_ref", "_url", "_uri")) or key in {"expression", "query"}
            limit = 2000 if literal else 180
            selected[key] = value[:limit] + ("…" if len(value) > limit else "")
            if len(value) > limit:
                selected[f"{key}_truncated"] = True
        elif isinstance(value, (list, tuple)):
            if key.endswith(("smiles", "_ids", "_refs")):
                selected[key] = [item[:2000] + ("…" if len(item) > 2000 else "")
                                 if isinstance(item, str) else _brief_value(item, _depth + 1)
                                 for item in value[:3]]
                if any(isinstance(item, str) and len(item) > 2000 for item in value[:3]):
                    selected[f"{key}_truncated"] = True
            else:
                selected[key] = [_brief_value(item, _depth + 1) for item in value[:3]]
            if len(value) > 3:
                selected[f"{key}_count"] = len(value)
        elif isinstance(value, Mapping):
            selected[key] = _brief_value(value, _depth + 1)
            if len(value) > 5:
                selected[f"{key}_field_count"] = len(value)
        else:
            selected[key] = value
    return selected

def _brief_value(value: Any, depth: int) -> Any:
    """Bound arbitrary nested payloads inside explicitly selected brief fields."""
    if isinstance(value, str):
        return value[:180] + ("…" if len(value) > 180 else "")
    if isinstance(value, (Mapping, list, tuple)) and depth >= 3:
        return {"summary_omitted": True, "item_count": len(value)}
    if isinstance(value, Mapping):
        keys = tuple(key for key in islice(value, 5) if isinstance(key, str) and len(key) <= 80)
        result = _brief_fields(value, keys, _depth=depth)
        if len(value) > len(keys):
            result["omitted_field_count"] = len(value) - len(keys)
        return result
    if isinstance(value, (list, tuple)):
        return {"preview": [_brief_value(item, depth + 1) for item in value[:3]], "item_count": len(value)}
    return value

def _brief_recipe(recipe: Mapping[str, Any]) -> dict[str, Any]:
    """Show ingredients and reported operating fields, leaving provenance on disk."""
    groups = ("catalysts", "ligands", "bases", "acids", "oxidants", "reductants",
              "condensation_agents", "additives", "solvents", "other_components")
    components = [(group, index, item) for group in groups
                  for index, item in enumerate(recipe.get(group, [])) if isinstance(item, Mapping)]
    # An unresolved ingredient late in a long recipe must not disappear behind
    # six ordinary components. This orders the preview, not the recipe itself.
    components.sort(key=lambda entry: entry[2].get("identity_status") == "resolved" and not entry[2].get("warnings"))
    result = _brief_fields(recipe, (
        "temperature_c", "time_h", "concentration_m", "atmosphere", "declared_absences", "warnings",
    ))
    result = {key: value for key, value in result.items() if value is not None and value not in ([], ())}
    result["component_count"] = len(components)
    result["components"] = []
    for group, index, component in components[:6]:
        compact = _brief_fields(component, (
            "primary_role", "identity_status", "role_status", "warnings",
        ))
        compact.update(_brief_fields({"name": component.get("canonical_name") or component.get("raw_identifier")}, ("name",)))
        compact["saved_path"] = [group, index]
        if component.get("amount") is not None:
            compact.update(_brief_fields(component, ("amount", "amount_unit", "quantity_status")))
        result["components"].append(compact)
    if len(components) > 6:
        result["omitted_component_count"] = len(components) - 6
    if recipe.get("stages"):
        result["stage_count"] = len(recipe["stages"])
        result["stage_details_omitted"] = True
    return result
