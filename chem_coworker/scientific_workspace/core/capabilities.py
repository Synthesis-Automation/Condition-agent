"""Local environment diagnostics, separate from observed agent-tool capabilities."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any, Mapping

from .baseline import environment_versions


def local_capabilities(baseline: Mapping[str, Any]) -> dict[str, Any]:
    """Probe RDKit, input presence, and index metadata without scanning datasets.

    File presence is deliberately not index compatibility, admission, stock
    availability, or a claim that the model can execute a particular tool.
    """
    try:
        from rdkit import Chem

        molecule = Chem.MolFromSmiles("CCO")
        rdkit = {"status": "available" if molecule is not None else "failed",
                 "probe": "parse_ethanol", "experimental_validation": False}
    except Exception as exc:
        rdkit = {"status": "failed", "error_type": type(exc).__name__, "message": str(exc)}
    # Read prerequisites from the operation catalog rather than duplicating them.
    from ..adapters.operations import ScientificOperations

    artifacts = baseline.get("artifacts", {})
    index_compatibility = {}
    condition_path = artifacts.get("condition_index", {}).get("path")
    if condition_path and Path(condition_path).is_file():
        if Path(condition_path).suffix.lower() in {".sqlite", ".sqlite3", ".db"}:
            from condition_recommender.sqlite_indexing import load_sqlite_generic_index

            try:
                # This loader validates metadata and returns lazy row/lookups;
                # it does not scan the corpus or bypass the normal checks.
                load_sqlite_generic_index(condition_path)
            except (ValueError, OSError) as exc:
                index_compatibility["condition_index"] = {
                    "status": "incompatible_or_invalid",
                    "error": {"type": type(exc).__name__, "message": str(exc)},
                    "scope": "metadata_contract_only",
                }
            else:
                index_compatibility["condition_index"] = {
                    "status": "compatible_metadata", "scope": "metadata_contract_only",
                }
        else:
            index_compatibility["condition_index"] = {"status": "not_checked", "scope": "non_sqlite_artifact"}
    operation_readiness = {}
    for definition in ScientificOperations.DEFINITIONS:
        missing = [name for name in definition.required_artifacts
                   if name not in artifacts or not Path(artifacts[name]["path"]).is_file()]
        incompatible = [name for name in definition.required_artifacts
                        if index_compatibility.get(name, {}).get("status") == "incompatible_or_invalid"]
        operation_readiness[definition.name] = {
            "required_artifacts": list(definition.required_artifacts),
            "missing_artifacts": missing,
            "incompatible_artifacts": incompatible,
            "status": ("missing_artifacts" if missing else "incompatible_artifacts" if incompatible
                       else "prerequisites_present_not_validated"),
            "content_validation": "checked_on_call",
        }
    if index_compatibility.get("condition_index", {}).get("status") == "incompatible_or_invalid":
        operation_readiness["inspect_step_precedents"]["optional_enrichment"] = {
            "condition_observations": "unavailable", "template_and_source_inspection": "available_if_configured",
        }
    return {
        "schema_version": "scientific_capabilities.v1",
        "python": {"executable": sys.executable, "status": "running"},
        "versions": environment_versions(),
        "rdkit": rdkit,
        "operation_readiness": operation_readiness,
        "index_compatibility": index_compatibility,
        "pdf_text_extraction": {
            "status": "available" if importlib.util.find_spec("pypdf") else "not_installed",
            "ocr": "not_available", "document_access": "not_checked",
        },
        "artifacts": {
            name: {"status": "file_present" if Path(value["path"]).is_file() else "missing",
                   "content_validation": "not_checked", "baseline_status": value.get("status")}
            for name, value in baseline.get("artifacts", {}).items()
        },
        "source_fetch": {"status": "implemented", "network_access": "not_checked"},
        "fragment_search": {
            "status": "implemented", "required_artifact": "fragment_index",
            "content_validation": "checked_on_call", "automatic_build": False,
        },
        "fragment_suggestions": {
            "status": "implemented", "requires_index": False, "automatic_search": False,
        },
        "fragment_investigation": {
            "status": "implemented", "query_operation": "propose_fragment_queries",
            "source_operation": "investigate_fragment_precedent",
            "requires_saved_target_search": True, "requires_retro_library": False,
        },
        "starting_material_assessment": {
            "status": "implemented", "operation": "assess_starting_material",
            "optional_artifact": "fragment_index", "exact_product_lookup": True,
            "availability_verified": False, "stops_are_planning_assumptions": True,
        },
        "weak_label_screening": {
            "status": "implemented",
            "operation": "generate_weak_label_screening_array",
            "required_artifacts": ["weak_label_records", "weak_label_recipe_catalog"],
            "content_validation": "checked_on_call", "requires_condition_index": False,
            "source_structures_verified": False,
        },
        "molecular_inspection": {
            "status": "implemented", "requires_index": False,
            "operations": ["compare_molecules", "inspect_reactive_sites"],
            "experimental_selectivity_prediction": False,
        },
        "retro_validity": {
            "status": "implemented", "operation": "assess_retro_validity",
            "required_artifact": "retro_library",
            "optional_artifacts": ["condition_index", "shared_core_index"],
            "forward_check": "join_saved_bounded_audit",
            "ranking_influence": "none_advisory_only", "calibrated_success_probability": False,
        },
        "agent_web_search": {"status": "not_checked",
                             "observation_location": "turns/<turn>/runtime-observations.json"},
        "model_and_reasoning": {"status": "not_confirmed",
                                "request_location": "turns/<turn>/runtime-request.json"},
    }
