"""Local environment diagnostics, separate from observed agent-tool capabilities."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any, Mapping

from .baseline import environment_versions


def local_capabilities(baseline: Mapping[str, Any]) -> dict[str, Any]:
    """Probe RDKit and input presence without loading datasets or contacting a provider.

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
    return {
        "schema_version": "scientific_capabilities.v1",
        "python": {"executable": sys.executable, "status": "running"},
        "versions": environment_versions(),
        "rdkit": rdkit,
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
