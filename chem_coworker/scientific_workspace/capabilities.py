"""Local environment diagnostics, separate from observed agent-tool capabilities."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
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
        "agent_web_search": {"status": "not_checked",
                             "observation_location": "turns/<turn>/runtime-observations.json"},
        "model_and_reasoning": {"status": "not_confirmed",
                                "request_location": "turns/<turn>/runtime-request.json"},
    }
