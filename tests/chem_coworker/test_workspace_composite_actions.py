"""Agent composite expansion uses pinned artefacts, saved evidence and replay."""

from pathlib import Path

import pytest

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.core.baseline import (
    artifact_identity,
    code_manifest,
    environment_versions,
)
from core_retrosynthesis.composite_actions import (
    build_composite_strategy_catalog,
    save_composite_strategy_catalog,
)
from core_retrosynthesis.generic_library import save_generic_library
from tests.core_retrosynthesis_tests import test_composite_actions as composite_fixtures
from tests.core_retrosynthesis_tests.test_composite_actions import (
    ACTIVATION,
    SUBSTITUTION,
    candidate,
    strategy,
)


# Explicitly expose the shared chemistry fixture to pytest.
composite_library = composite_fixtures.composite_library


@pytest.fixture
def workspace(tmp_path, composite_library):
    library = tmp_path / "operators.json.gz"
    catalog = tmp_path / "catalog.json"
    save_generic_library(composite_library, library)
    save_composite_strategy_catalog(
        build_composite_strategy_catalog(
            (strategy(candidate(ACTIVATION), candidate(SUBSTITUTION)),)
        ),
        catalog,
    )
    root = Path(__file__).resolve().parents[2]
    baseline = {
        "repository": str(root),
        "code_files": code_manifest(root),
        "environment": environment_versions(),
        "artifacts": {
            "composite_library": artifact_identity(library),
            "composite_catalog": artifact_identity(catalog),
        },
    }
    InvestigationStore.create(
        tmp_path / "workspace",
        objective="Expand a coupled substitution",
        baseline=baseline,
    )
    return ScientificWorkspace(tmp_path / "workspace")


def test_composite_tool_records_both_physical_steps_and_replays(workspace):
    operation = next(
        item
        for item in workspace.operations.catalog()
        if item["name"] == "disconnect_composite"
    )
    assert operation["required_artifacts"] == ["composite_library", "composite_catalog"]
    event = workspace.run("disconnect_composite", {"target_smiles": "CCN", "top_k": 1})
    payload = workspace.store.read_artifact(event.artifact_ref)
    assert payload["execution_status"] == "completed", payload
    action = payload["result"]["actions"][0]
    assert action["physical_step_count"] == action["physical_step_cost"] == 2
    assert action["dependency"]["admitted"]
    assert all(
        step["evidence_kind"] == "predicted" for step in action["physical_steps"]
    )
    for detailed in (True, False):
        result = workspace.call_summary(event, detailed=detailed)["result_summary"]
        assert result["actions"][0]["physical_step_cost"] == 2
        assert len(result["actions"][0]["physical_steps"]) == 2
        assert result["actions"][0]["dependency"]["status"] == "verified"
    replay = workspace.replay(event.artifact_ref)
    assert workspace.store.read_artifact(replay.artifact_ref)["matches"]


def test_catalogue_change_is_detected_before_agent_execution(workspace):
    path = Path(
        workspace.store.manifest["baseline"]["artifacts"]["composite_catalog"]["path"]
    )
    path.write_text(path.read_text() + " ", encoding="utf-8")
    event = workspace.run("disconnect_composite", {"target_smiles": "CCN"})
    payload = workspace.store.read_artifact(event.artifact_ref)
    assert payload["execution_status"] != "completed"
