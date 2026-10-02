"""Selected-step forward checks use saved evidence and killable, pinned inputs."""

from dataclasses import replace
import json
import sys
from time import monotonic

import pytest

from chem_coworker.scientific_workspace import InvestigationStore
from chem_coworker.scientific_workspace.core.baseline import artifact_identity
from chem_coworker.scientific_workspace.adapters import forward_check, route_investigation
from chem_coworker.scientific_workspace.adapters.operations import ScientificOperations
from chem_coworker.scientific_workspace.core.store import canonical_bytes
from core_retrosynthesis.external_proposal_assessment import load_external_proposal_admission_policy
from core_retrosynthesis.generic_library import build_generic_library, save_generic_library
from forward_synthesis import assess_proposed_step, build_forward_library, save_forward_library
from tests.core_retrosynthesis_tests.test_external_proposal_admission import (
    _row, _route_value, FIRST_REACTION, SECOND_REACTION,
)


@pytest.fixture(scope="module")
def libraries():
    generic = build_generic_library(tuple(_row(reaction, index) for index, reaction in enumerate(
        (FIRST_REACTION, SECOND_REACTION), 1,
    )), levels=("L0", "L1", "L2"))
    return generic, build_forward_library(generic)


@pytest.fixture
def operations(tmp_path, libraries):
    generic_path, forward_path = tmp_path / "retro.json.gz", tmp_path / "forward.json.gz"
    save_generic_library(libraries[0], generic_path)
    save_forward_library(libraries[1], forward_path)
    baseline = {"artifacts": {"retro_library": artifact_identity(generic_path),
                              "forward_library": artifact_identity(forward_path)}}
    return ScientificOperations(InvestigationStore.create(
        tmp_path / "investigation", objective="Check one proposed step", baseline=baseline,
    ))


def _source(operations, proposal=None, *, operation="assess_route_proposal", options=None):
    record = route_investigation.assess_route(
        operations, proposal or _route_value(), None, False, False, None,
    )
    if options is not None:
        record["assessment_options"] = options
    event = operations.store.append("call", {"operation": operation, "execution_status": "completed", "result": record})
    return event.artifact_ref, record


def _check(operations, source_ref, **kwargs):
    return forward_check.assess_step_forward(
        operations, source_ref, kwargs.pop("step_id", "step-1"),
        kwargs.pop("question", "Does the intended product compete with another predicted product?"), **kwargs,
    )


def test_forward_check_matches_canonical_assessment_and_preserves_source(operations, libraries):
    source, record = _source(operations)
    before = canonical_bytes(operations.store.read_artifact(source))
    result = _check(operations, source)
    assert result["execution_status"] == "completed", result
    selected = next(item["assessment"] for item in record["assessment"]["step_assessments"]
                    if item["external_step_id"] == "step-1")
    match = next(item for item in selected["operator_matches"]
                 if item["match_id"] == selected["selected_operator_match_id"])
    limits = load_external_proposal_admission_policy().limits
    direct = assess_proposed_step(
        selected["canonical_precursor_smiles"], selected["canonical_target_smiles"], libraries[1],
        operator_hint=match["operator_id"], top_k=limits.maximum_forward_products,
        max_operators_to_apply=limits.maximum_forward_operators,
    )
    assert canonical_bytes(result["assessment"]) == canonical_bytes(direct.to_dict())
    assert result["experimental_feasibility"] == "not_established"
    assert result["provenance"]["library_source_verified"] is True
    assert result["provenance"]["input_hashes_verified"] is True
    assert result["provenance"]["source_assessment_id"] == selected["assessment_id"]
    assert result["provenance"]["inputs"] == operations.store.manifest["baseline"]["artifacts"]
    assert canonical_bytes(operations.store.read_artifact(source)) == before
    assert {stage["stage"] for stage in result["execution"]["stages"]} == {
        "verify_inputs", "load_retro_library", "load_forward_library", "validate_library_source",
        "assess_selected_step", "verify_inputs_after_assessment",
    }


def test_missing_forward_library_does_not_build_or_start_worker(operations, monkeypatch):
    source, _ = _source(operations)
    del operations.store.manifest["baseline"]["artifacts"]["forward_library"]
    monkeypatch.setattr(forward_check.subprocess, "Popen", lambda *a, **k: pytest.fail("Worker must not start"))
    with pytest.raises(FileNotFoundError, match="forward_library"):
        _check(operations, source)


@pytest.mark.parametrize("step_id,precursors", [("missing-step", None), ("step-1", "CC.CN")])
def test_unknown_or_ineligible_step_does_not_start_worker(operations, monkeypatch, step_id, precursors):
    proposal = _route_value()
    if precursors:
        next(item for item in proposal["steps"] if item["external_step_id"] == "step-1")["precursor_smiles"] = precursors
    source, _ = _source(operations, proposal)
    monkeypatch.setattr(forward_check.subprocess, "Popen", lambda *a, **k: pytest.fail("Worker must not start"))
    with pytest.raises(ValueError, match="previously assessed|exact admitted"):
        _check(operations, source, step_id=step_id)


def test_source_must_be_recorded_completed_route_call(operations):
    source, _ = _source(operations)
    payload = operations.store.read_artifact(source)
    payload["untrusted_attachment"] = True
    attachment = operations.store.append("derived_file", payload)
    with pytest.raises(ValueError, match="recorded route"):
        _check(operations, attachment.artifact_ref)
    failed = operations.store.append("call", {**payload, "execution_status": "error"})
    with pytest.raises(ValueError, match="completed route"):
        _check(operations, failed.artifact_ref)


def test_deadline_kills_worker_during_load_and_retains_stage_without_assessment(operations, monkeypatch):
    source, _ = _source(operations)
    processes = []
    original = forward_check.subprocess.Popen

    def spawn(*args, **kwargs):
        process = original(*args, **kwargs)
        processes.append(process)
        return process

    monkeypatch.setattr(forward_check.subprocess, "Popen", spawn)
    monkeypatch.setattr(forward_check, "_worker_command", lambda directory: [sys.executable, "-c", (
        "import pathlib,sys,time; p=pathlib.Path(sys.argv[1]); "
        "(p/'stages.jsonl').write_text('{\"stage\":\"load_forward_library\",\"status\":\"running\"}\\n'); "
        "time.sleep(60); (p/'result.json').write_text('{}')"
    ), str(directory)])
    started = monotonic()
    result = _check(operations, source, timeout_seconds=1)
    assert monotonic() - started < 8
    assert processes and processes[0].poll() is not None
    assert result["execution_status"] == "timed_out"
    assert result["assessment"] is None
    assert result["error"]["type"] == "TimeoutExpired"
    assert result["execution"]["stages"] == [{"stage": "load_forward_library", "status": "running"}]
    assert "result.json" not in result["execution"]["diagnostics"]


def test_platform_cleanup_error_does_not_lose_timeout_or_leave_worker(operations, monkeypatch):
    from chem_coworker.scientific_workspace.core import process_utils

    source, _ = _source(operations)
    processes = []
    original = forward_check.subprocess.Popen

    def spawn(*args, **kwargs):
        process = original(*args, **kwargs)
        processes.append(process)
        return process

    def broken_cleanup(process):
        raise ProcessLookupError("Simulated platform cleanup race")

    monkeypatch.setattr(forward_check.subprocess, "Popen", spawn)
    monkeypatch.setattr(process_utils, "stop_process_tree", broken_cleanup)
    monkeypatch.setattr(forward_check, "_worker_command", lambda directory: [sys.executable, "-c", "import time; time.sleep(60)"])
    result = _check(operations, source, timeout_seconds=1)
    assert result["execution_status"] == "timed_out"
    assert result["error"]["type"] == "TimeoutExpired"
    assert result["execution"]["cleanup_error"]["type"] == "ProcessLookupError"
    assert processes[0].poll() is not None


@pytest.mark.parametrize("output", [None, "[]", "not-json"])
def test_worker_exit_without_valid_result_preserves_error_diagnostics(operations, monkeypatch, output):
    source, _ = _source(operations)
    monkeypatch.setattr(forward_check, "_worker_command", lambda directory: [sys.executable, "-c", (
        "import pathlib,sys; p=pathlib.Path(sys.argv[1]); print('worker diagnostic',file=sys.stderr); "
        + (f"(p/'result.json').write_text({output!r})" if output is not None else "pass")
    ), str(directory)])
    result = _check(operations, source)
    assert result["execution_status"] == "error"
    assert result["assessment"] is None
    assert result["error"]["type"] in {"FileNotFoundError", "ValueError", "JSONDecodeError"}
    stderr = operations.store.root / result["execution"]["diagnostics"]["stderr.log"]
    assert "worker diagnostic" in stderr.read_text()


def test_prebuilt_input_change_is_an_error_not_a_prediction(operations):
    source, _ = _source(operations)
    operations._path("forward_library").write_bytes(b"changed")
    result = _check(operations, source)
    assert result["execution_status"] == "error"
    assert "Pinned forward_library artifact changed" in result["error"]["message"]
    assert result["assessment"] is None
    assert result["execution"]["stages"][-1]["stage"] == "verify_inputs"


def test_library_with_wrong_source_is_rejected_before_prediction(operations, libraries):
    source, _ = _source(operations)
    forward_path = operations._path("forward_library")
    save_forward_library(replace(libraries[1], source_library_definition_id="wrong-source"), forward_path)
    operations.store.manifest["baseline"]["artifacts"]["forward_library"] = artifact_identity(forward_path)
    result = _check(operations, source)
    assert result["execution_status"] == "error"
    assert result["assessment"] is None
    assert result["execution"]["stages"][-1]["stage"] == "validate_library_source"


def test_worker_never_builds_a_forward_library(operations, monkeypatch):
    import forward_synthesis
    from core_retrosynthesis import forward_assessment

    source, _ = _source(operations)
    # Run the exact worker entry point in this interpreter to make the guard
    # effective, using the request prepared by the real parent path.
    monkeypatch.setattr(forward_synthesis, "build_forward_library", lambda *a, **k: pytest.fail("No runtime build"))
    monkeypatch.setattr(forward_assessment, "build_forward_library_from_generic", lambda *a, **k: pytest.fail("No runtime build"))
    monkeypatch.setattr(forward_check, "_worker_command", lambda directory: [sys.executable, "-c", "pass"])
    incomplete = _check(operations, source)
    directory = (operations.store.root / incomplete["execution"]["diagnostics"]["request.json"]).parent
    assert forward_check._run_worker(directory) == 0
    result = json.loads((directory / "result.json").read_text())
    assert result["execution_status"] == "completed"


@pytest.mark.parametrize("timeout", [0, 31, True, 1.2])
def test_deadline_requires_bounded_integer(operations, timeout):
    with pytest.raises(ValueError, match="integer between 1 and 30"):
        _check(operations, "unused", timeout_seconds=timeout)


@pytest.mark.parametrize("operation,proposal", [
    ("assess_route_step", {"target_smiles": "CCN", "precursor_smiles": "CC=O.N"}),
    ("assess_route_proposal", _route_value()),
])
def test_standard_assessment_rejects_forward_flag_before_loading(operations, monkeypatch, operation, proposal):
    monkeypatch.setattr(operations, "_external_route_library", lambda: pytest.fail("No library may load"))
    with pytest.raises(ValueError, match="assess_route_step_forward"):
        operations.invoke(operation, {"proposal": proposal, "include_forward": True})


def test_revision_does_not_inherit_old_expensive_forward_option(operations, monkeypatch):
    source, _ = _source(operations, options={"include_conditions": False, "include_forward": True})
    original = route_investigation.assess_route
    observed = []

    def spy(*args, **kwargs):
        observed.append(kwargs["include_forward"])
        return original(*args, **kwargs)

    monkeypatch.setattr(route_investigation, "assess_route", spy)
    result = route_investigation.revise_branch(
        operations, source, ["step-1"], [{
            "external_step_id": "step-1", "target_smiles": "CCN", "precursor_smiles": "CCBr.N",
        }], "Recheck alternative precursor", None, None, ["Feasibility remains unknown"],
    )
    assert observed == [False]
    assert result["assessment_options"]["include_forward"] is False
    assert all(step["assessment"]["forward_assessment"] is None for step in result["assessment"]["step_assessments"])
