"""Literature route provenance must never manufacture structural support."""

from copy import deepcopy
from pathlib import Path

import pytest

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.adapters import literature
from chem_coworker.scientific_workspace.core.baseline import (
    artifact_identity,
    code_manifest,
    environment_versions,
)
from chem_coworker.scientific_workspace.core.store import canonical_bytes
from core_retrosynthesis.external_proposal_assessment import (
    ExternalRetrosynthesisProposal, assess_external_retrosynthesis_proposal,
)
from core_retrosynthesis.generic_library import build_generic_library, save_generic_library
from tests.core_retrosynthesis_tests.test_external_proposal_admission import (
    FIRST_REACTION, SECOND_REACTION, _route_value, _row,
)


@pytest.fixture(scope="module")
def library():
    return build_generic_library(
        (_row(FIRST_REACTION, 1), _row(SECOND_REACTION, 2)), levels=("L0", "L1", "L2"),
    )


@pytest.fixture
def workspace(tmp_path: Path, library) -> ScientificWorkspace:
    root = Path(__file__).resolve().parents[2]
    path = tmp_path / "operators.json.gz"
    save_generic_library(library, path)
    baseline = {
        "repository": str(root), "code_files": code_manifest(root), "environment": environment_versions(),
        "artifacts": {"retro_library": artifact_identity(path)},
    }
    InvestigationStore.create(tmp_path / "investigation", objective="Literature provenance regression", baseline=baseline)
    return ScientificWorkspace(tmp_path / "investigation")


def _call(workspace, operation, **arguments):
    event = workspace.run(operation, arguments)
    record = workspace.store.read_artifact(event.artifact_ref)
    assert record["execution_status"] == "completed", record
    return event, record["result"]


def _captured(workspace):
    return workspace.capture_source(
        "The authors claim this synthesis works. Experimental details remain to be checked.",
        url="https://example.org/paper", locator="Reported experimental section",
    )


@pytest.mark.parametrize("kind", ["fetched", "captured", "excerpt"])
@pytest.mark.parametrize("precursors", ["CC=O.N", "CC.CN"])
def test_literature_is_provenance_and_cannot_change_step_gates(workspace, library, monkeypatch, kind, precursors):
    if kind == "fetched":
        monkeypatch.setattr(literature, "_request", lambda *args, **kwargs: {
            "status": 200, "headers": {"content-type": "text/plain"},
            "body": b"The authors claim this synthesis works.",
        })
        source = workspace.fetch_source("https://example.org/paper")
    else:
        source = _captured(workspace)
    reference = source
    if kind == "excerpt":
        reference = workspace.record_source_excerpt(
            source.artifact_ref, excerpt="The authors claim this synthesis works.",
        )
    proposal = {"target_smiles": "CCN", "precursor_smiles": precursors}
    event, result = _call(workspace, "assess_route_step", proposal=proposal,
                          evidence_refs=[reference.artifact_ref, reference.artifact_ref])
    direct = assess_external_retrosynthesis_proposal(ExternalRetrosynthesisProposal.from_dict(proposal), library)
    assert canonical_bytes(result["assessment"]) == canonical_bytes(direct.to_dict())
    assert result["evidence_refs"] == [reference.artifact_ref]
    assert reference.artifact_ref in event.evidence_refs
    context = result["evidence_provenance"][0]
    assert context["source_ref"] == source.artifact_ref
    assert context["role"] == "literature_provenance"
    assert context["claim_support"] == "not_assessed"
    assert any("does not verify" in item for item in result["evidence_warnings"])
    if kind != "fetched":
        assert any("not been independently" in item for item in result["evidence_warnings"])
    assert result["experimental_feasibility"] == "not_established"
    replay = ScientificWorkspace(workspace.store.root).replay(event.artifact_ref)
    assert workspace.store.read_artifact(replay.artifact_ref)["matches"] is True


@pytest.mark.parametrize("mode", ["failed", "empty"])
def test_failed_fetch_and_empty_extraction_are_not_route_support(workspace, monkeypatch, mode):
    def request(*args, **kwargs):
        if mode == "failed":
            raise OSError("Development fixture: source unavailable")
        return {"status": 200, "headers": {"content-type": "text/plain"}, "body": b""}

    monkeypatch.setattr(literature, "_request", request)
    source = workspace.fetch_source("https://example.org/unavailable")
    event = workspace.run("assess_route_proposal", {
        "proposal": _route_value(), "evidence_refs": [source.artifact_ref],
    })
    record = workspace.store.read_artifact(event.artifact_ref)
    assert record["execution_status"] == "error"
    assert "debugging records" in record["error"]["message"]
    assert "capture_source" in record["error"]["message"]
    assert source.artifact_ref in event.evidence_refs


def test_partial_extraction_limitations_remain_visible(workspace, monkeypatch):
    monkeypatch.setattr(literature, "MAX_TEXT_CHARACTERS", 12)
    monkeypatch.setattr(literature, "_request", lambda *args, **kwargs: {
        "status": 200, "headers": {"content-type": "text/plain"},
        "body": b"Reported procedure followed by text beyond the extraction limit.",
    })
    source = workspace.fetch_source("https://example.org/long-source")
    _, result = _call(workspace, "assess_route_proposal", proposal=_route_value(),
                      evidence_refs=[source.artifact_ref])
    assert result["evidence_provenance"][0]["extraction_status"] == "partial_text_limit"
    assert any("partial_text_limit" in item for item in result["evidence_warnings"])
    assert any("omits molecular drawings" in item for item in result["evidence_warnings"])


@pytest.mark.parametrize("tamper", ["text", "source_kind", "schema"])
def test_excerpt_requires_recorded_source_and_exact_captured_text(workspace, tamper):
    source = _captured(workspace)
    excerpt = workspace.record_source_excerpt(source.artifact_ref, excerpt="The authors claim this synthesis works.")
    value = deepcopy(workspace.store.read_artifact(excerpt.artifact_ref))
    if tamper == "text":
        value["text"] = "Fabricated experiment"
    elif tamper == "source_kind":
        value["source_ref"] = workspace.store.note("hypothesis", "Unverified idea").artifact_ref
    else:
        value["schema_version"] = "unknown.v1"
    forged = workspace.store.append("literature_excerpt", value)
    event = workspace.run("assess_route_proposal", {
        "proposal": _route_value(), "evidence_refs": [forged.artifact_ref],
    })
    assert workspace.store.read_artifact(event.artifact_ref)["execution_status"] == "error"


def test_revision_inspection_and_comparison_preserve_literature_context(workspace):
    source = _captured(workspace)
    original, before = _call(workspace, "assess_route_proposal", proposal=_route_value(),
                             evidence_refs=[source.artifact_ref])
    revised, after = _call(workspace, "revise_route_branch", source_ref=original.artifact_ref,
                           remove_step_ids=["step-1"], replacement_steps=[{
                               "external_step_id": "step-1", "target_smiles": "CCN", "precursor_smiles": "CC.CN",
                           }], reason="Investigate a proposed alternative",
                           risks=["The source claim and the alternative remain unverified"])
    assert after["evidence_refs"] == [original.artifact_ref, source.artifact_ref]
    assert after["evidence_warnings"] == before["evidence_warnings"]
    assert after["assessment"]["status"] == "partially_supported"
    _, inspected = _call(workspace, "inspect_route_step", source_ref=revised.artifact_ref, step_id="step-1")
    assert inspected["evidence_provenance"] == after["evidence_provenance"]
    _, comparison = _call(workspace, "compare_route_proposals", source_refs=[original.artifact_ref, revised.artifact_ref])
    assert all(item["evidence_warnings"] == before["evidence_warnings"] for item in comparison["alternatives"])
    replay = ScientificWorkspace(workspace.store.root).replay(revised.artifact_ref)
    assert workspace.store.read_artifact(replay.artifact_ref)["matches"] is True
