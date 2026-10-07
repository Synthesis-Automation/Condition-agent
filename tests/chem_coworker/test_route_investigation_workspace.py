"""Recorded proposed-route assessment, branch revision, recovery and comparison."""

from copy import deepcopy
from pathlib import Path

import pytest

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.core.baseline import (
    artifact_identity,
    code_manifest,
    environment_versions,
)
from chem_coworker.scientific_workspace.core.store import canonical_bytes
from core_retrosynthesis.external_proposal_assessment import (
    ExternalRetrosynthesisProposal,
    assess_external_retrosynthesis_proposal,
)
from core_retrosynthesis.external_route_admission import (
    ExternalRouteProposal,
    assess_external_route_proposal,
)
from core_retrosynthesis.generic_library import build_generic_library, save_generic_library
from tests.core_retrosynthesis_tests.test_external_proposal_admission import (
    _row,
    _route_value,
    FIRST_REACTION,
    SECOND_REACTION,
)


@pytest.fixture(scope="module")
def library():
    return build_generic_library(tuple(_row(reaction, ordinal) for ordinal, reaction in enumerate(
        (FIRST_REACTION, SECOND_REACTION, "CCBr.N>>CCN"), 1,
    )), levels=("L0", "L1", "L2"))


@pytest.fixture
def workspace(tmp_path: Path, library) -> ScientificWorkspace:
    root = Path(__file__).resolve().parents[2]
    path = tmp_path / "operators.json.gz"
    save_generic_library(library, path)
    baseline = {"repository": str(root), "code_files": code_manifest(root), "environment": environment_versions(),
                "artifacts": {"retro_library": artifact_identity(path)}}
    InvestigationStore.create(tmp_path / "investigation", objective="Development route revision", baseline=baseline)
    return ScientificWorkspace(tmp_path / "investigation")


def call(workspace, operation, **arguments):
    event = workspace.run(operation, arguments)
    saved = workspace.store.read_artifact(event.artifact_ref)
    assert saved["execution_status"] == "completed", saved
    return event, saved["result"]


@pytest.mark.parametrize("proposal", [
    {"target_smiles": "CCN", "precursor_smiles": "CC=O.N"},
    {"target_smiles": "CCN", "precursor_smiles": "N.CC=O"},
    {"target_smiles": "CCN", "precursor_smiles": "CC.CN"},
    {"target_smiles": "invalid", "precursor_smiles": "CCO"},
    {"target_smiles": "CCO", "precursor_smiles": "CC=O.N",
     "mapped_reaction_smiles": "[CH3:1][CH:2]=[O:3].[NH3:4]>>[CH3:1][CH2:2][NH2:4]"},
])
def test_step_adapter_preserves_positive_negative_ambiguous_and_conflicting_evidence(workspace, library, proposal) -> None:
    _, result = call(workspace, "assess_route_step", proposal=proposal)
    direct = assess_external_retrosynthesis_proposal(ExternalRetrosynthesisProposal.from_dict(proposal), library)
    assert canonical_bytes(result["assessment"]) == canonical_bytes(direct.to_dict())
    assert result["experimental_feasibility"] == "not_established"
    assert result["assessment"]["actionable"] is False


def test_route_assessment_and_inspection_preserve_domain_evidence(workspace, library) -> None:
    event, record = call(workspace, "assess_route_proposal", proposal=_route_value())
    direct = assess_external_route_proposal(ExternalRouteProposal.from_dict(_route_value()), library)
    assert canonical_bytes(record["assessment"]) == canonical_bytes(direct.to_dict())
    inspected, view = call(workspace, "inspect_route_step", source_ref=event.artifact_ref, step_id="step-1")
    assert view["downstream_step_ids"] == ["step-2"]
    assert "CCN" in view["molecule_audits"]
    assert view["assessment"]["precedent_matches"]
    assert inspected.evidence_refs == (event.artifact_ref,)


def test_branch_revision_rechecks_preserved_downstream_step_and_replays_after_reopen(workspace, monkeypatch) -> None:
    source, before = call(workspace, "assess_route_proposal", proposal=_route_value(), unavailable_starting_materials=["CC=O"])
    source_bytes = canonical_bytes(workspace.store.read_artifact(source.artifact_ref))
    from core_retrosynthesis import external_route_admission

    original = external_route_admission.assess_external_retrosynthesis_proposal
    checked = []

    def spy(proposal, *args, **kwargs):
        checked.append(proposal.target_smiles)
        return original(proposal, *args, **kwargs)

    monkeypatch.setattr(external_route_admission, "assess_external_retrosynthesis_proposal", spy)
    reopened = ScientificWorkspace(workspace.store.root)
    revision, after = call(reopened, "revise_route_branch", source_ref=source.artifact_ref,
                          remove_step_ids=["step-1"], replacement_steps=[{
                              "external_step_id": "step-1", "target_smiles": "CCN", "precursor_smiles": "CCBr.N",
                          }], reason="Declared aldehyde is unavailable; investigate an alternative amine preparation",
                          risks=["Overalkylation and experimental feasibility remain unverified"])
    assert len(checked) == 2 and _route_value()["target_smiles"] in checked
    assert before["material_constraints"]["status"] == "violated"
    assert after["material_constraints"]["status"] == "satisfied_for_declared_constraints"
    assert after["revision"]["preserved_step_ids"] == ["step-2"]
    assert after["revision"]["improvement_status"] == "not_automatically_established"
    assert after["assessment_options"] == before["assessment_options"]
    assert canonical_bytes(workspace.store.read_artifact(source.artifact_ref)) == source_bytes
    assert source.artifact_ref in revision.evidence_refs
    compare, table = call(reopened, "compare_route_proposals", source_refs=[source.artifact_ref, revision.artifact_ref])
    assert len(table["alternatives"]) == 2 and table["ranking"] == "not_performed"
    assert compare.evidence_refs == (source.artifact_ref, revision.artifact_ref)
    replay = ScientificWorkspace(workspace.store.root).replay(revision.artifact_ref)
    assert workspace.store.read_artifact(replay.artifact_ref)["matches"] is True


def test_unsupported_revision_remains_inspectable_and_noop_is_an_error(workspace) -> None:
    source, _ = call(workspace, "assess_route_proposal", proposal=_route_value())
    unknown, record = call(workspace, "revise_route_branch", source_ref=source.artifact_ref,
                          remove_step_ids=["step-1"], replacement_steps=[{
                              "external_step_id": "step-1", "target_smiles": "CCN", "precursor_smiles": "CC.CN",
                          }], reason="Investigate an unsupported alternative", risks=["No validated correspondence"])
    assert record["assessment"]["status"] == "partially_supported"
    assert "step-1" in record["assessment"]["unresolved_step_ids"]
    _, broken = call(workspace, "revise_route_branch", source_ref=unknown.artifact_ref,
                     remove_step_ids=["step-2"], replacement_steps=[], reason="Demonstrate disconnected topology",
                     risks=["Target no longer produced"])
    assert broken["assessment"]["status"] == "invalid"
    assert broken["material_constraints"]["status"] == "unknown"
    error = workspace.run("revise_route_branch", {
        "source_ref": source.artifact_ref, "remove_step_ids": [], "replacement_steps": [],
        "reason": "No change", "risks": ["No change"],
    })
    assert workspace.store.read_artifact(error.artifact_ref)["execution_status"] == "error"


def test_unknown_fields_bad_sources_and_mismatched_comparisons_are_rejected(workspace) -> None:
    value = _route_value()
    altered = deepcopy(value)
    altered["steps"][0]["confidence"] = 1
    for proposal in (altered, {**value, "admission_status": "verified"}):
        error = workspace.run("assess_route_proposal", {"proposal": proposal})
        assert workspace.store.read_artifact(error.artifact_ref)["execution_status"] == "error"
    good, _ = call(workspace, "assess_route_proposal", proposal=value)
    constrained, _ = call(workspace, "assess_route_proposal", proposal=value, unavailable_starting_materials=["CCN"])
    failed = workspace.run("compare_route_proposals", {"source_refs": [good.artifact_ref, constrained.artifact_ref]})
    assert "same options" in workspace.store.read_artifact(failed.artifact_ref)["error"]["message"]
    note = workspace.store.note("hypothesis", "Not scientific execution evidence")
    failed = workspace.run("assess_route_proposal", {"proposal": value, "evidence_refs": [note.artifact_ref]})
    assert workspace.store.read_artifact(failed.artifact_ref)["execution_status"] == "error"


def test_supplied_recipe_is_assessed_separately_and_forward_challenge_is_optional(workspace) -> None:
    _, recipe = call(workspace, "resolve_recipe", components=[{
        "raw_identifier": "ethanol", "source_field": "user", "identifier_type": "name",
    }])
    _, result = call(workspace, "assess_route_step", proposal={
        "target_smiles": "CCN", "precursor_smiles": "CC=O.N", "proposed_conditions": recipe,
    })
    assert result["proposed_recipe_assessment"]["schema_version"] == "reaction_recipe_assessment.v3"
    assert result["assessment"]["forward_assessment"] is None
    assert next(gate for gate in result["assessment"]["gates"] if gate["gate_id"] == "condition_support")["status"] == "unresolved"


def test_material_constraint_cannot_be_silently_removed_during_revision(workspace) -> None:
    source, _ = call(workspace, "assess_route_proposal", proposal=_route_value(), unavailable_starting_materials=["CC=O"])
    error = workspace.run("revise_route_branch", {
        "source_ref": source.artifact_ref, "remove_step_ids": [], "replacement_steps": [],
        "reason": "Attempt to change constraints", "risks": ["Unresolved"], "unavailable_starting_materials": [],
    })
    assert workspace.store.read_artifact(error.artifact_ref)["execution_status"] == "error"


def test_requested_condition_retrieval_requires_recorded_data(workspace) -> None:
    error = workspace.run("assess_route_step", {
        "proposal": {"target_smiles": "CCN", "precursor_smiles": "CC=O.N"}, "include_conditions": True,
    })
    saved = workspace.store.read_artifact(error.artifact_ref)
    assert saved["execution_status"] == "error"
    assert "condition_index" in saved["error"]["message"]


@pytest.mark.parametrize("target", ["CCN", "[He]", "invalid"])
def test_single_step_preserves_engine_evidence_without_stock_or_internal_review(workspace, library, target, monkeypatch) -> None:
    from chem_coworker.contracts import RetrosynthesisRequest
    from chem_coworker.retrosynthesis import RetrosynthesisCoworker
    from chem_coworker.multistep import MultistepRetrosynthesisCoworker

    def forbidden(*args, **kwargs):
        pytest.fail("Single-step workspace must not invoke a multistep planner or hidden reviewer")

    monkeypatch.setattr(MultistepRetrosynthesisCoworker, "plan", forbidden)
    event, result = call(workspace, "disconnect_target", target_smiles=target)
    direct = RetrosynthesisCoworker(
        library=library, library_path=workspace.operations._path("retro_library"),
    ).disconnect(RetrosynthesisRequest(
        target_smiles=target, top_k=3, max_templates_to_apply=40,
        max_candidates_to_validate=10, include_conditions=False,
    ))
    assert canonical_bytes(result) == canonical_bytes(direct.to_dict())
    assert result["request"]["review"]["mode"] == "off"
    assert result["review"] is None
    assert "stock_index" not in workspace.store.manifest["baseline"]["artifacts"]
    if target == "CCN":
        assert result["strategies"]
        summary = workspace.call_summary(event)["result_summary"]
        candidate = summary["strategies"][0]["representative"]
        assert candidate["precursor_smiles"] == result["strategies"][0]["representative"]["precursor_smiles"]
        assert candidate["forward_validation_status"] == "verified_signature"
        assert summary["warnings"] == result["warnings"]
    else:
        assert not result["strategies"]
    replay = ScientificWorkspace(workspace.store.root).replay(event.artifact_ref)
    assert workspace.store.read_artifact(replay.artifact_ref)["matches"] is True


def test_agent_selects_next_intermediate_and_assesses_assembled_route(workspace) -> None:
    from rdkit import Chem

    def canonical(smiles: str) -> str:
        return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    target = _route_value()["target_smiles"]
    first, initial = call(workspace, "disconnect_target", target_smiles=target)
    chosen = next(candidate for strategy in initial["strategies"]
                  for candidate in [strategy["representative"], *strategy["alternate_realizations"]]
                  if "CCN" in candidate["precursor_smiles"].split("."))
    # Caller explicitly chooses a leaf; the first call has not expanded it.
    second, expanded = call(workspace, "disconnect_target", target_smiles="CCN")
    amine = next(strategy["representative"] for strategy in expanded["strategies"]
                 if canonical(strategy["representative"]["precursor_smiles"]) == canonical("CC=O.N"))
    event, route = call(workspace, "assess_route_proposal", proposal={
        "target_smiles": target,
        "steps": [
            {"external_step_id": "amide", "target_smiles": chosen["target_smiles"],
             "precursor_smiles": chosen["precursor_smiles"]},
            {"external_step_id": "amine", "target_smiles": amine["target_smiles"],
             "precursor_smiles": amine["precursor_smiles"]},
        ],
    }, evidence_refs=[first.artifact_ref, second.artifact_ref])
    assert len(route["assessment"]["step_assessments"]) == 2
    assert route["assessment"]["status"] == "admitted_review_only"
    assert event.evidence_refs == (first.artifact_ref, second.artifact_ref)


@pytest.mark.parametrize("operation", ["plan_routes", "revise_routes", "prepare_route_proposal"])
def test_workspace_rejects_automatic_route_operations(workspace, operation) -> None:
    assert operation not in {item["name"] for item in workspace.operations.catalog()}
    error = workspace.run(operation, {})
    assert "Unknown scientific operation" in workspace.store.read_artifact(error.artifact_ref)["error"]["message"]


@pytest.mark.parametrize("extra", [{"max_depth": 3}, {"beam_width": 6}, {"review": {"mode": "always"}}, {"top_k": 0}])
def test_single_step_rejects_planner_review_and_invalid_options(workspace, extra) -> None:
    error = workspace.run("disconnect_target", {"target_smiles": "CCN", **extra})
    assert workspace.store.read_artifact(error.artifact_ref)["execution_status"] == "error"


def test_single_step_conditions_are_opt_in_and_require_recorded_data(workspace) -> None:
    error = workspace.run("disconnect_target", {"target_smiles": "CCN", "include_conditions": True})
    assert "condition_index" in workspace.store.read_artifact(error.artifact_ref)["error"]["message"]


def test_prompt_assigns_multistep_decisions_to_agent(workspace) -> None:
    from chem_coworker.scientific_workspace.agent_context.prompts import investigation_prompt

    prompt = investigation_prompt(workspace, "Propose a synthesis")
    assert "Do not invoke the built-in multistep planner, including through custom Python scripts." in prompt
    assert "assess_route_proposal" in prompt
    assert "single-step" in prompt
    assert "Do not read answer-schema.json at startup" in prompt
    assert "decision-critical field" in prompt
    assert "Repeat a search only with a" in " ".join(prompt.split())
    for retired in ("plan_routes", "revise_routes", "prepare_route_proposal", "beam_width"):
        assert retired not in prompt


def test_prompt_compact_route_example_validates_against_real_contract(workspace) -> None:
    import ast

    from chem_coworker.scientific_workspace.answers.answer_contracts import ScientificAnswer
    from chem_coworker.scientific_workspace.agent_context.prompts import investigation_prompt

    prompt = investigation_prompt(workspace, "Propose a synthesis")
    snippet = next(block.split("\n```", 1)[0] for block in prompt.split("```python\n")
                   if block.startswith("draft = w.route_answer("))
    source = workspace.capture_source("Example 1: test illustration only.",
                                      url="https://example.org/paper", locator="Example 1")
    # Evaluate only the trusted documentation's input expression, not its file handoff.
    expression = ast.Expression(ast.parse(snippet).body[0].value.args[0])
    request = eval(compile(expression, "prompt-example", "eval"), {"__builtins__": {}}, {
        "target_smiles": "CC=O", "reactants": "CCO", "product": "CC=O",
        "captured_source_ref": source.artifact_ref,
    })
    answer = ScientificAnswer.model_validate(workspace.route_answer(request))
    assert answer.steps[0].conditions[0].text == "Conditions to develop"
    assert answer.steps[0].conditions[0].basis == "unknown"
    assert answer.routes[0].step_ids == [answer.steps[0].id]
    assert answer.sources[0].artifact_ref == source.artifact_ref
    assert not answer.steps[0].literature_reactions
