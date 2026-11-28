"""Corpus support is recorded for the exact saved route step the agent selects."""

from pathlib import Path

import pytest

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.core.baseline import artifact_identity
from chem_coworker.scientific_workspace.core.store import canonical_bytes
from condition_recommender import assess_step_plausibility
from condition_recommender.generic_indexing import build_generic_index
from condition_recommender.shared_core_index import build_shared_core_index
from condition_recommender.sqlite_indexing import load_sqlite_generic_index, save_sqlite_generic_index
from tests.chem_coworker import test_route_investigation_workspace as route_tests
from tests.chem_coworker.test_route_investigation_workspace import call
from tests.condition_recommender.test_shared_core_retrieval import record
from tests.core_retrosynthesis_tests.test_external_proposal_admission import _route_value


library = route_tests.library
workspace = route_tests.workspace


def test_replay_preserves_selected_route_step_and_complete_domain_evidence(
    workspace: ScientificWorkspace, tmp_path: Path,
) -> None:
    reaction = "CC=O.N>>CCN"
    index_path = tmp_path / "index.sqlite"
    shared_path = tmp_path / "shared.sqlite"
    save_sqlite_generic_index(build_generic_index([record(1, reaction)]), index_path)
    build_shared_core_index(load_sqlite_generic_index(index_path), shared_path)
    baseline = dict(workspace.store.manifest["baseline"])
    baseline["artifacts"] = {
        **baseline["artifacts"], "condition_index": artifact_identity(index_path),
        "shared_core_index": artifact_identity(shared_path),
    }
    store = InvestigationStore.create(tmp_path / "with-corpus", objective="Check a route step", baseline=baseline)
    study = ScientificWorkspace(store.root)
    route, _ = call(study, "assess_route_proposal", proposal=_route_value())
    event, result = call(study, "assess_step_plausibility", source_ref=route.artifact_ref, step_id="step-1",
                         candidate_limit=8, match_limit=1)
    selected = result["proposal"]
    assert selected["external_step_id"] == result["step_id"] == "step-1"
    assert result["source_ref"] == route.artifact_ref
    assert event.evidence_refs == (route.artifact_ref,)
    engine = study.operations._conditions()
    expected = assess_step_plausibility(
        f"{selected['precursor_smiles']}>>{selected['target_smiles']}", engine.index, engine.shared_core_index,
        candidate_limit=8, match_limit=1,
    )
    assert canonical_bytes(result["assessment"]) == canonical_bytes(expected.to_dict())
    assert result["assessment"]["matches"][0]["precedent"]["reaction_id"] == "reaction-1"
    brief = study.call_summary(event)["result_summary"]
    detailed = study.call_summary(event, detailed=True)["result_summary"]["assessment"]
    assert brief["grade"] == detailed["grade"] == expected.grade
    assert brief["first_match"]["reaction_id"] == "reaction-1"
    assert detailed["matches"][0]["outcome"] == expected.matches[0].outcome
    source_field = ("result", "assessment", "matches", 0, "precedent", "reference_id")
    assert study.inspect_artifact(event.artifact_ref, source_field)["preview"] == brief["first_match"]["reference_id"]
    dossier = study.step_dossier(route.artifact_ref, step_id="step-1", evidence_refs=(event.artifact_ref,))
    attached, = store.read_artifact(dossier.artifact_ref)["assessments"]["step_plausibility"]
    assert attached["grade"] == expected.grade and attached["first_match"]["reaction_id"] == "reaction-1"
    replay = ScientificWorkspace(store.root).replay(event.artifact_ref)
    assert store.read_artifact(replay.artifact_ref)["matches"]

    # A supplied map must remain the query for this particular proposal.
    proposal = {"target_smiles": "CCN", "precursor_smiles": "CC=O.N",
                "mapped_reaction_smiles": "[CH3:1][CH:2]=[O:3].[NH3:4]>>[CH3:1][CH2:2][NH2:4]"}
    saved, structural = call(study, "assess_route_step", proposal=proposal)
    _, mapped = call(study, "assess_step_plausibility", source_ref=saved.artifact_ref)
    assert mapped["assessment"]["query_reaction_smiles"] == structural["assessment"]["normalized_mapped_reaction_smiles"]
    assert mapped["proposal"]["mapped_reaction_smiles"] == proposal["mapped_reaction_smiles"]

    # Direct assessment needs only the bound indexes and preserves maps and middle-field coordinate bonds.
    left, right = proposal["mapped_reaction_smiles"].split(">>")
    direct_query = f"{left}>N->[Cu+2]>{right}"
    direct_event, direct = call(study, "assess_step_plausibility", reaction_smiles=direct_query,
                                candidate_limit=8, match_limit=1)
    assert direct_event.evidence_refs == ()
    assert direct["source_ref"] is None and direct["step_id"] is None
    assert direct["proposal"] == {"precursor_smiles": left, "target_smiles": right}
    assert direct["assessment"]["query_reaction_smiles"] == direct_query
    direct_expected = assess_step_plausibility(direct_query, engine.index, engine.shared_core_index,
                                                candidate_limit=8, match_limit=1)
    assert canonical_bytes(direct["assessment"]) == canonical_bytes(direct_expected.to_dict())


def test_selector_and_mapping_conflicts_cannot_produce_another_steps_grade(workspace: ScientificWorkspace) -> None:
    route = _route_value()
    route["steps"][0]["external_step_id"] = "step-1"
    source, _ = call(workspace, "assess_route_proposal", proposal=route)
    for step_id, message in (("step-1", "ambiguous"), ("unknown", "not present")):
        with pytest.raises(ValueError, match=message):
            workspace.operations.assess_step_plausibility(source.artifact_ref, step_id=step_id)
    conflicting, _ = call(workspace, "assess_route_step", proposal={
        "target_smiles": "CCO", "precursor_smiles": "CC=O.N",
        "mapped_reaction_smiles": "[CH3:1][CH:2]=[O:3].[NH3:4]>>[CH3:1][CH2:2][NH2:4]",
    })
    failed = workspace.run("assess_step_plausibility", {"source_ref": conflicting.artifact_ref})
    failure = workspace.store.read_artifact(failed.artifact_ref)
    assert failure["execution_status"] == "error"
    assert "no validated mapping" in failure["error"]["message"]
    assert "result" not in failure


@pytest.mark.parametrize(("arguments", "error"), [
    ({}, "exactly one"),
    ({"source_ref": "saved", "reaction_smiles": "C>>C"}, "exactly one"),
    ({"reaction_smiles": "C>>C", "step_id": "step-1"}, "requires a saved source_ref"),
    ({"reaction_smiles": "C>>"}, "must contain precursors"),
])
def test_direct_plausibility_requires_one_unambiguous_input(tmp_path: Path, arguments: dict, error: str) -> None:
    store = InvestigationStore.create(tmp_path / "invalid-input", objective="Check tool input", baseline={})
    with pytest.raises(ValueError, match=error):
        ScientificWorkspace(store.root).operations.assess_step_plausibility(**arguments)
