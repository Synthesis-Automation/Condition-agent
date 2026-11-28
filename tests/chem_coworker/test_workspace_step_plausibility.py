"""Corpus support is recorded for the exact saved route step the agent selects."""

from pathlib import Path

import pytest

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.core.baseline import artifact_identity
from chem_coworker.scientific_workspace.core.store import canonical_bytes
from chem_coworker.scientific_workspace.views.call_summaries import summarize_call
from chem_coworker.scientific_workspace.views.brief_summaries import summarize_call_brief
from condition_recommender import assess_step_plausibility
from condition_recommender.generic_indexing import build_generic_index
from condition_recommender.shared_core_index import build_shared_core_index
from condition_recommender.sqlite_indexing import load_sqlite_generic_index, save_sqlite_generic_index
from core_retrosynthesis.chemistry import canonical_smiles
from tests.chem_coworker import test_route_investigation_workspace as route_tests
from tests.chem_coworker.test_route_investigation_workspace import call
from tests.condition_recommender.test_shared_core_retrieval import record
from tests.core_retrosynthesis_tests.test_external_proposal_admission import _route_value


library = route_tests.library
workspace = route_tests.workspace


@pytest.mark.parametrize(("oxime", "status"), (
    ("[CH2:3][C:4](=[N:5][OH:6])", "unspecified"),
    (r"[CH2:3]/[C:4](=[N:5]\[OH:6])", "anti_aligned"),
    ("[CH2:3]/[C:4](=[N:5]/[OH:6])", "syn_aligned"),
))
def test_step_and_route_surface_migration_geometry(workspace, oxime, status):
    ring = "[c:7]2[cH:8][cH:9][cH:10][cH:11][c:12]21"
    precursor = "[CH2:1]1[CH2:2]" + oxime + ring
    product = "[CH2:1]1[CH2:2][CH2:3][C:4](=[O:6])[NH:5]" + ring
    step = {"precursor_smiles": canonical_smiles(precursor), "target_smiles": canonical_smiles(product),
            "mapped_reaction_smiles": precursor + ">>" + product}
    _, single = call(workspace, "assess_route_step", proposal=step)
    event, route = call(workspace, "assess_route_proposal", proposal={
        "target_smiles": step["target_smiles"], "steps": [{"external_step_id": "closure", **step}],
    })
    check = single["step_review"]["compatibility"]
    assert route["step_reviews"][0]["compatibility"] == check
    geometry = check["checks"]["migration_geometry"]
    assert geometry["status"] == status
    assert geometry["alternative_product_smiles"] == "O=C1NCCCc2ccccc21"
    assert check["checks"]["functional_group_selectivity"]["status"] == "out_of_scope"
    flagged = status in {"unspecified", "syn_aligned"}
    assert any(alert.get("code") == geometry["audit_id"] for alert in check["alerts"]) == flagged
    assert ("closure" in route["weakest_steps"]["compatibility_alert_step_ids"]) == flagged
    summary = workspace.call_summary(event)["result_summary"]["step_reviews"][0]
    assert summary["compatibility"]["checks"]["migration_geometry"] == geometry
    validity_event, _ = call(workspace, "assess_retro_validity", source_ref=event.artifact_ref, step_id="closure")
    payload = workspace.store.read_artifact(validity_event.artifact_ref)
    for project in (summarize_call, summarize_call_brief):
        assert project(payload)["result_summary"]["validity"]["migration_geometry"] == geometry


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
    route, route_result = call(study, "assess_route_proposal", proposal=_route_value())
    reviews = {item["step_id"]: item for item in route_result["step_reviews"]}
    supported = reviews["step-1"]
    assert supported["recipe_status"] == "not_supplied"
    assert supported["precedent"]["status"] == "assessed"
    assert supported["precedent"]["closest_match"]["reference_id"] == "REF1:1"
    assert supported["precedent"]["closest_match"]["level"] == "whole_reaction"
    assert reviews["step-2"]["precedent"]["evidence_rank"] < supported["precedent"]["evidence_rank"]
    assert route_result["weakest_steps"]["step_ids"] == ["step-2"]
    _, inspected = call(study, "inspect_route_step", source_ref=route.artifact_ref, step_id="step-1")
    assert inspected["step_review"] == supported
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
    assert structural["step_review"]["precedent"]["query_reaction_smiles"] == mapped["assessment"]["query_reaction_smiles"]
    assert structural["step_review"]["precedent"]["query_basis"] == "validated_supplied_mapping"
    assert structural["step_review"]["precedent"]["grade"] == mapped["assessment"]["grade"]
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
    conflicting, conflict_result = call(workspace, "assess_route_step", proposal={
        "target_smiles": "CCO", "precursor_smiles": "CC=O.N",
        "mapped_reaction_smiles": "[CH3:1][CH:2]=[O:3].[NH3:4]>>[CH3:1][CH2:2][NH2:4]",
    })
    assert conflict_result["step_review"]["precedent"]["status"] == "unresolved_structure"
    assert "reaction_side_consistency" in conflict_result["step_review"]["precedent"]["unverified_gates"]
    failed = workspace.run("assess_step_plausibility", {"source_ref": conflicting.artifact_ref})
    failure = workspace.store.read_artifact(failed.artifact_ref)
    assert failure["execution_status"] == "error"
    assert "no validated mapping" in failure["error"]["message"]
    assert "result" not in failure


def test_absent_optional_corpus_preserves_seen_gates_and_unknown_support(workspace):
    _, result = call(workspace, "assess_route_proposal", proposal=_route_value())
    weak = result["weakest_steps"]
    assert weak["step_ids"] == [] and weak["basis"] == "insufficient_evidence"
    assert set(weak["unassessed_precedent_step_ids"]) == {"step-1", "step-2"}
    for review in result["step_reviews"]:
        assert review["precedent"] == {"status": "not_configured", "missing_indexes": ["condition_index", "shared_core_index"]}
        assert review["compatibility"]["checks"]["reaction_compatibility"]["status"] == "pass"


@pytest.mark.parametrize("paired,expected", [(False, "incomplete"), (True, "unavailable")])
def test_optional_corpus_failure_preserves_step_review(workspace, tmp_path, monkeypatch, paired, expected):
    dummy = tmp_path / "invalid-index.sqlite"
    dummy.write_bytes(b"incompatible development fixture")
    baseline = dict(workspace.store.manifest["baseline"])
    baseline["artifacts"] = {**baseline["artifacts"], "condition_index": artifact_identity(dummy)}
    if paired:
        baseline["artifacts"]["shared_core_index"] = artifact_identity(dummy)
    store = InvestigationStore.create(tmp_path / "incompatible", objective="Retain structural review", baseline=baseline)
    study = ScientificWorkspace(store.root)

    def refuse_index():
        assert paired
        raise ValueError("Incompatible indices")

    monkeypatch.setattr(study.operations, "_conditions", refuse_index)
    _, result = call(study, "assess_route_step", proposal={"precursor_smiles": "CC=O.N", "target_smiles": "CCN"})
    assert result["step_review"]["precedent"]["status"] == expected
    assert result["assessment"]["canonical_target_smiles"] == "CCN"


def test_review_distinguishes_conflicts_ties_unassessed_and_ambiguous_steps():
    from chem_coworker.scientific_workspace.adapters.route_step_review import weakest_steps

    first = {"step_id": "make-amine", "structural_status": "verified_with_cautions", "recipe_conflicts": [],
             "precedent": {"evidence_rank": 2}, "compatibility": {"alerts": [{"rule_id": "advisory"}]}}
    second = {**first, "step_id": "make-amide", "compatibility": {"alerts": []}}
    unknown = {**second, "step_id": "unsearched", "precedent": {"status": "unavailable"}}
    tied = weakest_steps([first, second, unknown])
    assert tied["step_ids"] == ["make-amine", "make-amide"]
    assert tied["basis"] == "lowest_assessed_precedent_grade"
    assert tied["unassessed_precedent_step_ids"] == ["unsearched"]
    assert tied["compatibility_alert_step_ids"] == ["make-amine"]
    conflict = weakest_steps([first, {**second, "recipe_conflicts": ["documented-rule"]}])
    assert conflict["basis"] == "structural_or_recipe_conflict" and conflict["step_ids"] == ["make-amide"]
    ambiguous = weakest_steps([first, first])
    assert ambiguous["step_ids"] == [] and ambiguous["ambiguous_step_ids"] == ["make-amine"]


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
