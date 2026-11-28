"""Attribution and rendering of scoped step evidence, without a corpus search."""

from copy import deepcopy

import pytest

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.answers.answer_contracts import ScientificAnswer, validate_answer_evidence
from chem_coworker.scientific_workspace.views.step_dossier import step_plausibility_summary


@pytest.fixture
def workspace(tmp_path):
    store = InvestigationStore.create(tmp_path / "study", objective="Inspect a route step", baseline={})
    return ScientificWorkspace(store.root)


def save(workspace, operation, result, *, arguments=None):
    return workspace.store.append("call", {
        "operation": operation, "execution_status": "completed", "arguments": arguments or {}, "result": result,
    }).artifact_ref


def plausibility(source_ref, proposal):
    """An assessment-shaped display fixture with an adverse source experiment."""
    return {
        "schema_version": "route_step_plausibility.v1", "source_ref": source_ref,
        "step_id": "s1", "proposal": proposal,
        "assessment": {
            "query_reaction_smiles": f"{proposal['precursor_smiles']}>>{proposal['target_smiles']}",
            "grade": "context_supported", "evidence_rank": 2, "measured_evidence_status": "only_low_or_zero",
            "qualified_count": 2, "candidate_count": 4, "candidate_limit": 4,
            "candidate_truncated": True, "matches_truncated": True, "warnings": ["BOUNDED_SEARCH"],
            "support_counts": {"supporting_source_group_count": 0, "comparable_observation_count": 1,
                               "zero_yield_count": 1, "low_yield_count": 0, "unreported_yield_count": 0,
                               "unknown_yield_count": 0},
            "matches": [{
                "precedent": {"reaction_id": "reaction-a", "reference_id": "paper", "observation_id": "well",
                              "reaction_smiles": "CCBr.N>>CCN", "level": "whole_reaction",
                              "yield_pct": 0.0, "differences": []},
                "context": {"status": "covered", "reverse_status": "covered", "minimum_partner_similarity": 1.0,
                            "missing": [], "unresolved": [], "source_only": [], "source_unresolved": []},
                "outcome": "zero_yield", "contributes_to_support": False,
            }],
        },
    }


def test_cited_plausibility_is_bound_to_exact_source_and_route_step(workspace):
    proposal = {"precursor_smiles": "CCBr.N", "target_smiles": "CCN", "external_step_id": "s1"}
    steps = [proposal, {**proposal, "external_step_id": "s2"}]
    route = {"proposal": {"steps": steps}, "assessment": {"step_assessments": [
        {"external_step_id": item["external_step_id"], "assessment": {}} for item in steps
    ]}}
    source = save(workspace, "assess_route_proposal", route)
    revised = save(workspace, "revise_route_branch", route)
    ref = save(workspace, "assess_step_plausibility", plausibility(source, proposal))
    for parent, step_id, links in ((source, "s1", ()), (source, "s2", (ref,)), (revised, "s1", (ref,))):
        event = workspace.step_dossier(parent, step_id=step_id, evidence_refs=links)
        assert not workspace.store.read_artifact(event.artifact_ref)["assessments"]["step_plausibility"]

    event = workspace.step_dossier(source, step_id="s1", evidence_refs=(ref,))
    dossier = workspace.store.read_artifact(event.artifact_ref)
    [check] = dossier["assessments"]["step_plausibility"]
    assert check["artifact_ref"] == ref
    assert check["grade"] == "context_supported"
    assert check["support_counts"]["supporting_source_group_count"] == 0
    assert check["support_counts"]["zero_yield_count"] == 1
    assert check["candidate_truncated"] and check["matches_truncated"]
    assert check["first_match"]["source_yield_pct"] == 0
    assert not check["first_match"]["contributes_to_support"]
    assert dossier["experimental_feasibility"] == "not_established"


@pytest.mark.parametrize(("grade", "rank"), [("unresolved", None), ("no_qualified_precedent", 0)])
def test_empty_plausibility_keeps_the_engine_grade(grade, rank):
    result = plausibility("source", {"precursor_smiles": "CCBr.N", "target_smiles": "CCN"})
    result["assessment"].update(grade=grade, evidence_rank=rank, matches=[], qualified_count=0)
    summary = step_plausibility_summary(result)
    assert (summary["grade"], summary["evidence_rank"], summary["first_match"]) == (grade, rank, None)


def test_plausibility_citation_requires_exact_proposed_step(workspace):
    source = save(workspace, "assess_step_plausibility", plausibility("source", {
        "precursor_smiles": "CCBr.N", "target_smiles": "CCN",
    }))
    draft = workspace.route_answer({"target_smiles": "CCN", "routes": [{"steps": [{
        "reaction_smiles": "N.CCBr>>NCC", "support": [{"ref": source, "locator": "assessment"}],
    }]}]})
    assert source in validate_answer_evidence(ScientificAnswer.model_validate(draft), workspace.store)
    unrelated = deepcopy(draft)
    product_id = unrelated["steps"][0]["product_ids"][0]
    next(mol for mol in unrelated["molecules"] if mol["id"] == product_id)["smiles"] = "C[13CH2]N"
    with pytest.raises(ValueError, match="this exact reaction"):
        validate_answer_evidence(ScientificAnswer.model_validate(unrelated), workspace.store)
    draft["steps"][0]["basis"] = "reported"
    with pytest.raises(ValueError, match="reported step needs source evidence"):
        validate_answer_evidence(ScientificAnswer.model_validate(draft), workspace.store)


def test_dative_reaction_fields_survive_route_authoring_recipe_and_evidence_links(workspace):
    left, right = "N->[Ru].C", "CN->[Ru]"
    query = f"{left}>O>{right}"
    recipe = save(workspace, "assess_proposed_recipe", {
        "reaction_smiles": query, "proposed_recipe": {"recipe_id": "recipe"},
        "compatibility": {"status": "unknown"},
    }, arguments={"reaction_smiles": query})
    condition = save(workspace, "inspect_condition_precedents", {
        "schema_version": "condition_evidence_comparison.v2", "query_reaction_smiles": query,
    }, arguments={"reaction_smiles": query})
    draft = workspace.route_answer({"target_smiles": right, "routes": [{"steps": [{
        "reaction_smiles": f"{left}>>{right}",
        "support": [{"ref": ref, "locator": "recorded reaction"} for ref in (recipe, condition)],
    }]}]})
    assert set(validate_answer_evidence(ScientificAnswer.model_validate(draft), workspace.store)) == {recipe, condition}
    event = workspace.step_dossier(recipe)
    dossier = workspace.store.read_artifact(event.artifact_ref)
    assert dossier["proposal"] == {"precursor_smiles": left, "target_smiles": right}
    assert dossier["assessments"]["recipe_assessments"][0]["artifact_ref"] == recipe
