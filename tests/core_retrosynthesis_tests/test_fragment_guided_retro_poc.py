"""Chemistry evidence and deterministic orchestration contracts for the POC."""

from copy import deepcopy
from dataclasses import asdict
import hashlib
import json
from types import SimpleNamespace

import pytest

from condition_recommender.fragment_index import build_fragment_index
from condition_recommender.fragment_search import search_fragment_precedents
from core_retrosynthesis import build_generic_library, save_generic_library
from core_retrosynthesis.fragment_guidance import (
    DEFAULT_POLICY, evaluate_transfers, load_policy, project_construction_guidance,
    restore_source_smiles, scientific_projection, summarize_transfers,
)
from examples.ai_native.fragment_guided_retrosynthesis_poc import render_report
from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.search_fragments import suggest_search_fragments


REACTION = (
    "[CH3:1][CH2:2][Br:5].[NH2:3][CH3:4]>>"
    "[CH3:1][CH2:2][NH:3][CH3:4]"
)


def search_evidence(tmp_path, reaction=REACTION, target="CCNC"):
    """Build a tiny canonical fixture through the production discovery contracts."""
    row = {"observation_id": "obs-1", "reaction_id": "source-1", "reference_id": "ref-1",
           "reaction_smiles": reaction, "admission_tier": "review",
           "reaction_observation": asdict(featurize_reaction(reaction).observation)}
    source, index = tmp_path / "source.jsonl", tmp_path / "fragments.sqlite"
    source.write_text(json.dumps(row) + "\n", "utf-8")
    build_fragment_index(source, index)
    suggestions = suggest_search_fragments(target).to_dict()
    candidate = next(c for c in suggestions["candidates"] if c["query"] == target)
    result = search_fragment_precedents(index, candidate["query"], target_smiles=target)
    searches = [{"candidate_id": candidate["candidate_id"], "artifact_ref": "search-1",
                 "execution_status": "completed", "result": result}]
    return row, suggestions, searches


@pytest.fixture
def evidence(tmp_path):
    return search_evidence(tmp_path)


def test_observed_construction_guides_admitted_source_transfer(tmp_path, monkeypatch):
    reaction = "[CH2:1]=[CH2:2].[CH2:3]=[CH2:4]>>[CH2:1]1[CH2:2][CH2:3][CH2:4]1"
    row, suggestions, searches = search_evidence(tmp_path, reaction, "C1CCC1")
    policy = load_policy(DEFAULT_POLICY)
    guidance = project_construction_guidance(suggestions, searches, policy)
    assert guidance["focus_bonds"]
    assert guidance["projection_basis"] == "query_embedding_alignment_not_reaction_mapping"
    library = build_generic_library([row], admission_mode="data_driven", levels=("L2", "L1", "L0"))
    path = tmp_path / "library.json"
    save_generic_library(library, path)
    monkeypatch.chdir(tmp_path)
    result = evaluate_transfers({"guidance": guidance, "policy": policy, "library": str(path)})
    candidates = result["guided"][0]["direct_source_transfer"]["candidates"]
    assert any(c["precursor_smiles"] == "C=C.C=C" for c in candidates)
    assert all(c["forward_validation_status"] == "verified_signature" for c in candidates)
    assert all(c["bond_focus_check"]["status"] == "verified" for c in candidates)
    assert result["repeat_scientific_results_identical"]


def test_local_witness_does_not_override_whole_source_compiler_rejection(evidence, tmp_path, monkeypatch):
    row, suggestions, searches = evidence
    policy = load_policy(DEFAULT_POLICY)
    guidance = project_construction_guidance(suggestions, searches, policy)
    assert guidance["focus_bonds"][0]["target_atom_ids"] == [1, 2]
    # A separate known source supplies the baseline. The mapped leaving-atom
    # fixture still has a local witness but fails the current whole-core gate.
    known = {**row, "reaction_smiles": REACTION.replace("[Br:5]", "Br")}
    library = build_generic_library([known], admission_mode="data_driven", levels=("L2", "L1", "L0"))
    path = tmp_path / "library.json"
    save_generic_library(library, path)
    monkeypatch.chdir(tmp_path)
    result = evaluate_transfers({"guidance": guidance, "policy": policy, "library": str(path)})
    assert result["compiled_source_template_count"] == 0
    assert result["source_admissions"][0]["reason"] == "materialized_core_not_verified"
    assert result["guided"][0]["direct_source_transfer"]["candidates"] == []
    assert result["guided"][0]["witness_directed_library"]["candidates"]


@pytest.mark.parametrize("case", ["partial", "timed_out", "mismatch", "unresolved", "truncated", "retention"])
def test_uncertain_or_nonconstruction_evidence_cannot_seed_a_transfer(evidence, case):
    _, suggestions, searches = deepcopy(evidence)
    search, result = searches[0], searches[0]["result"]
    if case == "partial":
        result["search_status"] = "partial"
    elif case == "timed_out":
        search["execution_status"] = "timed_out"
    elif case == "mismatch":
        result["target_validation"]["target_smiles"] = "CCO"
    else:
        for hit in result["hits"]:
            for match in hit["matches"]:
                if case == "unresolved":
                    match["relationships"].append("unresolved")
                elif case == "truncated":
                    match["embedding_truncated"] = True
                else:
                    match["relationships"] = ["carried_through"]
                    match["witnesses"] = []
    guidance = project_construction_guidance(suggestions, searches, load_policy(DEFAULT_POLICY))
    assert guidance["focus_bonds"] == []
    assert guidance["source_records"] == []


def test_ambiguous_embeddings_remain_separate_hypotheses(tmp_path):
    reaction = "[CH2:1]=[CH2:2].[CH2:3]=[CH2:4]>>[CH2:1]1[CH2:2][CH2:3][CH2:4]1"
    _, suggestions, searches = search_evidence(tmp_path, reaction, "C1CCC1")
    hit = searches[0]["result"]["hits"][0]
    alternate = deepcopy(hit["matches"][0])
    # Rotation is a valid alternate query embedding in a symmetric four-ring.
    original = alternate["query_to_original_product_atoms"]
    alternate["query_to_original_product_atoms"] = original[1:] + original[:1]
    hit["matches"].append(alternate)
    guidance = project_construction_guidance(suggestions, searches, load_policy(DEFAULT_POLICY))
    assert guidance["eligible_source_count"] == 1
    assert guidance["eligible_bond_count"] == 4
    assert {s["match_index"] for b in guidance["focus_bonds"] for s in b["supports"]} == {0, 1}


def test_conflicting_projected_bond_is_an_error(evidence):
    _, suggestions, searches = deepcopy(evidence)
    candidate = next(c for c in suggestions["candidates"]
                     if c["candidate_id"] == searches[0]["candidate_id"])
    candidate["query_atom_target_ids"] = [0, 1, 3, 2]
    with pytest.raises(ValueError, match="not a target bond"):
        project_construction_guidance(suggestions, searches, load_policy(DEFAULT_POLICY))


def test_policy_cannot_disable_evidence_gates(tmp_path):
    value = load_policy(DEFAULT_POLICY)
    value["require_resolved_embedding"] = False
    path = tmp_path / "policy.json"
    path.write_text(json.dumps(value), "utf-8")
    with pytest.raises(ValueError, match="resolved embeddings"):
        load_policy(path)


def test_repeat_comparison_ignores_only_execution_measurements():
    left = {"search_status": "complete", "hits": [], "execution": {"elapsed_seconds": 1}}
    right = {**left, "execution": {"elapsed_seconds": 2}}
    assert scientific_projection(left) == scientific_projection(right)
    right["search_status"] = "partial"
    assert scientific_projection(left) != scientific_projection(right)


def test_chunked_source_structures_are_losslessly_restored(evidence):
    _, suggestions, searches = deepcopy(evidence)
    text = searches[0]["result"]["hits"][0]["record"]["reaction_smiles"]
    chunks = {"truncated": False, "total_characters": len(text),
              "source_sha256": hashlib.sha256(text.encode()).hexdigest(),
              "chunks": [{"start": 0, "end": 30, "text": text[:30]},
                         {"start": 30, "end": len(text), "text": text[30:]}]}
    assert restore_source_smiles(chunks) == text
    searches[0]["result"]["hits"][0]["record"]["reaction_smiles"] = chunks
    guidance = project_construction_guidance(suggestions, searches, load_policy(DEFAULT_POLICY))
    assert guidance["source_records"][0]["reaction_smiles"] == text
    chunks["source_sha256"] = "wrong"
    with pytest.raises(ValueError, match="hash conflicts"):
        project_construction_guidance(suggestions, searches, load_policy(DEFAULT_POLICY))


def test_truncated_source_structure_is_not_compiled(evidence):
    _, suggestions, searches = deepcopy(evidence)
    searches[0]["result"]["hits"][0]["record"]["reaction_smiles"] = {"truncated": True}
    guidance = project_construction_guidance(suggestions, searches, load_policy(DEFAULT_POLICY))
    assert not guidance["source_records"]
    assert guidance["exclusions"][0]["reason"] == "missing_or_truncated_source_structure"


def test_review_report_draws_proposals_and_source_records(tmp_path):
    render_report({"target_smiles": "CCNC", "guidance": {"source_records": [
        {"reaction_id": "source-1", "reaction_smiles": REACTION},
    ]}, "transfers": {"baseline": {"candidates": [{
        "proposed_reaction_smiles": "CCBr.CN>>CCNC", "forward_validation_status": "verified_signature",
        "abstraction_level": "L2", "precedent_reaction_ids": ["source-1"],
    }]}, "guided": [], "source_admissions": []}, "limitations": ["Development only"]}, tmp_path)
    page = (tmp_path / "report.html").read_text("utf-8")
    assert page.count("<svg") == 2
    assert "Development only" in page and "source-1" in page


@pytest.mark.parametrize("signature, focus", [
    ("unresolved", {"status": "verified"}),
    ("verified_signature", None),
    ("verified_signature", {"status": "rejected"}),
])
def test_transfer_results_cannot_bypass_signature_or_focused_formation_gate(monkeypatch, signature, focus):
    import core_retrosynthesis.generic_search as search

    def invalid_result(*args, **kwargs):
        candidate = {"forward_validation_status": signature, "bond_focus_check": focus}
        if "required_disconnection_bond" not in kwargs:
            candidate = {"forward_validation_status": "verified_signature", "bond_focus_check": None}
        return [SimpleNamespace(to_dict=lambda: candidate)], SimpleNamespace(to_dict=lambda: {})

    monkeypatch.setattr(search, "disconnect_operator_ladder_detailed", invalid_result)
    parameters = {
        "guidance": {"target_smiles": "CC", "source_records": [],
                     "focus_bonds": [{"target_atom_ids": [0, 1]}]},
        "policy": load_policy(DEFAULT_POLICY),
    }
    with pytest.raises(ValueError, match="signature|formation validation"):
        evaluate_transfers(parameters, library=SimpleNamespace(templates=[1]),
                           repeat=False, source_library_path=None)


def test_coverage_counts_deduplicate_precursor_sets_and_sum_actual_work():
    def arm(precursors, templates, validations):
        return {"candidates": [{"precursor_smiles": s} for s in precursors],
                "diagnostics": {"validation_attempt_count": validations,
                                "level_diagnostics": {"L2": {"applied_template_count": templates}}}}
    result = summarize_transfers({
        "baseline": arm(["CC", "CC"], 10, 3),
        "guided": [{"direct_source_transfer": arm(["CC", "CO"], 5, 2),
                    "witness_directed_library": arm(["CO", "CN"], 20, 4)}],
    })
    assert result["baseline_unique_precursor_count"] == 1
    assert result["guided_unique_precursor_count"] == 3
    assert result["additional_guided_precursor_sets"] == ["CN", "CO"]
    assert result["guided_work"] == {"template_applications": 25, "validation_attempts": 6}
