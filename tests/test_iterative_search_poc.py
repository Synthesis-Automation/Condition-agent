"""Pilot lineage, evidence preservation, and budget checks on a tiny real index."""

from dataclasses import asdict
import json

import pytest

from condition_recommender.fragment_index import build_fragment_index
from examples.ai_native.iterative_search_poc import SearchPilot, search_card
from examples.ai_native.iterative_search_replay import scientific_projection
from reactive_taxonomy import featurize_reaction


@pytest.fixture
def pilot(tmp_path):
    reaction = "[CH3:1][Br:2].[OH:3][CH3:4]>>[CH3:1][O:3][CH3:4]"
    source = tmp_path / "records.jsonl"
    source.write_text(json.dumps({
        "observation_id": "obs", "reaction_id": "rxn", "reference_id": "ref",
        "reaction_smiles": reaction,
        "reaction_observation": asdict(featurize_reaction(reaction).observation),
    }) + "\n", "utf-8")
    index = tmp_path / "index.sqlite"
    build_fragment_index(source, index)
    return SearchPilot.create(tmp_path / "pilot", {
        "objective": "Test recorded query selection", "target_smiles": "COC",
        "artifacts": {"fragment_index": str(index)}, "budget_seconds": 90,
    })


def test_saved_query_search_and_inspection_preserve_evidence(pilot):
    prepared = pilot.dispatch({"action": "prepare", "query": "COC", "reason": "Find ether formation"})
    result = pilot.dispatch({"action": "search", "query_ref": prepared["query_ref"],
                             "reason": "Search validated core"})
    assert result["cards"][0]["observation_id"] == "obs"
    assert "constructed" in result["cards"][0]["relationships"]
    opened = pilot.dispatch({"action": "inspect", "search_ref": result["search_ref"],
                             "observation_id": "obs", "reason": "Check actual formed bond"})
    assert opened["hit"]["record"]["reaction_smiles"].startswith("[CH3:1]")
    assert pilot.records("pilot_inspection")[0]["usefulness"] == "not_yet_adjudicated"
    child = pilot.dispatch({"action": "prepare", "query": "CO", "reason": "Broaden after inspection",
                            "parent_search_ref": result["search_ref"]})
    assert child["payload"]["execution_status"] == "completed"
    event = next(e for e in pilot.workspace.store.events() if e.artifact_ref == child["query_ref"])
    assert result["search_ref"] in event.evidence_refs


def test_invalid_target_fragment_cannot_be_searched(pilot):
    prepared = pilot.dispatch({"action": "prepare", "query": "c1ccccc1", "reason": "Negative control"})
    assert prepared["payload"]["execution_status"] == "error"
    with pytest.raises(ValueError, match="completed"):
        pilot.dispatch({"action": "search", "query_ref": prepared["query_ref"], "reason": "Must reject"})
    assert pilot.records("pilot_run") == []


def test_budget_refuses_more_searches(pilot):
    pilot.workspace.store.append("pilot_run", {"arm": "iterative", "operation_seconds": 90})
    with pytest.raises(ValueError, match="budget exhausted"):
        pilot.dispatch({"action": "search", "query_ref": "unused", "reason": "Must reject"})


def test_compact_cards_do_not_upgrade_unresolved_or_hide_warnings():
    card = search_card({"observation_id": "unmapped", "relationships": ["unresolved"],
                        "warnings": ["mapping_missing"], "matches": []})
    assert card["relationships"] == ["unresolved"]
    assert card["warnings"] == ["mapping_missing"]
    assert card["witnesses"] == []


def test_replay_comparison_preserves_scientific_status_and_evidence():
    result = {"execution": {"elapsed_seconds": 3}, "execution_status": "completed",
              "search_status": "partial", "counts": {"precision": "at_least"},
              "hits": [{"warnings": ["mapping_conflict"], "matches": []}]}
    projected = scientific_projection(result)
    assert projected["search_status"] == "partial"
    assert projected["counts"]["precision"] == "at_least"
    assert projected["hits"][0]["warnings"] == ["mapping_conflict"]
    assert "execution" not in projected
