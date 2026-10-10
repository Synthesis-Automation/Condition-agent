"""Automatic discovery retains construction evidence through query refinement."""

from dataclasses import asdict
from copy import deepcopy
import json

import pytest

from condition_recommender.fragment_index import build_fragment_index
from condition_recommender.precedent_discovery import find_synthesis_precedents
from reactive_taxonomy import featurize_reaction


@pytest.fixture
def discovery_index(tmp_path):
    # Synthetic graph fixtures test observation/ranking, not experimental feasibility.
    ring = "[CH2:1]1[CH2:2][NH:3][CH2:4][CH2:5]1"
    reactions = [
        ("construction", "ref-build", "[CH3:1][CH2:2][NH:3][CH2:4][CH3:5]>>" + ring),
        ("retention", "ref-acyl", ring + ".[CH3:6][C:7](=[O:8])[Cl:9]>>[CH2:1]1[CH2:2][N:3]([C:7]([CH3:6])=[O:8])[CH2:4][CH2:5]1"),
        ("unknown", "ref-unknown", "NCCCC>>C1CCNC1"),
    ]
    source = tmp_path / "source.jsonl"
    source.write_text("\n".join(json.dumps({"observation_id": oid, "reaction_id": oid,
        "reference_id": ref, "reaction_smiles": reaction, "admission_tier": "review",
        "reaction_observation": asdict(featurize_reaction(reaction).observation)})
        for oid, ref, reaction in reactions), encoding="utf-8")
    index = tmp_path / "fragments.sqlite"
    build_fragment_index(source, index)
    return index


def test_construction_survives_acetyl_refinement_and_library_is_loaded_once(discovery_index):
    result = find_synthesis_precedents(discovery_index, "CC(=O)N1CCCC1")
    assert result["search_status"] == "complete"
    assert result["execution"]["library_loads"] == 1
    assert result["execution"]["candidate_reuses"] > 0
    assert result["hits"][0]["observation_id"] == "construction"
    assert result["hits"][0]["discovery"]["core_relationship"] == "constructed"
    assert result["hits"][0]["discovery"]["nitrogen_hydrogen_differences"]
    retention = next(h for h in result["hits"] if h["observation_id"] == "retention")
    assert retention["discovery"]["core_relationship"] != "constructed"
    unknown = next(h for h in result["hits"] if h["observation_id"] == "unknown")
    assert unknown["discovery"]["core_relationship"] == "unresolved"
    assert result["exact_target"]["status"] == "matched"


def test_no_matches_stops_a_ladder_without_claiming_literature_absence(discovery_index):
    result = find_synthesis_precedents(discovery_index, "c1ccncc1")
    assert result["hits"] == []
    assert len(result["attempts"]) == 1
    assert result["attempts"][0]["decision"] == "no_matches_try_next_core"
    assert any("absence" in text for text in result["limitations"])


def test_discovery_hit_inspection_preserves_source_and_partial_scope(discovery_index):
    from core_retrosynthesis.fragment_investigation import investigate_fragment_precedent
    result = find_synthesis_precedents(discovery_index, "CC(=O)N1CCCC1")
    original = deepcopy(result)
    inspected = investigate_fragment_precedent(result, "construction", result["target_smiles"])
    assert inspected["source"] == result["hits"][0]
    assert inspected["query"]["expression"] == result["hits"][0]["discovery"]["query"]
    assert result == original
    result["search_status"] = "partial"
    result["stop_reason"] = "deadline"
    inspected = investigate_fragment_precedent(result, "construction", result["target_smiles"])
    assert inspected["transfer"]["status"] == "incomplete_search"
    assert inspected["search_scope"]["search_status"] == "partial"


@pytest.mark.parametrize("mutation", ["forged_id", "missing_attempt", "ambiguous_attempt", "wrong_query", "wrong_target"])
def test_discovery_inspection_rejects_unbound_evidence(discovery_index, mutation):
    from condition_recommender.fragment_investigation import inspect_fragment_precedent
    result = find_synthesis_precedents(discovery_index, "CC(=O)N1CCCC1")
    target, observation = result["target_smiles"], "construction"
    if mutation == "forged_id":
        observation = "not-returned"
    elif mutation == "missing_attempt":
        result["attempts"] = []
    elif mutation == "ambiguous_attempt":
        result["attempts"] *= 2
    elif mutation == "wrong_query":
        result["hits"][0]["discovery"]["query"] = "[#6]"
    else:
        target = "c1ccccc1"
    with pytest.raises(ValueError):
        inspect_fragment_precedent(result, observation, target)


def test_partial_candidates_are_not_reused(discovery_index, monkeypatch):
    import condition_recommender.fragment_search as search
    original = search.fragment_search_policy
    monkeypatch.setattr(search, "fragment_search_policy", lambda: {**original(), "max_matched_products": 1})
    result = find_synthesis_precedents(discovery_index, "CC(=O)N1CCCC1")
    assert result["attempts"][0]["search_status"] == "too_broad"
    assert not result["attempts"][1]["execution"]["candidate_set_reused"]


def test_complete_candidate_reuse_matches_fresh_search(discovery_index):
    from condition_recommender.fragment_search import FragmentSearchSession, search_fragment_precedents
    from reactive_taxonomy.precedent_queries import plan_precedent_queries
    steps = plan_precedent_queries("CC(=O)N1CCCC1")["ladders"][0]
    with FragmentSearchSession(discovery_index) as session:
        parent = search_fragment_precedents(discovery_index, steps[0]["query"], "smarts", "subgraph", _session=session)
        ids = session.candidates[parent["query"]["query_id"]]
        reused = search_fragment_precedents(discovery_index, steps[-1]["query"], "smarts", "subgraph", _session=session, _candidate_ids=ids)
    fresh = search_fragment_precedents(discovery_index, steps[-1]["query"], "smarts", "subgraph")
    assert reused["counts"] == fresh["counts"]
    assert reused["hits"] == fresh["hits"]


@pytest.mark.parametrize("include_ester", [True, False])
@pytest.mark.parametrize("product_cap", [1, 500])
def test_sibling_contexts_do_not_prune_each_other(tmp_path, include_ester, product_cap, monkeypatch):
    import condition_recommender.fragment_search as search
    from condition_recommender.fragment_search import search_fragment_precedents

    original_policy = search.fragment_search_policy
    monkeypatch.setattr(search, "fragment_search_policy", lambda: {
        **original_policy(), "max_matched_products": product_cap,
    })

    products = ["COc1ccccn1"]
    if include_ester:
        products.append("CCOC(=O)c1ccncc1")
    source = tmp_path / "source.jsonl"
    source.write_text("\n".join(json.dumps({
        "observation_id": str(i), "reaction_id": str(i), "reference_id": str(i),
        "reaction_smiles": "CC>>" + product, "admission_tier": "review",
        "reaction_observation": asdict(featurize_reaction("CC>>" + product).observation),
    }) for i, product in enumerate(products)), encoding="utf-8")
    index = tmp_path / "fragments.sqlite"
    build_fragment_index(source, index)
    result = find_synthesis_precedents(index, "COc1cc(C(=O)OCC)ccn1")
    focused = [a for a in result["attempts"] if a["level"] == "focused_context"]
    assert len(focused) == 2
    assert focused[0]["counts"]["products"]["value"] == int(include_ester)
    assert focused[1]["counts"]["products"]["value"] == 1
    if not include_ester:
        assert focused[0]["decision"] == "no_matches_try_next_context"
    for attempt in focused:
        fresh = search_fragment_precedents(index, attempt["query"], "smarts", "subgraph")
        assert attempt["counts"] == fresh["counts"]
        assert attempt["returned_observation_ids"] == [hit["observation_id"] for hit in fresh["hits"]]
    assert result["execution"]["library_loads"] == 1
    assert result["search_status"] == ("partial" if include_ester and product_cap == 1 else "complete")
    assert all(hit["discovery"]["core_relationship"] == "unresolved" for hit in result["hits"])


def test_focused_schedule_preserves_deadline_and_saved_hits(discovery_index, monkeypatch):
    from itertools import chain, repeat
    import condition_recommender.precedent_discovery as discovery

    clock = chain([0.0, 0.0, 0.0], repeat(2.0))
    monkeypatch.setattr(discovery, "monotonic", lambda: next(clock))
    result = discovery.find_synthesis_precedents(discovery_index, "CC(=O)N1CCCC1", timeout_seconds=1)
    assert len(result["attempts"]) == 1
    assert result["hits"]
    assert result["search_status"] == "partial"
    assert result["stop_reason"] == "deadline"
