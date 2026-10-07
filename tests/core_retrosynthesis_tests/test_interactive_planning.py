"""Regression checks for occurrence-preserving human route edits."""

from copy import deepcopy

import pytest

from core_retrosynthesis.interactive_planning import (
    PlanningChoice,
    edit_session,
    export_route,
    normalize_settings,
    record_search,
    restore_session,
    session_response,
    start_session,
)
from core_retrosynthesis.route_contract import (
    ReactionRouteTree,
    assert_valid_route_tree,
)


def search_result(target, *precursor_sets):
    """Build reported proposals; these tests assess topology, not their chemistry."""
    candidates = [
        dict(
            target_smiles=target,
            precursor_smiles=precursors,
            proposed_reaction_smiles=f"{precursors}>>{target}",
            forward_validation_status="verified_signature",
            score=0.8,
            template_id=f"template-{index}",
            operator_id=f"operator-{index}",
            disconnection_site_key="site",
            independent_reference_support=2,
            condition_precedent_reaction_ids=["source-1"],
        )
        for index, precursors in enumerate(precursor_sets)
    ]
    return dict(
        target_smiles=target,
        strategies=[
            dict(
                representative=candidate,
                alternate_realizations=[],
                strategy_id=f"strategy-{index}",
            )
            for index, candidate in enumerate(candidates)
        ],
        warnings=[],
        schema_version="2.0",
    )


def searched(session, node_id, *precursors):
    from core_retrosynthesis.interactive_planning import find_node

    target = find_node(session.root, node_id).smiles
    return record_search(session, node_id, {}, search_result(target, *precursors))


def choose(session, node_id="root", index=0):
    return edit_session(
        session,
        "select",
        node_id,
        PlanningChoice(session.searches[-1].search_id, index, 0),
    )


def test_branch_replacement_undo_redo_and_canonical_export():
    plan = choose(searched(start_session("CCOC"), "root", "CCO.CI", "CCBr.CO"))
    root_search = plan.searches[0].search_id
    branch = plan.root.children[0]
    sibling = plan.root.children[1]
    plan = edit_session(plan, "stop", sibling.node_id)
    plan = choose(searched(plan, branch.node_id, "CCBr.O"), branch.node_id)
    child = plan.root.children[0].children[0]
    plan = choose(searched(plan, child.node_id, "C.CBr"), child.node_id)
    assert export_route(plan).reaction_count == 3
    assert plan.root.children[1].stopped
    previous = plan
    replaced = edit_session(plan, "select", "root", PlanningChoice(root_search, 1, 0))
    assert export_route(replaced).reaction_count == 1
    undone = edit_session(replaced, "undo")
    assert undone.root == previous.root
    assert edit_session(undone, "redo").root == replaced.root
    restored = restore_session(undone.to_dict())
    assert restored == undone
    tree = export_route(restored)
    assert_valid_route_tree(ReactionRouteTree.from_dict(tree.to_dict()))
    assert tree.root.reaction.evidence.evidence_kind == "predicted"
    assert session_response(restored)["summary"]["unresolved_count"] == 3


def test_repeated_molecules_are_independent_occurrences_with_shared_search():
    plan = choose(searched(start_session("CC"), "root", "C.C"))
    first, second = plan.root.children
    assert first.smiles == second.smiles
    assert first.node_id != second.node_id
    plan = edit_session(plan, "stop", first.node_id)
    assert plan.root.children[0].stopped and not plan.root.children[1].stopped
    plan = searched(plan, second.node_id, "[CH3]I")
    selected = choose(plan, second.node_id)
    assert selected.root.children[0].stopped
    assert selected.root.children[1].choice


def test_cycle_rejected_by_graph_identity_and_siblings_preserved():
    plan = choose(searched(start_session("CCO"), "root", "CCBr.O"))
    branch = plan.root.children[0]
    plan = searched(plan, branch.node_id, "OCC")
    with pytest.raises(ValueError, match="cycle"):
        choose(plan, branch.node_id)
    assert plan.root.children[0].choice is None


@pytest.mark.parametrize("target", ["", "bad smiles", "CC.O", "CC>>CO"])
def test_target_must_be_one_valid_connected_structure(target):
    with pytest.raises(ValueError):
        start_session(target)


@pytest.mark.parametrize(
    "field,value",
    [
        ("target_smiles", "CCN"),
        ("proposed_reaction_smiles", "C.O>>CCO"),
        ("forward_validation_status", "unresolved"),
        ("precursor_smiles", "not-a-molecule"),
    ],
)
def test_conflicting_candidate_evidence_is_rejected(field, value):
    result = search_result("CCO", "CCBr.O")
    result["strategies"][0]["representative"][field] = value
    with pytest.raises(ValueError):
        record_search(start_session("CCO"), "root", {}, result)


def test_import_rejects_missing_precursors_duplicates_and_invalid_history():
    plan = choose(searched(start_session("CCO"), "root", "CCBr.O"))
    raw = plan.to_dict()
    raw["root"]["children"].pop()
    with pytest.raises(ValueError, match="every precursor"):
        restore_session(raw)
    raw = plan.to_dict()
    raw["root"]["children"][1]["node_id"] = raw["root"]["children"][0]["node_id"]
    with pytest.raises(ValueError, match="Duplicate molecule"):
        restore_session(raw)
    raw = deepcopy(plan.to_dict())
    raw["past"][0]["smiles"] = "N"
    with pytest.raises(ValueError, match="preserve the target"):
        restore_session(raw)


def test_empty_results_are_unresolved_and_user_stop_is_not_stock():
    plan = searched(start_session("CCO"), "root")
    assert session_response(plan)["summary"]["unresolved_count"] == 1
    stopped = edit_session(plan, "stop")
    assert (
        export_route(stopped).root.terminal_evidence
        == "user_designated_starting_material"
    )
    assert session_response(stopped)["summary"]["unresolved_count"] == 0
    assert not edit_session(stopped, "reopen").root.stopped


def test_route_identity_excludes_volatile_search_metadata():
    result = search_result("CCO", "CCBr.O")
    first = choose(record_search(start_session("CCO"), "root", {}, result))
    result = {**result, "elapsed_seconds": 10}
    second = choose(record_search(start_session("CCO"), "root", {}, result))
    assert first.root.children[0].node_id != second.root.children[0].node_id
    assert export_route(first).tree_id == export_route(second).tree_id


def test_malformed_saved_evidence_is_rejected_before_replacing_a_session():
    raw = start_session("CCO").to_dict()
    raw["conditions"] = {
        "saved-choice": {"status": "recommended_direct", "recommendations": "corrupt"}
    }
    with pytest.raises(ValueError, match="recommendations must be a list"):
        restore_session(raw)
    plan = searched(start_session("CCO"), "root", "CCBr.O")
    raw = plan.to_dict()
    raw["searches"][0]["result"]["warnings"] = "corrupt"
    with pytest.raises(ValueError, match="warnings must be a list"):
        restore_session(raw)


def test_settings_and_history_are_deterministic_and_strict():
    assert normalize_settings({})["library_mode"] == "full"
    with pytest.raises(ValueError):
        normalize_settings({"top_k": True})
    with pytest.raises(ValueError):
        normalize_settings({"arbitrary_rule": True})
    session = start_session("[CH3:1][OH:2]")
    assert session.root.smiles == "CO"
    for _ in range(20):
        session = edit_session(edit_session(session, "stop"), "reopen")
    assert len(session.past) == 30
    undone = edit_session(session, "undo")
    changed = edit_session(undone, "reopen")
    assert not changed.future
    assert (
        export_route(restore_session(session.to_dict())).tree_id
        == export_route(session).tree_id
    )


@pytest.mark.parametrize(
    "reaction",
    [
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        "Brc1ccccc1.CN>>CNc1ccccc1",
        "Brc1ccccc1.CO>>COc1ccccc1",
        "Brc1ccccc1.CS>>CSc1ccccc1",
    ],
)
def test_actual_grouped_operator_results_can_be_selected_and_restored(reaction):
    from core_retrosynthesis.generic_library import build_generic_library
    from core_retrosynthesis.strategy_search import disconnect_strategies_detailed

    library = build_generic_library(
        [
            {
                "reaction_id": "source-1",
                "reference_id": "reference-1",
                "reaction_smiles": reaction,
            }
        ],
        levels=("L2", "L1", "L0"),
        admission_mode="data_driven",
    )
    target = reaction.split(">>")[1]
    result = disconnect_strategies_detailed(
        target, library, top_k_strategies=2
    ).to_dict()
    assert result["strategies"]
    result["target_smiles"] = target
    session = record_search(start_session(target), "root", {}, result)
    selected = choose(session)
    assert export_route(selected).reaction_count == 1
    assert restore_session(selected.to_dict()) == selected
