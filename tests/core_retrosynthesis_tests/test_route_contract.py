"""Shared evidence-neutral route-tree contract regressions."""

from __future__ import annotations

from core_retrosynthesis import (
    MoleculeOccurrenceNode,
    ReactionRouteTree,
    RouteReactionNode,
    RouteStepEvidence,
    validate_route_tree,
)


def _tree(*, duplicate_child: bool = False) -> ReactionRouteTree:
    first = MoleculeOccurrenceNode(
        occurrence_id="leaf-a",
        smiles="C",
        depth=1,
        terminal=True,
        terminal_evidence="source_leaf",
        unresolved_reason=None,
    )
    second = first if duplicate_child else MoleculeOccurrenceNode(
        occurrence_id="leaf-b",
        smiles="O",
        depth=1,
        terminal=True,
        terminal_evidence="source_leaf",
        unresolved_reason=None,
    )
    reaction = RouteReactionNode(
        reaction_node_id="reaction-1",
        step_id="step-1",
        depth=1,
        reaction_smiles="C.O>>CO",
        evidence=RouteStepEvidence(evidence_kind="observed"),
        children=(first, second),
    )
    root = MoleculeOccurrenceNode(
        occurrence_id="root",
        smiles="CO",
        depth=0,
        terminal=False,
        terminal_evidence="expanded",
        unresolved_reason=None,
        reaction=reaction,
    )
    return ReactionRouteTree(
        tree_id="tree-1",
        route_kind="observed",
        target_smiles="CO",
        root=root,
        reaction_count=1,
        maximum_depth=1,
        fingerprint_tokens=("reaction:1",),
    )


def test_shared_route_tree_validates_and_serializes_neutral_evidence() -> None:
    tree = _tree()

    report = validate_route_tree(tree)
    serialized = tree.to_dict()

    assert report.valid is True
    assert report.molecule_occurrence_count == 3
    assert serialized["route_kind"] == "observed"
    assert serialized["schema_version"] == "2.0"
    assert serialized["root"]["occurrence_id"] == "root"
    reaction = serialized["root"]["reaction"]
    assert reaction["evidence"]["evidence_kind"] == "observed"
    assert reaction["planned_action"] is None
    assert reaction["operator_id"] is None
    assert reaction["proposed_reaction_smiles"] == reaction["reaction_smiles"]
    restored = ReactionRouteTree.from_dict(serialized)
    assert restored == tree


def test_route_tree_validator_rejects_reused_occurrence_identity() -> None:
    report = validate_route_tree(_tree(duplicate_child=True))

    assert report.valid is False
    assert any(
        issue.startswith("duplicate_molecule_occurrence_id:")
        for issue in report.issues
    )
