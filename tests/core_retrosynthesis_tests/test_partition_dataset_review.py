"""Observed-route K-way dataset review regressions."""

from __future__ import annotations

import gzip
import json

from core_retrosynthesis import (
    MoleculeOccurrenceNode,
    ReactionRouteTree,
    RouteReactionNode,
    RouteStepEvidence,
    build_partition_dataset_review,
    render_partition_dataset_review_html,
    write_partition_dataset_review,
)


REACTION = (
    "[CH3:1][C:2](=[O:3])[OH:6]."
    "[NH2:4][CH2:5][CH3:7]>O>"
    "[CH3:1][C:2](=[O:3])[NH:4][CH2:5][CH3:7]"
)


def _tree(suffix: str) -> ReactionRouteTree:
    acid = MoleculeOccurrenceNode(
        occurrence_id=f"acid:{suffix}",
        smiles="CC(=O)O",
        depth=1,
        terminal=True,
        terminal_evidence="observed_leaf",
        unresolved_reason=None,
    )
    amine = MoleculeOccurrenceNode(
        occurrence_id=f"amine:{suffix}",
        smiles="CCN",
        depth=1,
        terminal=True,
        terminal_evidence="observed_leaf",
        unresolved_reason=None,
    )
    reaction = RouteReactionNode(
        reaction_node_id=f"reaction:{suffix}",
        step_id=f"step:{suffix}",
        depth=1,
        reaction_smiles=REACTION,
        evidence=RouteStepEvidence(
            evidence_kind="observed",
            source_reaction_id=f"source-reaction:{suffix}",
        ),
        children=(acid, amine),
    )
    root = MoleculeOccurrenceNode(
        occurrence_id=f"target:{suffix}",
        smiles="CCNC(C)=O",
        depth=0,
        terminal=False,
        terminal_evidence="expanded",
        unresolved_reason=None,
        reaction=reaction,
    )
    return ReactionRouteTree(
        tree_id=f"tree:{suffix}",
        route_kind="observed",
        target_smiles="CCNC(C)=O",
        root=root,
        reaction_count=1,
        maximum_depth=1,
        fingerprint_tokens=(),
        source_dataset_id="test-routes",
        source_route_id=f"route:{suffix}",
        patent_id=f"patent:{suffix}",
        split="test",
    )


def _write_trees(path, trees) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="\n") as handle:
        for tree in trees:
            handle.write(json.dumps(tree.to_dict(), sort_keys=True) + "\n")


def test_partition_dataset_review_selects_and_renders_deterministically(
    tmp_path,
) -> None:
    source = tmp_path / "routes.jsonl.gz"
    _write_trees(source, (_tree("a"), _tree("b")))

    review = build_partition_dataset_review(
        source,
        sample_size=2,
        seed=17,
        minimum_k=2,
    )
    repeated = build_partition_dataset_review(
        source,
        sample_size=2,
        seed=17,
        minimum_k=2,
    )

    assert review.to_dict() == repeated.to_dict()
    assert review.source_route_count == 2
    assert review.eligible_route_count == 2
    assert review.projection_count_by_maximum_k == ((2, 2),)
    assert all(case.fully_projected for case in review.cases)
    assert all(
        [frontier.partition.k for frontier in case.projection.frontiers]
        == [1, 2]
        for case in review.cases
    )
    html = render_partition_dataset_review_html(review)
    assert "2-route graphical fragmentation review" in html
    assert "maximum successfully projected k=2" in html
    assert html.count("<article class='case'") == 2
    assert html.count("<section class='frontier-panel'") == 4
    assert html.count("<svg") >= 12
    assert "Export review JSON" in html


def test_partition_dataset_review_writes_json_and_html(tmp_path) -> None:
    source = tmp_path / "routes.jsonl.gz"
    output_json = tmp_path / "review.json"
    output_html = tmp_path / "review.html"
    _write_trees(source, (_tree("one"),))
    review = build_partition_dataset_review(
        source,
        sample_size=1,
        minimum_k=2,
    )

    write_partition_dataset_review(review, output_json, output_html)

    assert json.loads(output_json.read_text("utf-8")) == review.to_dict()
    assert "partition-dataset-review-data" in output_html.read_text("utf-8")
