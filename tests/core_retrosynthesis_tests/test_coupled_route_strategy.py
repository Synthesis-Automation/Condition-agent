"""Coupled two-step route-strategy mining regressions."""

from __future__ import annotations

import json
from dataclasses import replace

import pytest

from core_retrosynthesis.cli import _parser
from core_retrosynthesis.coupled_route_strategy import (
    extract_coupled_route_strategies,
    mine_coupled_route_strategy_poc,
    write_coupled_route_strategy_report,
)
from core_retrosynthesis.coupled_route_strategy_review import (
    load_coupled_route_strategy_report,
    render_coupled_route_strategy_html,
    write_coupled_route_strategy_html,
)
from core_retrosynthesis.route_core import (
    RouteAtomLineageCandidate,
    RouteCoreLineageLink,
    RouteCoreProjection,
    RouteCoreStep,
)


def _atom(map_number: int) -> dict[str, int]:
    return {"atom_map_number": map_number}


def _step(
    node: str,
    reaction: str,
    *,
    active_maps: tuple[int, ...],
    formed_maps: tuple[int, ...] = (),
    key: str,
) -> RouteCoreStep:
    edits = []
    if formed_maps:
        edits.append(
            {
                "edit_type": "formed",
                "atom_1": _atom(formed_maps[0]),
                "atom_2": _atom(formed_maps[1]),
            }
        )
    for map_number in active_maps:
        if map_number not in formed_maps:
            edits.append(
                {
                    "edit_type": "hydrogen_change",
                    "atom_1": _atom(map_number),
                    "atom_2": None,
                }
            )
    return RouteCoreStep(
        reaction_node_id=node,
        step_id=f"step-{node}",
        source_reaction_id=f"source-{node}",
        retrosynthetic_depth=1,
        observed_remaining_steps=1,
        product_occurrence_id=f"product-{node}",
        precursor_occurrence_ids=(),
        internal_precursor_occurrence_ids=(),
        terminal_precursor_occurrence_ids=(),
        reaction_smiles=reaction,
        reagents_smiles="",
        chemistry_valid=True,
        evidence_quality="validated_atom_mapping",
        quality_status="pass",
        reaction_signature_id=f"RS3:{key}",
        reaction_core_id=f"RCP2:{key}",
        exact_core_key=f"RCX2:{key}",
        typed_core_key=f"RCT2:{key}",
        shape_core_key=f"RSH2:{key}",
        center_transition_key=f"RCS2:{key}",
        transformation_class="generic_graph_transformation",
        named_family=None,
        event_count=1,
        active_atom_count=len(active_maps),
        edit_tokens=tuple(f"active:{value}" for value in active_maps),
        reaction_signature={"edits": edits},
        reaction_core=None,
        minimum_reaction_smiles=reaction,
        render_reaction_smiles=reaction,
        display_definition_id="test",
        display_retained_atom_count=3,
        display_removed_substituent_count=0,
        display_r_group_count=0,
        abstracted_reaction_smiles=reaction,
        abstraction_status="available",
        warnings=(),
    )


def _projection(
    *,
    producer_active: tuple[int, ...],
    consumer_active: tuple[int, ...],
    candidates: tuple[tuple[tuple[int, int], ...], ...],
    status: str = "unique",
    formed_maps: tuple[int, ...] = (),
    route_id: str = "route-1",
    patent_id: str = "US-1",
) -> RouteCoreProjection:
    first = _step(
        "producer",
        "[CH3:1][CH2:2][CH2:3][Br:5].[NH2:4]>>[CH3:1][CH2:2][CH2:3][NH2:4]",
        active_maps=producer_active,
        formed_maps=formed_maps,
        key="first",
    )
    second = _step(
        "consumer",
        "[CH3:10][CH2:11][CH2:12][NH2:13]>>[CH3:10][CH2:11][CH:12]=[NH:13]",
        active_maps=consumer_active,
        key="second",
    )
    lineage = RouteCoreLineageLink(
        lineage_id="lineage-1",
        intermediate_occurrence_id="intermediate-1",
        intermediate_smiles="CCCN",
        producer_reaction_node_id="producer",
        consumer_reaction_node_id="consumer",
        producer_product_component_index=0,
        consumer_reactant_component_index=0,
        status=status,
        candidate_count=len(candidates),
        candidate_limit_reached=False,
        candidates=tuple(
            RouteAtomLineageCandidate(
                candidate_id=f"candidate-{index}", atom_map_pairs=pairs
            )
            for index, pairs in enumerate(candidates, 1)
        ),
    )
    return RouteCoreProjection(
        route_core_id=f"route-core-{route_id}",
        source_tree_id=f"tree-{route_id}",
        source_route_id=route_id,
        patent_id=patent_id,
        split="test",
        target_smiles="CC=NC",
        route_kind="observed",
        reaction_count=2,
        maximum_depth=2,
        exact_route_core_key="exact-route",
        typed_route_core_key="typed-route",
        shape_route_core_key="shape-route",
        steps=(first, second),
        lineage_links=(lineage,),
        motifs=(),
        typed_motif_keys=(),
        shape_motif_keys=(),
        fully_chemistry_resolved=True,
        fully_lineage_connected=True,
        warnings=(),
    )


def _positive_projection(**kwargs) -> RouteCoreProjection:
    return _projection(
        producer_active=(3, 4),
        consumer_active=(12, 13),
        formed_maps=(3, 4),
        candidates=(((1, 10), (2, 11), (3, 12), (4, 13)),),
        **kwargs,
    )


def test_formed_handle_consumed_at_same_site_is_strict_strategy() -> None:
    occurrence = extract_coupled_route_strategies(_positive_projection())[0]

    assert occurrence.relationship_class == "handle_progression"
    assert occurrence.admission_class == "strict"
    assert occurrence.coupling_score == 1.0
    assert occurrence.evidence.overlap_counts == (2,)
    assert occurrence.evidence.formed_overlap_counts == (2,)
    assert occurrence.overall_reaction_smiles
    assert occurrence.typed_strategy_id.startswith("CRST1:")


def test_distant_active_sites_are_rejected_as_independent() -> None:
    projection = _projection(
        producer_active=(1,),
        consumer_active=(13,),
        candidates=(((1, 10), (2, 11), (3, 12), (4, 13)),),
    )

    occurrence = extract_coupled_route_strategies(projection)[0]

    assert occurrence.relationship_class == "independent_sites"
    assert occurrence.admission_class == "rejected"
    assert occurrence.evidence.minimum_distances == (3,)


def test_lineage_dependent_site_overlap_is_retained_for_review() -> None:
    projection = _projection(
        producer_active=(1,),
        consumer_active=(13,),
        candidates=(
            ((1, 13), (2, 12), (3, 11), (4, 10)),
            ((1, 10), (2, 11), (3, 12), (4, 13)),
        ),
        status="symmetry_ambiguous",
    )

    occurrence = extract_coupled_route_strategies(projection)[0]

    assert occurrence.relationship_class == "lineage_ambiguous"
    assert occurrence.admission_class == "review"
    assert occurrence.evidence.overlap_counts == (1, 0)
    assert not occurrence.evidence.ambiguity_invariant
    assert "LINEAGE_RELATIONSHIP_NOT_INVARIANT" in occurrence.warnings


def test_strategy_identity_is_deterministic() -> None:
    projection = _positive_projection()

    first = extract_coupled_route_strategies(projection)[0]
    second = extract_coupled_route_strategies(projection)[0]

    assert first == second
    assert first.occurrence_id.startswith("CRSO1:")


def test_mining_aggregates_patent_disjoint_support_and_balanced_sample(
    tmp_path, monkeypatch
) -> None:
    first = _positive_projection(route_id="positive-a", patent_id="US-A")
    second = replace(
        first,
        route_core_id="route-core-positive-b",
        source_tree_id="tree-positive-b",
        source_route_id="positive-b",
        patent_id="US-B",
    )
    independent = _projection(
        producer_active=(1,),
        consumer_active=(13,),
        candidates=(((1, 10), (2, 11), (3, 12), (4, 13)),),
        route_id="negative",
    )
    ambiguous = _projection(
        producer_active=(1,),
        consumer_active=(13,),
        candidates=(
            ((1, 13), (2, 12), (3, 11), (4, 10)),
            ((1, 10), (2, 11), (3, 12), (4, 13)),
        ),
        status="symmetry_ambiguous",
        route_id="ambiguous",
    )
    monkeypatch.setattr(
        "core_retrosynthesis.coupled_route_strategy.iter_route_core_projections",
        lambda _: iter((first, second, independent, ambiguous)),
    )
    source = tmp_path / "source.jsonl.gz"
    source.write_bytes(b"test")

    report = mine_coupled_route_strategy_poc(
        source,
        strict_sample_size=1,
        review_sample_size=1,
        rejected_sample_size=1,
        required_route_ids=(),
    )

    assert report["route_count"] == 4
    assert report["lineage_pair_count"] == 4
    assert report["admission_counts"] == {
        "rejected": 1,
        "review": 1,
        "strict": 2,
    }
    assert report["recurrent_strict_strategy_count"] == 1
    assert report["sample_counts"] == {
        "rejected": 1,
        "review": 1,
        "strict": 1,
    }


def test_report_and_review_artifacts_are_deterministic(tmp_path) -> None:
    occurrence = extract_coupled_route_strategies(_positive_projection())[0]
    report = {
        "artifact_type": "coupled_route_strategy_poc",
        "algorithm_version": "test",
        "source_sha256": "abc123",
        "sample_seed": 7,
        "route_count": 1,
        "lineage_pair_count": 1,
        "admission_counts": {"strict": 1},
        "sample": [occurrence.to_dict()],
    }
    json_path = tmp_path / "report.json"
    html_path = tmp_path / "review.html"

    json_summary = write_coupled_route_strategy_report(report, json_path)
    html_summary = write_coupled_route_strategy_html(report, html_path)
    document = html_path.read_text(encoding="utf-8")

    assert json_summary["sample_count"] == 1
    assert html_summary["sample_count"] == 1
    assert "Logical overall transformation" in document
    assert "Physical step 1" in document
    assert "Carried intermediate" in document
    assert "Structural coupling evidence" in document
    assert "Export review JSON" in document
    assert "localStorage" in document
    assert 'class="render-error"' not in document
    assert render_coupled_route_strategy_html(report) == document


def test_loader_and_cli_contract(tmp_path) -> None:
    source = tmp_path / "wrong.json"
    source.write_text(json.dumps({"artifact_type": "other"}), encoding="utf-8")
    with pytest.raises(ValueError, match="unexpected.*artifact type"):
        load_coupled_route_strategy_report(source)

    arguments = _parser().parse_args(
        [
            "mine-coupled-route-strategies",
            "source.jsonl.gz",
            "report.json",
            "review.html",
            "--strict-sample-size",
            "12",
        ]
    )
    assert arguments.command == "mine-coupled-route-strategies"
    assert arguments.strict_sample_size == 12
