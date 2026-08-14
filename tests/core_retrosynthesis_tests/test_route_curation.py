"""Deterministic route-corpus curation regressions."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from core_retrosynthesis import (
    RouteQualityError,
    RouteSubsetPolicy,
    curate_route_subset,
    normalize_linear_route_record,
)


def _route(
    patent_id: str,
    route_index: int,
    reaction_offset: int,
    *,
    reduction: int = 1,
) -> dict:
    route_id = f"{patent_id}_{route_index}"
    subtree_id = f"{route_id}_0"
    source_ids = [
        f"{reaction_offset + index}_0" for index in range(3)
    ]
    root_step = {
        "_id": f"{subtree_id}_{source_ids[0]}",
        "reaction_smiles": (
            "[CH3:1][CH2:2][CH3:3].[OH2:4]>O>"
            "[CH3:1][CH2:2][CH2:3][OH:4]"
        ),
        "abstracted_reaction_smiles": (
            "[CH3:1][CH2:2][CH3:3].[OH2:4]>>"
            "[CH3:1][CH2:2][CH2:3][OH:4]"
        ),
    }
    middle_step = {
        "_id": f"{subtree_id}_{source_ids[1]}",
        "reaction_smiles": (
            "[CH3:1][CH3:2].[CH4:3]>N>"
            "[CH3:1][CH2:2][CH3:3]"
        ),
        "abstracted_reaction_smiles": "" if reduction else (
            "[CH3:1][CH3:2].[CH4:3]>>[CH3:1][CH2:2][CH3:3]"
        ),
    }
    leaf_step = {
        "_id": f"{subtree_id}_{source_ids[2]}",
        "reaction_smiles": (
            "[CH4:1].[CH4:2]>[Na+]>[CH3:1][CH3:2]"
        ),
        "abstracted_reaction_smiles": (
            "[CH4:1].[CH4:2]>>[CH3:1][CH3:2]"
        ),
    }
    higher_level_count = 3 - reduction
    return {
        "route_id": route_id,
        "original_tree": {
            "depth": 3,
            "num_reactions": 3,
            "reaction_ids": source_ids,
        },
        "subtrees": [
            {
                "subtree_id": subtree_id,
                "num_reactions": [3, higher_level_count],
                "depth": [3, higher_level_count],
                # Source array order is deliberately not synthesis order.
                "reactions": [leaf_step, root_step, middle_step],
            }
        ],
    }


def _write_routes(path: Path, routes: list[dict]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for route in routes:
            handle.write(json.dumps(route))
            handle.write("\n")


def _read_routes(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle]


def test_normalization_reconstructs_order_and_retains_reagents() -> None:
    normalized = normalize_linear_route_record(
        _route("US00000001A1", 0, 100),
        split="train",
    )

    assert normalized["target_smiles"] == "CCCO"
    assert [step["product_smiles"] for step in normalized["steps"]] == [
        "CCCO",
        "CCC",
        "CC",
    ]
    assert [step["reagents_smiles"] for step in normalized["steps"]] == [
        "O",
        "N",
        "[Na+]",
    ]
    assert normalized["quality"]["source_array_order_used"] is False
    assert normalized["abstraction_reduction"] == 1


def test_normalization_rejects_an_ambiguous_target() -> None:
    route = _route("US00000001A1", 0, 100)
    subtree = route["subtrees"][0]
    subtree["reactions"][2]["reaction_smiles"] = (
        "[CH4:1].[OH2:2]>N>[CH3:1][OH:2]"
    )

    with pytest.raises(RouteQualityError) as error:
        normalize_linear_route_record(route, split="train")

    assert error.value.reason == "ambiguous_route_target"


def test_curator_is_deterministic_and_patent_disjoint(tmp_path: Path) -> None:
    source = tmp_path / "source.jsonl.gz"
    first_output = tmp_path / "first.jsonl.gz"
    second_output = tmp_path / "second.jsonl.gz"
    routes = [
        _route("US00000001A1", 0, 100),
        _route("US00000001A1", 1, 200),
        _route("US00000002A1", 0, 300),
        _route("US00000003A1", 0, 400, reduction=0),
    ]
    _write_routes(source, routes)
    policy = RouteSubsetPolicy(maximum_routes=2)

    first = curate_route_subset(source, first_output, policy=policy)
    second = curate_route_subset(source, second_output, policy=policy)

    first_records = _read_routes(first_output)
    assert len(first_records) == 2
    assert len({record["patent_id"] for record in first_records}) == 2
    assert all(record["abstraction_reduction"] == 1 for record in first_records)
    assert first["counts"]["selected_patents"] == 2
    assert first["counts"]["scan"]["rejected_no_abstraction_reduction"] == 1
    assert (
        first["output"]["route_file_sha256"]
        == second["output"]["route_file_sha256"]
    )


def test_curator_excludes_test_targets(tmp_path: Path) -> None:
    source = tmp_path / "source.jsonl.gz"
    output = tmp_path / "output.jsonl.gz"
    testset = tmp_path / "testset.smi"
    _write_routes(
        source,
        [
            _route("US00000001A1", 0, 100),
            _route("US00000002A1", 0, 200),
        ],
    )
    testset.write_text("CCC.O>>CCCO\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="Only 0 routes passed"):
        curate_route_subset(
            source,
            output,
            testset_path=testset,
            policy=RouteSubsetPolicy(maximum_routes=1),
        )
