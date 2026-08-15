"""Route-action replay corpus conversion regressions."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from reactive_taxonomy import featurize_reaction

from core_retrosynthesis.generic_library import (
    build_generic_library,
    save_generic_library,
)
from core_retrosynthesis.route_action_conversion import (
    convert_route_action_corpus,
    iter_route_action_evaluations,
    merge_route_action_shards,
    stable_route_shard,
)
from core_retrosynthesis.route_action_evaluation import RouteActionEvaluationConfig
from core_retrosynthesis.route_conversion import build_observed_route_tree

from .test_route_action_evaluation import REACTION, _record


def _write_inputs(directory: Path) -> tuple[Path, Path]:
    tree_path = directory / "trees.jsonl.gz"
    library_path = directory / "library.json.gz"
    tree = build_observed_route_tree(_record())
    with gzip.open(tree_path, "wt", encoding="utf-8") as handle:
        handle.write(json.dumps(tree.to_dict()) + "\n")
    analysis = featurize_reaction(REACTION).to_dict()
    library = build_generic_library(
        (
            {
                "reaction_id": "training-reaction",
                "reference_id": "US99999999A1",
                "reaction_smiles": REACTION,
                "reaction_core": analysis["reaction_core"],
                "reaction_observation": analysis["observation"],
                "reaction_completeness": analysis["reaction_completeness"],
            },
        ),
        levels=("L0", "L1", "L2"),
        admission_mode="data_driven",
    )
    save_generic_library(library, library_path)
    return tree_path, library_path


def _write_many_inputs(directory: Path, count: int = 6) -> tuple[Path, Path]:
    tree_path, library_path = _write_inputs(directory)
    trees = []
    for index in range(count):
        record = _record()
        record["route_id"] = f"US00000001A1_{index:02d}"
        trees.append(build_observed_route_tree(record))
    with gzip.open(tree_path, "wt", encoding="utf-8") as handle:
        for tree in trees:
            handle.write(json.dumps(tree.to_dict()) + "\n")
    return tree_path, library_path


def test_route_action_corpus_is_deterministic_and_reloadable(
    tmp_path: Path,
) -> None:
    trees, library = _write_inputs(tmp_path)
    first = tmp_path / "first.jsonl.gz"
    second = tmp_path / "second.jsonl.gz"
    config = RouteActionEvaluationConfig(
        top_k=5,
        max_templates_to_apply=10,
        max_candidates_to_validate=10,
    )

    first_report = convert_route_action_corpus(
        trees, library, first, config=config
    )
    second_report = convert_route_action_corpus(
        trees, library, second, config=config
    )

    assert first.read_bytes() == second.read_bytes()
    assert first_report["output"]["sha256"] == second_report["output"]["sha256"]
    assert first_report["metrics"]["route_count"] == 1
    assert first_report["metrics"]["step_count"] == 1
    assert first_report["metrics"]["eligibility_facet_counts"]["strategy"] == 1
    assert first_report["metrics"]["eligibility_facet_counts"][
        "operator_roundtrip"
    ] == 1
    assert first_report["metrics"]["recall"]["exact_precursor"]["top1"] == 1.0
    restored = tuple(iter_route_action_evaluations(first))
    assert len(restored) == 1
    assert restored[0].steps[0].exact_precursor_rank == 1


def test_labels_only_conversion_reports_facets_without_false_zero_recall(
    tmp_path: Path,
) -> None:
    trees, library = _write_inputs(tmp_path)
    output = tmp_path / "labels.jsonl.gz"

    report = convert_route_action_corpus(
        trees,
        library,
        output,
        config=RouteActionEvaluationConfig(run_search=False),
    )

    metrics = report["metrics"]
    assert metrics["search_eligible_step_count"] == 1
    assert metrics["searched_step_count"] == 0
    assert metrics["candidate_coverage"] is None
    assert metrics["recall"]["strategy"]["denominator"] == 0
    assert metrics["recall"]["strategy"]["top1"] is None
    restored = tuple(iter_route_action_evaluations(output))
    assert restored[0].steps[0].search_status == "not_run"
    assert restored[0].steps[0].observed_action.strategy_verified


def test_stable_route_shards_are_deterministic_and_bounded() -> None:
    first = [stable_route_shard(f"route-{index}", 7) for index in range(100)]
    second = [stable_route_shard(f"route-{index}", 7) for index in range(100)]

    assert first == second
    assert set(first) == set(range(7))
    with pytest.raises(ValueError, match="positive"):
        stable_route_shard("route", 0)


def test_shards_resume_and_merge_to_direct_deterministic_output(
    tmp_path: Path,
) -> None:
    trees, library = _write_many_inputs(tmp_path)
    config = RouteActionEvaluationConfig(run_search=False)
    shards = []
    for shard_index in range(3):
        output = tmp_path / f"part-{shard_index}.jsonl.gz"
        report = convert_route_action_corpus(
            trees,
            library,
            output,
            config=config,
            shard_count=3,
            shard_index=shard_index,
        )
        resumed = convert_route_action_corpus(
            trees,
            library,
            output,
            config=config,
            shard_count=3,
            shard_index=shard_index,
            resume=True,
        )
        assert resumed["output"]["sha256"] == report["output"]["sha256"]
        shards.append(output)

    merged = tmp_path / "merged.jsonl.gz"
    merged_report = merge_route_action_shards(shards, merged)
    direct = tmp_path / "direct.jsonl.gz"
    direct_report = convert_route_action_corpus(
        trees,
        library,
        direct,
        config=config,
    )

    assert merged.read_bytes() == direct.read_bytes()
    assert merged_report["output"]["sha256"] == direct_report["output"]["sha256"]
    assert merged_report["metrics"]["route_count"] == 6
    assert merged_report["selection"]["merged_shard_indices"] == [0, 1, 2]


def test_merge_rejects_incomplete_and_duplicate_shard_sets(tmp_path: Path) -> None:
    trees, library = _write_many_inputs(tmp_path, count=2)
    shards = []
    for shard_index in range(2):
        output = tmp_path / f"part-{shard_index}.jsonl.gz"
        convert_route_action_corpus(
            trees,
            library,
            output,
            config=RouteActionEvaluationConfig(run_search=False),
            shard_count=2,
            shard_index=shard_index,
        )
        shards.append(output)

    with pytest.raises(ValueError, match="incomplete"):
        merge_route_action_shards(shards[:1], tmp_path / "missing.jsonl.gz")
    with pytest.raises(ValueError, match="duplicate shard index"):
        merge_route_action_shards(
            [shards[0], shards[0]], tmp_path / "duplicate.jsonl.gz"
        )
