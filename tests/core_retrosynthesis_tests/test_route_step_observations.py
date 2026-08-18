"""Canonical multistep-route physical observation extraction regressions."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from condition_recommender.conversion.intermediate import (
    iter_intermediate_records,
)
from core_retrosynthesis.cli import _parser
from core_retrosynthesis.route_conversion import build_observed_route_tree
from core_retrosynthesis.route_core import build_route_core_projection
from core_retrosynthesis.route_step_observations import (
    extract_route_step_observations,
)

from .test_route_core import _two_step_record


def _sources(
    directory: Path,
    *,
    splits: tuple[str, ...] = ("train", "validation", "test"),
) -> tuple[Path, Path]:
    trees = []
    projections = []
    for ordinal, split in enumerate(splits, 1):
        record = _two_step_record()
        record["route_id"] = f"route-{ordinal}"
        record["patent_id"] = f"PATENT-{ordinal}"
        record["split"] = split
        for step in record["steps"]:
            step["source_reaction_id"] = (
                f"{step['source_reaction_id']}-{ordinal}"
            )
        if ordinal == 1:
            reaction = record["steps"][0]["reaction_smiles"]
            reactants, product = reaction.split(">>")
            record["steps"][0]["reaction_smiles"] = (
                f"{reactants}>[Pd].O>{product}"
            )
            record["steps"][0]["reagents_smiles"] = "[Pd].O"
        tree = build_observed_route_tree(record)
        trees.append(tree)
        projections.append(build_route_core_projection(tree))
    tree_source = directory / "routes.tree.v2.jsonl.gz"
    core_source = directory / "routes.core.v1.jsonl.gz"
    with gzip.open(tree_source, "wt", encoding="utf-8") as stream:
        for tree in trees:
            stream.write(json.dumps(tree.to_dict()) + "\n")
    with gzip.open(core_source, "wt", encoding="utf-8") as stream:
        for projection in projections:
            stream.write(json.dumps(projection.to_dict()) + "\n")
    return tree_source, core_source


def _rows(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        return [json.loads(line) for line in stream]


def test_route_steps_are_patent_split_source_observations(
    tmp_path: Path,
) -> None:
    trees, cores = _sources(tmp_path)
    output = tmp_path / "observations"

    report = extract_route_step_observations(trees, cores, output)

    assert report["multistep_route_count"] == 3
    assert report["step_occurrence_count"] == 6
    assert report["unique_step_count"] == 6
    assert {
        split: report["artifacts"][split]["row_count"]
        for split in ("train", "validation", "test")
    } == {"train": 2, "validation": 2, "test": 2}
    train_path = output / "route_steps.train.observations.jsonl.gz"
    train = _rows(train_path)
    assert all(item["schema_version"] == "source_observation.v1" for item in train)
    assert all(item["source"]["source_groups"]["split"] == "train" for item in train)
    assert {item["source"]["reference"] for item in train} == {"PATENT-1"}
    with_agents = next(
        item for item in train if item["raw_fields"]["reagents_smiles"]
    )
    assert with_agents["reaction"]["reaction_smiles"].count(">>") == 1
    assert with_agents["reaction"]["reaction_smiles"].count(">") == 2
    assert {
        component["identifiers"][0]["value"]
        for component in with_agents["conditions"]["components"]
    } == {"[Pd]", "O"}
    assert all(
        component["source_role_hint"] is None
        for component in with_agents["conditions"]["components"]
    )
    converted = tuple(iter_intermediate_records(train_path))
    assert len(converted) == 2
    assert all(record.reaction_smiles.count(">>") == 1 for record in converted)
    assert all(record.reference == "PATENT-1" for record in converted)


def test_route_step_extraction_is_byte_deterministic(tmp_path: Path) -> None:
    trees, cores = _sources(tmp_path, splits=("train",))
    first = tmp_path / "first"
    second = tmp_path / "second"

    first_report = extract_route_step_observations(trees, cores, first)
    second_report = extract_route_step_observations(trees, cores, second)

    for split in ("train", "validation", "test"):
        first_path = first / f"route_steps.{split}.observations.jsonl.gz"
        second_path = second / f"route_steps.{split}.observations.jsonl.gz"
        assert first_path.read_bytes() == second_path.read_bytes()
        assert (
            first_report["artifacts"][split]["sha256"]
            == second_report["artifacts"][split]["sha256"]
        )


def test_source_reaction_cannot_cross_patent_splits(tmp_path: Path) -> None:
    trees, cores = _sources(tmp_path, splits=("train", "test"))
    tree_rows = _rows(trees)
    core_rows = _rows(cores)
    tree_rows[1]["root"]["reaction"]["evidence"]["source_reaction_id"] = (
        tree_rows[0]["root"]["reaction"]["evidence"]["source_reaction_id"]
    )
    core_rows[1]["steps"][0]["source_reaction_id"] = core_rows[0]["steps"][0][
        "source_reaction_id"
    ]
    with gzip.open(trees, "wt", encoding="utf-8") as stream:
        for row in tree_rows:
            stream.write(json.dumps(row) + "\n")
    with gzip.open(cores, "wt", encoding="utf-8") as stream:
        for row in core_rows:
            stream.write(json.dumps(row) + "\n")

    with pytest.raises(ValueError, match="across patent data splits"):
        extract_route_step_observations(trees, cores, tmp_path / "output")


def test_cli_exposes_route_step_extraction() -> None:
    arguments = _parser().parse_args(
        [
            "extract-route-step-observations",
            "trees.jsonl.gz",
            "cores.jsonl.gz",
            "output",
        ]
    )
    assert arguments.minimum_reaction_count == 2
