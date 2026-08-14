"""Observed route-to-tree conversion regressions."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

from core_retrosynthesis import (
    build_observed_route_tree,
    convert_observed_route_corpus,
    iter_molecule_occurrences,
    iter_route_reactions,
    iter_route_trees,
    validate_route_tree,
)


def _step(source_id: str, product: str, precursors: list[str]) -> dict:
    reactants = ".".join(precursors)
    reaction = f"{reactants}>>{product}"
    return {
        "retrosynthetic_position": 0,
        "source_reaction_id": source_id,
        "reaction_smiles": reaction,
        "abstracted_reaction_smiles": reaction,
        "reactants_smiles": reactants,
        "reagents_smiles": "",
        "product_smiles_mapped": product,
        "product_smiles": product,
        "precursor_smiles": precursors,
        "internal_precursor_smiles": [],
        "terminal_precursor_smiles": [],
    }


def _branched_record(route_id: str = "US00000001A1_0") -> dict:
    return {
        "schema_version": "1.0",
        "route_id": route_id,
        "patent_id": route_id.rsplit("_", 1)[0],
        "split": "train",
        "target_smiles": "CCOCO",
        "original_reaction_count": 3,
        "original_depth": 2,
        "higher_level_reaction_count": 2,
        "higher_level_depth": 2,
        "abstraction_reduction": 1,
        # Array order is deliberately unrelated to topology.
        "steps": [
            _step("leaf-b", "CO", ["C", "O"]),
            _step("root", "CCOCO", ["CCO", "CO"]),
            _step("leaf-a", "CCO", ["CC", "O"]),
        ],
    }


def test_observed_converter_builds_a_branched_occurrence_tree() -> None:
    tree = build_observed_route_tree(_branched_record())

    assert tree.route_kind == "observed"
    assert tree.maximum_depth == 2
    assert tree.reaction_count == 3
    assert validate_route_tree(tree).valid is True
    assert tree.root.reaction is not None
    assert len(tree.root.reaction.children) == 2
    assert all(child.reaction is not None for child in tree.root.reaction.children)
    assert len(tuple(iter_route_reactions(tree))) == 3
    occurrences = tuple(iter_molecule_occurrences(tree))
    assert len(occurrences) == 7
    assert len({node.occurrence_id for node in occurrences}) == 7
    assert all(
        reaction.evidence.connectivity_method
        == "reconstructed_from_canonical_molecule_identity"
        for reaction in iter_route_reactions(tree)
    )


def test_duplicate_terminal_molecules_remain_distinct_occurrences() -> None:
    record = _branched_record()
    record["target_smiles"] = "CC"
    record["original_reaction_count"] = 1
    record["steps"] = [_step("root", "CC", ["C", "C"])]

    tree = build_observed_route_tree(record)
    children = tree.root.reaction.children

    assert [child.smiles for child in children] == ["C", "C"]
    assert children[0].occurrence_id != children[1].occurrence_id


def _write_source(path: Path, records: list[dict]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record) + "\n")


def test_corpus_conversion_is_deterministic(tmp_path: Path) -> None:
    source = tmp_path / "curated.jsonl.gz"
    first = tmp_path / "first.jsonl.gz"
    second = tmp_path / "second.jsonl.gz"
    records = [
        _branched_record("US00000001A1_0"),
        _branched_record("US00000002A1_0"),
        _branched_record("US00000003A1_0"),
    ]
    _write_source(source, list(reversed(records)))

    first_report = convert_observed_route_corpus(
        source,
        first,
        sample_size=2,
        seed=11,
    )
    second_report = convert_observed_route_corpus(
        source,
        second,
        sample_size=2,
        seed=11,
    )

    assert first_report["conversion"]["converted_count"] == 2
    assert first_report["conversion"]["rejected_count"] == 0
    assert (
        first_report["output"]["sha256"]
        == second_report["output"]["sha256"]
    )
    with gzip.open(first, "rt", encoding="utf-8") as handle:
        converted = [json.loads(line) for line in handle]
    assert all(row["schema_version"] == "2.0" for row in converted)
    assert all(row["route_kind"] == "observed" for row in converted)
    restored = tuple(iter_route_trees(first))
    assert len(restored) == 2
    assert all(validate_route_tree(tree).valid for tree in restored)
