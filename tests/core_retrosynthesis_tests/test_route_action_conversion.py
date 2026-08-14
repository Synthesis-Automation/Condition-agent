"""Route-action replay corpus conversion regressions."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

from reactive_taxonomy import featurize_reaction

from core_retrosynthesis.generic_library import (
    build_generic_library,
    save_generic_library,
)
from core_retrosynthesis.route_action_conversion import (
    convert_route_action_corpus,
    iter_route_action_evaluations,
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
    assert first_report["metrics"]["eligible_step_count"] == 1
    assert first_report["metrics"]["recall"]["exact_precursor"]["top1"] == 1.0
    restored = tuple(iter_route_action_evaluations(first))
    assert len(restored) == 1
    assert restored[0].steps[0].exact_precursor_rank == 1
