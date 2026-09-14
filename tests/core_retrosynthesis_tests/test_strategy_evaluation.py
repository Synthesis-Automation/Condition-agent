"""Strategy panels freeze split membership before any operator search."""

import json
from types import SimpleNamespace

import pytest

from core_retrosynthesis import strategy_evaluation


def test_freeze_strategy_panel_preserves_grouping_and_refuses_overwrite(
    tmp_path, monkeypatch
):
    records = [
        SimpleNamespace(
            observation_id=f"observation-{i}",
            reaction_id=f"reaction-{i}",
            canonical_reaction_id=f"canonical-{i}",
            reference_id=f"reference-{i}",
            reaction_core={},
            scaffold_tokens=(),
            publication_year=None,
            recipe_id="recipe",
            source_dataset="fixture",
            reaction_smiles="C>>C",
            transformation_class="fixture",
        )
        for i in range(20)
    ]
    records[1].reference_id = records[0].reference_id
    index = SimpleNamespace(
        rows=records, select=lambda positions: [records[i] for i in positions]
    )
    monkeypatch.setattr(strategy_evaluation, "load_generic_index", lambda path: index)
    monkeypatch.setattr(
        strategy_evaluation,
        "build_generic_library",
        lambda *args, **kwargs: pytest.fail("freeze must not compile operators"),
    )
    path = strategy_evaluation.freeze_strategy_panel(
        "source", tmp_path, size=20, query_limit=3
    )
    panel = json.loads(path.read_text())
    train = {row["observation_id"] for row in panel["train"]}
    test = {row["observation_id"] for row in panel["test"]}
    assert not train & test
    assert set(panel["query_ids"]) <= test
    assert len(panel["query_ids"]) <= 3
    assert panel["leakage"]["reference_overlap_count"] == 0
    assert panel["leakage"]["canonical_reaction_overlap_count"] == 0
    assert ("observation-0" in train) == ("observation-1" in train)
    with pytest.raises(FileExistsError):
        strategy_evaluation.freeze_strategy_panel("source", tmp_path, size=20)


def test_source_compilation_failure_does_not_suppress_valid_product_queries(
    tmp_path, monkeypatch,
) -> None:
    panel = tmp_path / "panel.json"
    panel.write_text(json.dumps({
        "definition_id": "fixture", "train": [], "leakage": {},
        "query_ids": ["unresolved", "invalid"],
        "test": [
            {"observation_id": key, "reaction_smiles": reaction,
             "transformation_class": ""}
            for key, reaction in (("unresolved", "C>O>CO"), ("invalid", "C>>invalid"))
        ],
    }), encoding="utf-8")
    monkeypatch.setattr(
        strategy_evaluation, "build_generic_library",
        lambda *args, **kwargs: SimpleNamespace(
            source_row_count=0, accepted_observation_count=0,
        ),
    )
    monkeypatch.setattr(strategy_evaluation, "save_generic_library", lambda *args: None)
    monkeypatch.setattr(
        strategy_evaluation, "compile_generic_templates",
        lambda *args, **kwargs: SimpleNamespace(
            templates=(), rejection_reason="unresolved_correspondence",
        ),
    )
    candidate = SimpleNamespace(
        strategy_id="STRAT1:test", precursor_smiles="C=O",
        forward_validation_status="verified_signature", to_dict=lambda: {},
    )
    diagnostics = SimpleNamespace(validation_attempt_count=1, to_dict=lambda: {})
    queried = []

    def flat(target, *args, **kwargs):
        queried.append(("flat", target))
        return (candidate,), diagnostics

    def grouped(target, *args, **kwargs):
        queried.append(("strategy", target))
        return SimpleNamespace(
            strategies=(SimpleNamespace(realizations=(candidate,)),),
            diagnostics=diagnostics, to_dict=lambda: {},
        )

    monkeypatch.setattr(strategy_evaluation, "disconnect_operator_ladder_detailed", flat)
    monkeypatch.setattr(strategy_evaluation, "disconnect_strategies_detailed", grouped)
    report = strategy_evaluation.evaluate_strategy_panel(panel)
    assert queried == [("flat", "CO"), ("strategy", "CO")]
    assert report["query_count"] == 2
    assert report["source_compilation_failures"] == 2
    assert report["target_query_failures"] == 1
    for metrics in report["summary"].values():
        assert metrics["coverage"] == 0.5
        assert metrics["strategy_recovered"] is None
        assert metrics["recovery_denominators"]["strategy_recovered"] == 0
    with pytest.raises(FileExistsError):
        strategy_evaluation.evaluate_strategy_panel(panel)
