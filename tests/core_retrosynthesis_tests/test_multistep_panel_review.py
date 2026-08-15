"""Chemist-facing multistep comparison panel regressions."""

from __future__ import annotations

from pathlib import Path

from core_retrosynthesis.multistep import plan_multistep_routes
from core_retrosynthesis.multistep_panel_review import (
    MultistepPanelCase,
    MultistepPanelTarget,
    load_multistep_panel_targets,
    render_multistep_panel_html,
    write_multistep_panel_artifacts,
)

from .test_multistep import _LiteratureIndex, _candidate, _expander


def _case() -> MultistepPanelCase:
    target = MultistepPanelTarget(
        target_id="octane",
        name="Octane",
        category="test",
        smiles="CCCCCCCC",
    )
    result = plan_multistep_routes(
        target.smiles,
        object(),
        _LiteratureIndex(),
        max_depth=2,
        molecular_weight_threshold=100.0,
        expander=_expander(
            {target.smiles: (_candidate(target.smiles, "CCCC.CCCC"),)}
        ),
    )
    return MultistepPanelCase(target=target, baseline=result, policy=result)


def test_familiar_target_definition_is_valid_and_diverse() -> None:
    targets = load_multistep_panel_targets()

    assert len(targets) == 12
    assert len({target.target_id for target in targets}) == 12
    assert {target.category for target in targets} == {
        "amide",
        "biaryl",
        "ester",
        "ether",
        "sulfonamide",
    }


def test_route_policy_test_panel_is_valid_and_fixed() -> None:
    definition = (
        Path(__file__).resolve().parents[2]
        / "core_retrosynthesis"
        / "definitions"
        / "multistep_route_policy_test_panel.v1.json"
    )

    targets = load_multistep_panel_targets(definition)

    assert len(targets) == 8
    assert len({target.target_id for target in targets}) == 8
    assert {target.category for target in targets} == {
        "heldout_observed_3_step"
    }


def test_route_policy_validation_panel_is_valid_and_fixed() -> None:
    definition = (
        Path(__file__).resolve().parents[2]
        / "core_retrosynthesis"
        / "definitions"
        / "multistep_route_policy_validation_panel.v1.json"
    )

    targets = load_multistep_panel_targets(definition)

    assert len(targets) == 8
    assert len({target.target_id for target in targets}) == 8
    assert {target.category for target in targets} == {
        "validation_observed_3_step"
    }


def test_multistep_panel_html_contains_routes_and_review_controls() -> None:
    base_case = _case()
    case = MultistepPanelCase(
        target=base_case.target,
        baseline=base_case.baseline,
        policy=base_case.policy,
        reference_route_id="observed-route",
        reference_reactions=("CCCC.CCCC>>CCCCCCCC",),
        reference_maximum_depth=1,
    )
    document = render_multistep_panel_html((case,), top_k=1)

    assert "Octane" in document
    assert "CCCC.CCCC&gt;&gt;CCCCCCCC" not in document
    assert "Same ranked routes" in document
    assert "Export review JSON" in document
    assert "Observed precedent" in document
    assert "observed-route" in document
    assert "<svg" in document


def test_multistep_panel_artifacts_are_written(tmp_path: Path) -> None:
    output_html = tmp_path / "panel.html"
    output_json = tmp_path / "panel.json"

    report = write_multistep_panel_artifacts(
        (_case(),),
        output_html,
        output_json=output_json,
        top_k=1,
    )

    assert report["target_count"] == 1
    assert report["changed_ranking_count"] == 0
    assert output_html.is_file()
    assert output_json.is_file()


def test_reference_only_panel_does_not_claim_policy_ranking_change() -> None:
    base_case = _case()
    case = MultistepPanelCase(
        target=base_case.target,
        baseline=base_case.baseline,
        reference_route_id="observed-route",
        reference_reactions=("CCCC.CCCC>>CCCCCCCC",),
    )

    document = render_multistep_panel_html((case,), top_k=1)

    assert "Observed vs planner" in document
    assert "Ranking changed" not in document
