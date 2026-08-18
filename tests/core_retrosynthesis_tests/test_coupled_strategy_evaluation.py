"""Held-out v1 operator-pair evaluation regressions."""

from __future__ import annotations

import hashlib
from pathlib import Path

from core_retrosynthesis.cli import _parser
from core_retrosynthesis.coupled_strategy_evaluation import (
    CoupledStrategyEvaluationCase,
    CoupledStrategyEvaluationConfig,
    FrozenV1HeldoutPanel,
    PromotedV1OperatorPair,
    evaluate_v1_case,
    load_frozen_v1_heldout_panel,
    run_v1_coupled_strategy_evaluation,
    write_frozen_v1_heldout_panel,
)
from core_retrosynthesis.coupled_strategy_evaluation_review import (
    render_v1_coupled_strategy_evaluation_html,
)
from core_retrosynthesis.generic_compiler import (
    generic_operator_identity_from_observation,
)
from core_retrosynthesis.generic_models import (
    GenericDisconnectionCandidate,
    GenericSearchDiagnostics,
    GenericTemplateLibrary,
)


def _candidate(
    target: str,
    precursors: str,
    operator: str,
    template: str,
    score: float,
) -> GenericDisconnectionCandidate:
    return GenericDisconnectionCandidate(
        target_smiles=target,
        precursor_smiles=precursors,
        proposed_reaction_smiles=f"{precursors}>>{target}",
        transformation_kind=None,
        abstraction_level="L1",
        compiler_engine="test",
        template_id=template,
        score=score,
        context_similarity=0.0,
        product_similarity=1.0,
        precursor_similarity=1.0,
        template_specificity=1.0,
        independent_reference_support=2,
        forward_validation_status="verified",
        center_transition_key="center",
        disconnection_site_key="",
        precedent_reaction_ids=(f"LITERATURE:{template}",),
        operator_id=operator,
    )


def test_operator_identity_reuses_route_core_observation() -> None:
    observation = {
        "edits": [
            {
                "edit_type": "order_changed",
                "atom_1": {
                    "atom_map_number": 1,
                    "element": "C",
                    "formal_charge": 0,
                    "aromatic": False,
                    "hybridization": "SP3",
                },
                "atom_2": {
                    "atom_map_number": 2,
                    "element": "O",
                    "formal_charge": 0,
                    "aromatic": False,
                    "hybridization": "SP2",
                },
                "old_order": "SINGLE",
                "new_order": "DOUBLE",
            }
        ]
    }

    two_part = generic_operator_identity_from_observation(
        "[CH3:1][OH:2]>>[CH2:1]=[O:2]", observation
    )
    three_part = generic_operator_identity_from_observation(
        "[CH3:1][OH:2]>O>[CH2:1]=[O:2]", observation
    )

    assert two_part is not None
    assert three_part == two_part
    assert two_part.operator_id.startswith("OP1:")


def test_promoted_v1_pair_executes_two_validated_steps_with_fallback() -> None:
    strategy = PromotedV1OperatorPair(
        strategy_id="strategy:1",
        relationship_class="handle_progression",
        first_operator_id="op:first",
        second_operator_id="op:second",
        training_patent_ids=("TRAIN1", "TRAIN2"),
        training_occurrence_count=3,
        v2_dependency_counts=(("created_handle_consumed", 3),),
    )
    case = CoupledStrategyEvaluationCase(
        case_id="case:1",
        occurrence_id="occurrence:1",
        strategy_id="strategy:1",
        patent_id="HELDOUT1",
        split="test",
        relationship_class="handle_progression",
        v2_dependency_class="created_handle_consumed",
        target_smiles="CN",
        expected_intermediate_smiles="C=O",
        expected_terminal_precursor_smiles="CO",
        observed_first_reaction_smiles="CO>>C=O",
        observed_second_reaction_smiles="C=O>>CN",
        exact_target_seen_in_training=False,
        target_scaffold_seen_in_training=False,
    )
    second = _candidate("CN", "C=O", "op:second", "second", 0.8)
    first = _candidate("C=O", "CO", "op:first", "first", 0.9)

    def searcher(target, _library, *, operator_ids=(), **_kwargs):
        candidates = (second,) if target == "CN" else (first,) if target == "C=O" else ()
        if operator_ids:
            candidates = tuple(
                item for item in candidates if item.operator_id in operator_ids
            )
        return candidates, GenericSearchDiagnostics(
            validation_attempt_count=len(candidates),
            valid_candidate_count=len(candidates),
        )

    result = evaluate_v1_case(
        case,
        strategy,
        GenericTemplateLibrary((), 0, 0, {}, {}),
        config=CoupledStrategyEvaluationConfig(
            panel_size=1,
            top_k=3,
            max_templates_to_apply=3,
            max_candidates_to_validate=3,
        ),
        searcher=searcher,
    )

    assert result.baseline_intermediate_rank == 1
    assert result.baseline_operator_pair_rank == 1
    assert result.promoted_operator_pair_rank == 1
    assert result.one_step_fallback_preserved
    assert len(result.promoted_actions) == 1
    action = result.promoted_actions[0]
    assert action.first_reaction_smiles == "CO>>C=O"
    assert action.second_reaction_smiles == "C=O>>CN"
    assert action.first_forward_validation_status == "verified"
    assert action.second_forward_validation_status == "verified"

    report = {
        "artifact_type": "v1_coupled_strategy_heldout_evaluation",
        "panel_case_count": 1,
        "metrics": {
            "baseline_pair_recall": 1.0,
            "promoted_pair_recall": 1.0,
            "promoted_validation_attempt_count": 2,
        },
        "results": [result.to_dict()],
    }
    page = render_v1_coupled_strategy_evaluation_html(report)
    assert "two independently forward-validated physical reactions" in page
    assert "HELDOUT1" in page
    assert "Promoted logical actions" in page


def test_v1_promotion_rejects_independent_steps() -> None:
    try:
        PromotedV1OperatorPair(
            strategy_id="strategy:bad",
            relationship_class="independent_sites",
            first_operator_id="one",
            second_operator_id="two",
            training_patent_ids=("P1", "P2"),
            training_occurrence_count=2,
            v2_dependency_counts=(),
        )
    except ValueError as error:
        assert "v1 structural relationships" in str(error)
    else:
        raise AssertionError("independent steps must not become macro actions")


def test_cli_exposes_heldout_v1_evaluation() -> None:
    arguments = _parser().parse_args(
        [
            "evaluate-v1-coupled-strategies",
            "routes.jsonl.gz",
            "operators.json.gz",
            "report.json",
            "report.html",
        ]
    )
    assert arguments.panel_size == 12
    assert arguments.top_k == 5


def test_frozen_panel_evaluates_missing_library_without_reselection(
    tmp_path: Path,
) -> None:
    route_source = tmp_path / "routes.core.jsonl.gz"
    route_source.write_bytes(b"fixed-route-source")
    strategy = PromotedV1OperatorPair(
        strategy_id="strategy:frozen",
        relationship_class="handle_progression",
        first_operator_id="op:first",
        second_operator_id="op:second",
        training_patent_ids=("TRAIN1", "TRAIN2"),
        training_occurrence_count=3,
        v2_dependency_counts=(("created_handle_consumed", 3),),
    )
    case = CoupledStrategyEvaluationCase(
        case_id="case:frozen",
        occurrence_id="occurrence:frozen",
        strategy_id=strategy.strategy_id,
        patent_id="HELDOUT1",
        split="test",
        relationship_class="handle_progression",
        v2_dependency_class="created_handle_consumed",
        target_smiles="CN",
        expected_intermediate_smiles="C=O",
        expected_terminal_precursor_smiles="CO",
        observed_first_reaction_smiles="CO>>C=O",
        observed_second_reaction_smiles="C=O>>CN",
        exact_target_seen_in_training=False,
        target_scaffold_seen_in_training=False,
    )
    config = CoupledStrategyEvaluationConfig(
        panel_size=1,
        top_k=3,
        max_templates_to_apply=3,
        max_candidates_to_validate=3,
    )
    panel = FrozenV1HeldoutPanel(
        panel_id="CRV1PANEL2:test",
        route_core_source=str(route_source),
        route_core_sha256=hashlib.sha256(route_source.read_bytes()).hexdigest(),
        config=config,
        strategies=(strategy,),
        cases=(case,),
        required_strategy_ids=(strategy.strategy_id,),
    )
    panel_path = tmp_path / "panel.json"
    write_frozen_v1_heldout_panel(panel, panel_path)

    loaded = load_frozen_v1_heldout_panel(panel_path)
    report = run_v1_coupled_strategy_evaluation(
        route_source,
        GenericTemplateLibrary((), 0, 0, {}, {}),
        config=None,
        frozen_panel=loaded,
    )

    assert loaded == panel
    assert report["frozen_panel_id"] == panel.panel_id
    assert report["panel_selection"] == "library_independent_frozen"
    assert report["panel_case_count"] == 1
    assert len(report["capability_gaps"]) == 1
    assert report["metrics"]["promoted_pair_hit_count"] == 0


def test_cli_exposes_library_independent_panel() -> None:
    panel = _parser().parse_args(
        [
            "build-v1-coupled-panel",
            "routes.jsonl.gz",
            "panel.json",
            "--include-strategy-id",
            "strategy:required",
        ]
    )
    evaluation = _parser().parse_args(
        [
            "evaluate-v1-coupled-strategies",
            "routes.jsonl.gz",
            "operators.json.gz",
            "report.json",
            "report.html",
            "--frozen-panel",
            "panel.json",
        ]
    )

    assert panel.include_strategy_id == ["strategy:required"]
    assert evaluation.frozen_panel == "panel.json"
