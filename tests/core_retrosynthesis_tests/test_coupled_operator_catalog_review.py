"""Graphical two-step operator-pair catalog regressions."""

from __future__ import annotations

from pathlib import Path

from core_retrosynthesis.coupled_operator_catalog_review import (
    render_v1_operator_pair_catalog_html,
    write_v1_operator_pair_catalog_html,
)
from core_retrosynthesis.coupled_strategy_evaluation import (
    CoupledStrategyEvaluationCase,
    PromotedV1OperatorPair,
)
from core_retrosynthesis.generic_models import (
    GenericGraphOperator,
    GenericTemplateLibrary,
)


def _strategy() -> PromotedV1OperatorPair:
    return PromotedV1OperatorPair(
        strategy_id="CRV1OP1:test",
        relationship_class="handle_progression",
        first_operator_id="OP1:first",
        second_operator_id="OP1:second",
        training_patent_ids=("TRAIN1", "TRAIN2"),
        training_occurrence_count=3,
        v2_dependency_counts=(("created_handle_consumed", 3),),
    )


def _case() -> CoupledStrategyEvaluationCase:
    return CoupledStrategyEvaluationCase(
        case_id="CRV1CASE1:test",
        occurrence_id="occurrence:test",
        strategy_id="CRV1OP1:test",
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


def _operator(operator_id: str, edit: str) -> GenericGraphOperator:
    return GenericGraphOperator(
        operator_id=operator_id,
        operator_signature=f"signature:{operator_id}",
        edit_tokens=(edit,),
        realization_ids=(f"realization:{operator_id}",),
        abstraction_levels=("L0", "L1"),
        observation_support=3,
        independent_reference_support=2,
    )


def _library() -> GenericTemplateLibrary:
    return GenericTemplateLibrary(
        templates=(),
        source_row_count=2,
        accepted_observation_count=2,
        rejection_counts={},
        definition={},
        operators=(
            _operator("OP1:first", "order_changed:C-O:SINGLE>DOUBLE"),
            _operator("OP1:second", "formed:C-N:NONE>SINGLE"),
        ),
    )


def test_pair_catalog_renders_both_physical_reaction_graphics() -> None:
    page = render_v1_operator_pair_catalog_html(
        (_strategy(),),
        (_case(),),
        _library(),
        capability_gap_count=3,
    )

    assert page.count('<article class="pair-card"') == 1
    assert page.count("<svg") == 2
    assert "Physical step 1" in page
    assert "Physical step 2" in page
    assert "C=O" in page
    assert "OP1:first" in page
    assert "OP1:second" in page
    assert "Handle Progression" in page
    assert "created_handle_consumed" in page
    assert "3</strong><small>coverage gaps" in page


def test_pair_catalog_writer_reports_pair_and_graphic_counts(
    tmp_path: Path,
) -> None:
    output = tmp_path / "pairs.html"

    summary = write_v1_operator_pair_catalog_html(
        (_strategy(),),
        (_case(),),
        _library(),
        output,
    )

    assert output.is_file()
    assert summary["rendered_pair_count"] == 1
    assert summary["rendered_reaction_graphic_count"] == 2
    assert summary["html_bytes"] == output.stat().st_size
