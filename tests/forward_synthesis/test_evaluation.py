"""Product-hidden forward replay evaluation regressions."""

from __future__ import annotations

import pytest

from core_retrosynthesis import build_generic_library
from forward_synthesis import (
    build_forward_library,
    evaluate_product_hidden_replay,
)


def _library():
    generic = build_generic_library(
        (
            {
                "reaction_id": "train-1",
                "reference_id": "training-reference",
                "reaction_smiles": "CCBr.N>>CCN",
            },
        ),
        levels=("L1", "L2"),
        admission_mode="data_driven",
    )
    return build_forward_library(generic)


def test_product_hidden_replay_reports_exact_rank() -> None:
    report = evaluate_product_hidden_replay(
        (
            {
                "reaction_id": "test-1",
                "reference_id": "held-out-reference",
                "reaction_smiles": "CCBr.N>>CCN",
            },
        ),
        _library(),
    )

    assert report.case_count == 1
    assert report.exact_top_1_count == 1
    assert report.exact_any_count == 1
    assert report.reciprocal_rank_mean == 1.0


def test_product_hidden_replay_rejects_reference_leakage() -> None:
    with pytest.raises(ValueError, match="overlaps"):
        evaluate_product_hidden_replay(
            (
                {
                    "reaction_id": "leaked-1",
                    "reference_id": "training-reference",
                    "reaction_smiles": "CCBr.N>>CCN",
                },
            ),
            _library(),
        )
