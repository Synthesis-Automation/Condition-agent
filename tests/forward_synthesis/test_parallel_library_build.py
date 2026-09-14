"""Parallel source verification preserves admission, support and JSON identity."""

from dataclasses import replace

import pytest

from core_retrosynthesis.generic_library import build_generic_library
from forward_synthesis.library import build_forward_library


def test_parallel_forward_build_matches_serial_with_failed_and_duplicate_sources() -> None:
    generic = build_generic_library(
        [
            {"reaction_id": str(i), "reference_id": str(i), "reaction_smiles": reaction}
            for i, reaction in enumerate(("CCBr.N>>CCN", "CC=O>>CCO"))
        ],
        levels=("L1", "L2"), admission_mode="data_driven",
    )
    first = generic.templates[0]
    failing_source = replace(
        first, template_id="failing-source",
        precedents=(replace(first.precedents[0], product_smiles="C#N"),),
    )
    malformed = replace(first, precursor_smarts="invalid", template_id="malformed")
    duplicate = replace(
        first, observation_support=17, independent_reference_support=9,
        precedents=(replace(first.precedents[0], reaction_id="duplicate", reference_id="other"),),
    )
    source = {
        "templates": (*generic.templates, failing_source, malformed, duplicate),
        "definition": generic.definition,
    }
    serial = build_forward_library(source)
    parallel = build_forward_library(source, workers=2)
    assert serial.admitted_operator_count > 0
    assert serial.rejection_counts == {
        "invalid_operator_contract": 1, "source_forward_round_trip_failed": 1,
    }
    assert serial.to_dict() == parallel.to_dict()
    assert any(operator.observation_support == 17 for operator in serial.operators)


@pytest.mark.parametrize("workers", [0, -1])
def test_forward_build_rejects_nonpositive_worker_count(workers) -> None:
    with pytest.raises(ValueError, match="workers must be positive"):
        build_forward_library({"templates": []}, workers=workers)
