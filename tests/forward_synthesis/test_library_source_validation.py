"""Prepared forward operators stay compatible with their pinned source data."""

from dataclasses import replace

import pytest

from core_retrosynthesis import build_generic_library
from forward_synthesis import build_forward_library, validate_forward_library_source
from forward_synthesis import library as library_module


@pytest.fixture(scope="module")
def libraries():
    generic = build_generic_library(
        [{"reaction_id": "r1", "reference_id": "ref1", "reaction_smiles": "CCBr.N>>CCN"}],
        levels=("L1", "L2"), admission_mode="data_driven",
    )
    return generic, build_forward_library(generic)


def test_source_validation_accepts_objects_and_structural_mappings_without_chemistry(
    libraries, monkeypatch,
) -> None:
    generic, forward = libraries

    def forbidden(*args, **kwargs):
        raise AssertionError("source compatibility must not rebuild or execute chemistry")

    monkeypatch.setattr(library_module, "build_forward_library", forbidden)
    monkeypatch.setattr(library_module, "_source_round_trip", forbidden)
    monkeypatch.setattr(library_module.rdChemReactions, "ReactionFromSmarts", forbidden)
    validate_forward_library_source(forward, generic)
    validate_forward_library_source(forward, generic.to_dict())


@pytest.mark.parametrize("change, message", [
    ({"source_library_definition_id": "different.v1"}, "source definition"),
    ({"source_template_count": 100}, "template count"),
    ({"definition_id": "forward_operator_library.v999"}, "schema or definition"),
])
def test_source_validation_rejects_mismatched_metadata(libraries, change, message) -> None:
    generic, forward = libraries
    with pytest.raises(ValueError, match=message):
        validate_forward_library_source(replace(forward, **change), generic)


@pytest.mark.parametrize("change", [
    {"operator_id": "OP1:other"},
    {"realization_id": "REAL2:other"},
    {"template_id": "other"},
    {"edit_tokens": ("invented edit",)},
    {"operator_signature": "invented signature"},
    {"stereo_policy": "relaxed"},
    {"observation_support": 100},
    {"independent_reference_support": 100},
    {"named_annotations": ("invented annotation",)},
])
def test_source_validation_rejects_changed_operator_content(libraries, change) -> None:
    generic, forward = libraries
    original = forward.operators[0]
    assert any(getattr(original, key) != value for key, value in change.items())
    altered = replace(original, **change)
    with pytest.raises(ValueError, match="forward operator"):
        validate_forward_library_source(
            replace(forward, operators=(altered, *forward.operators[1:])), generic,
        )


def test_source_validation_rejects_changed_smarts_and_provenance(libraries) -> None:
    generic, forward = libraries
    original = forward.operators[0]
    altered_smarts = replace(
        original, precursor_smarts="[C:1]", product_smarts="[C:1]",
        forward_smarts="[C:1]>>[C:1]", reverse_smarts="[C:1]>>[C:1]",
    )
    altered_precedent = replace(
        original, precedents=(replace(original.precedents[0], reference_id="invented"),),
    )
    for altered in (altered_smarts, altered_precedent):
        with pytest.raises(ValueError, match="forward operator"):
            validate_forward_library_source(
                replace(forward, operators=(altered, *forward.operators[1:])), generic,
            )


def test_source_validation_accepts_merged_sources_and_rejected_duplicate(libraries) -> None:
    generic, _ = libraries
    first = generic.templates[0]
    duplicate = replace(
        first, template_id="duplicate", observation_support=17,
        independent_reference_support=9, named_annotations=("additional source label",),
        precedents=(replace(first.precedents[0], reaction_id="r2", reference_id="ref2"),),
    )
    rejected = replace(
        first, template_id="failed-round-trip",
        precedents=(replace(first.precedents[0], product_smiles="C#N"),),
    )
    source = {
        "templates": (*generic.templates, duplicate, rejected),
        "definition": generic.definition,
    }
    forward = build_forward_library(source)
    assert forward.rejection_counts["source_forward_round_trip_failed"] == 1
    assert any(operator.observation_support == 17 for operator in forward.operators)
    validate_forward_library_source(forward, source)

    # Removing a merged precedent from the pinned source must not pass validation.
    changed = {**source, "templates": (*generic.templates, first, rejected)}
    with pytest.raises(ValueError, match="provenance or support"):
        validate_forward_library_source(forward, changed)


def test_source_validation_rejects_missing_source_identity(libraries) -> None:
    generic, forward = libraries
    with pytest.raises(ValueError, match="source definition"):
        validate_forward_library_source(forward, {"templates": generic.templates})
