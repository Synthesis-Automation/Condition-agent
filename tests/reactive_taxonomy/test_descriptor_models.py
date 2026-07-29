import pytest

from reactive_taxonomy.descriptors import (
    DescriptorEvidence,
    ElectronicContribution,
    StericContribution,
)
from reactive_taxonomy.descriptors.registry import (
    descriptor_definition_versions,
    load_reactivity_descriptor_definitions,
)


def test_descriptor_evidence_validates_confidence() -> None:
    with pytest.raises(ValueError, match="confidence"):
        DescriptorEvidence(source="graph", method="test", confidence=1.01)


def test_contributions_validate_normalized_ranges() -> None:
    with pytest.raises(ValueError, match="between 0 and 1"):
        StericContribution(
            origin_atom_index=0,
            relation="branch",
            heavy_atom_count=1,
            branch_count=0,
            score=1.1,
        )
    with pytest.raises(ValueError, match="between -1 and 1"):
        ElectronicContribution(
            source_id="test",
            effect="withdrawing",
            pathway="inductive",
            positional_relation="alpha",
            contribution=1.1,
        )


def test_descriptor_definitions_are_validated_and_content_versioned() -> None:
    definitions = load_reactivity_descriptor_definitions()
    assert set(definitions) == {
        "aromatic_systems.v1.json",
        "reactivity_descriptor_rules.v1.json",
        "reactivity_rendering.v1.json",
    }
    versions = dict(descriptor_definition_versions())
    assert set(versions) == set(definitions)
    assert all("@sha256:" in value for value in versions.values())
