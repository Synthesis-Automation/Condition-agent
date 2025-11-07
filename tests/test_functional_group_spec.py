"""
Lightweight integrity tests for functional group definitions sourced from
chemtools/featurizers/calculable_features.json.
"""

import pytest

from chemtools.util import functional_groups
from chemtools.util.rdkit_helpers import rdkit_available


def test_functional_group_definitions_are_unique_and_well_formed() -> None:
    """Ensure every functional group entry has valid metadata."""
    definitions = list(functional_groups.iter_group_definitions())
    assert definitions, "functional_groups section should not be empty"

    names = [definition.name for definition in definitions]
    assert len(names) == len(set(names)), "functional group names must be unique"

    for definition in definitions:
        assert definition.smarts, f"{definition.name} must define at least one SMARTS"
        for smarts in definition.smarts:
            assert isinstance(smarts, str) and smarts.strip(), f"{definition.name} has invalid SMARTS"
        # Category tags are optional but must be strings when present
        for tag in definition.category_tags:
            assert isinstance(tag, str) and tag.strip(), f"{definition.name} has invalid category tag"


@pytest.mark.skipif(not rdkit_available(), reason="RDKit not available")
def test_category_summary_uses_spec_metadata() -> None:
    """Round-trip a SMILES string to ensure category metadata is respected."""
    categories = functional_groups.get_group_categories("CC(=O)O")
    assert categories["oxygen"], "oxygen category should list detected groups"
    assert "carboxylic_acid" in categories["oxygen"]
