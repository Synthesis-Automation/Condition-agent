"""Graph audit scope, invalid structures, material forms and stereo uncertainty."""

import pytest

from reactive_taxonomy.structure_audit import audit_structure


def test_audit_preserves_salt_isotopes_and_specified_stereo() -> None:
    audit = audit_structure("[13CH3][C@H](O)C.[Na+].[Cl-]")
    assert audit.valid
    assert audit.component_count == 3
    assert audit.formal_charge == 0
    assert "13C" in audit.formula
    assert "@" in audit.canonical_smiles
    assert audit.unspecified_stereo_count == 0


def test_map_labels_do_not_change_identity_or_invent_correspondence() -> None:
    mapped = audit_structure("[CH3:2][OH:8]")
    assert mapped.canonical_smiles == audit_structure("CO").canonical_smiles
    assert mapped.formula == "CH4O"
    assert "ATOM_MAP_LABELS_NOT_REACTION_CORRESPONDENCE" in mapped.warnings


@pytest.mark.parametrize("smiles", ["C1", "C(C)(C)(C)(C)C", "*CC", "[CH3]"])
def test_invalid_or_unresolved_graphs_are_retained(smiles: str) -> None:
    audit = audit_structure(smiles)
    assert audit.valid is False
    assert audit.input_smiles == smiles
    assert audit.warnings


def test_unspecified_stereo_and_ring_sizes_are_observations() -> None:
    audit = audit_structure("CC(O)C1CCCCC1")
    assert audit.valid
    assert audit.unspecified_stereo_count == 1
    assert audit.ring_sizes == (6,)
    assert "UNSPECIFIED_STEREOCHEMISTRY" in audit.warnings
