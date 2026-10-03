"""Material identity preserves chemical form and descriptor provenance."""

import pytest

from reactive_taxonomy.material_identity import identify_material


def test_material_order_and_atom_maps_do_not_change_identity():
    ordinary = identify_material("CCO")
    mapped = identify_material("[OH:3][CH2:2][CH3:1]")
    assert ordinary.canonical_smiles == mapped.canonical_smiles == "CCO"
    assert ordinary.molecular_weight == mapped.molecular_weight == pytest.approx(46.069)
    assert ordinary.molecular_weight_unit == "g/mol"
    assert ordinary.rdkit_version
    assert ordinary.to_dict()["schema_version"] == "material_identity.v1"


@pytest.mark.parametrize("left,right", [
    ("C[C@H](O)F", "C[C@@H](O)F"), ("C[C@H](O)F", "CC(O)F"),
    ("[13CH3]CO", "CCO"), ("C[NH3+]", "CN"),
    ("[Na+].[O-]C(=O)C", "CC(=O)O"), ("CC=O", "C=CO"),
    ("C/C=C/C", "C/C=C\\C"),
])
def test_material_identity_does_not_relax_chemistry(left, right):
    assert identify_material(left).canonical_smiles != identify_material(right).canonical_smiles


def test_complete_salt_form_and_unspecified_stereo_are_visible():
    salt = identify_material("[Na+].[O-]C(=O)C")
    assert salt.component_count == 2
    assert salt.molecular_weight == pytest.approx(82.034, abs=0.01)
    assert identify_material("CC(=O)[O-].[Na+]").canonical_smiles == salt.canonical_smiles
    assert salt.warnings
    assert "Unspecified stereochemistry" in identify_material("CC(O)F").warnings[0]


@pytest.mark.parametrize("smiles", [None, "", " ", "invalid", "C(C)(C)(C)(C)C", "*", "[CH3]"])
def test_invalid_or_unsupported_material_cannot_receive_a_mass_stop(smiles):
    with pytest.raises(ValueError):
        identify_material(smiles)
