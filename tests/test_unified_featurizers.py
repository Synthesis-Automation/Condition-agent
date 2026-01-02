import pytest

from chemtools.featurizers.unified import featurize_molecule, featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


@pytest.mark.skipif(not rdkit_available(), reason="RDKit not available")
def test_unified_molecule_bundle_shape():
    result = featurize_molecule("c1ccccc1Br")
    assert result["kind"] == "molecule"
    assert result["schema_version"] == "v1"
    molecule = result["molecule"]
    assert molecule["smiles"]
    assert "rdkit_props" in molecule
    assert molecule["rdkit_props"]["molecular_weight"] > 0
    assert isinstance(molecule["functional_groups"], list)


@pytest.mark.skipif(not rdkit_available(), reason="RDKit not available")
def test_unified_reaction_bundle_counts():
    result = featurize_reaction("c1ccccc1Br.Oc1ccccc1>>")
    reaction = result["reaction"]
    assert reaction["reactants"]
    assert reaction["aggregates"]["reactant_count"] == len(reaction["reactants"])
    assert "reaction_type" in reaction
    assert reaction["agent_roles"] is not None
    assert reaction["agent_roles"]["agent_count"] == 0
    assert "flags" in reaction["agent_roles"]
