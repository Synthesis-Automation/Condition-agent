import pytest

molpipeline = pytest.importorskip("molpipeline")

from chemtools import precedent


def test_attach_molpipeline_features_on_pack():
    pack = {
        "precedents": [
            {
                "reagents": [
                    {"role": "ligand", "smiles": "c1ccccc1"},
                    {"role": "base", "smiles": "CCO"},
                ],
            }
        ]
    }
    cfg = {
        "include_role_features": True,
        "include_concat": True,
        "suppress_errors": False,
        "query_role_smiles": {
            "ligand": ["c1ccccc1"],
            "base": ["CCO"],
        },
    }
    updated = precedent._attach_molpipeline_features(pack, cfg)
    first = updated["precedents"][0]
    assert "molpipeline_role_features" in first
    assert "molpipeline_feature_vector" in first
    assert isinstance(first["molpipeline_feature_vector"], list)
    assert "molpipeline_query_vector" in updated
    assert isinstance(updated["molpipeline_query_vector"], list)
