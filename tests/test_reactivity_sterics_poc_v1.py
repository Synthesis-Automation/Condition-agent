import pytest

from pathlib import Path

from chemtools.util.rdkit_helpers import rdkit_available


pytestmark = pytest.mark.skipif(not rdkit_available(), reason="RDKit required for steric reactivity POC v1 tests")


def test_reactivity_sterics_spec_bromotoluene_ortho1():
    from chemtools.taxonomy.v1.reactivity_sterics_computed_v1 import compute_reactivity_sterics_poc_v1

    res = compute_reactivity_sterics_poc_v1("Cc1ccccc1Br")  # o-bromotoluene
    assert res.get("error") is None
    assert res["aryl_anchor_site_count"] == 1
    assert res["aryl_ortho_sub_count_list"] == [1]
    assert res["aryl_ortho_sub_count_max"] == 1
    assert res["aryl_ortho_sub_count_sum"] == 1
    assert res["sites"][0]["label"] == "Ar-Br"


def test_reactivity_sterics_spec_para_bromotoluene_ortho0():
    from chemtools.taxonomy.v1.reactivity_sterics_computed_v1 import compute_reactivity_sterics_poc_v1

    res = compute_reactivity_sterics_poc_v1("Cc1ccc(Br)cc1")  # p-bromotoluene
    assert res.get("error") is None
    assert res["aryl_anchor_site_count"] == 1
    assert res["aryl_ortho_sub_count_list"] == [0]


def test_reactivity_sterics_spec_aromatic_aldehyde_ortho0():
    from chemtools.taxonomy.v1.reactivity_sterics_computed_v1 import compute_reactivity_sterics_poc_v1

    res = compute_reactivity_sterics_poc_v1("O=Cc1ccccc1")  # benzaldehyde
    assert res.get("error") is None
    assert res["aryl_anchor_site_count"] == 1
    assert res["aryl_ortho_sub_count_list"] == [0]
    assert res["sites"][0]["label"] == "Ar-CHO"


def test_reactivity_sterics_direct_matches_spec_for_key_fields():
    from chemtools.taxonomy.v1.ar_context_sterics_v1 import analyze_smiles_ortho_counts
    from chemtools.taxonomy.v1.reactivity_sterics_computed_v1 import compute_reactivity_sterics_poc_v1

    smiles = "Brc1ccc(cc1)C"

    direct = analyze_smiles_ortho_counts(
        smiles,
        Path("calculable_features.compiled.v1.json"),
        Path("organic_groups.v1.json"),
    )
    spec = compute_reactivity_sterics_poc_v1(smiles)

    assert direct.get("error") is None
    assert spec.get("error") is None
    for key in (
        "aryl_anchor_site_count",
        "aryl_ortho_sub_count_list",
        "aryl_ortho_sub_count_max",
        "aryl_ortho_sub_count_sum",
    ):
        assert direct[key] == spec[key]


def test_reactivity_sterics_invalid_smiles_returns_error():
    from chemtools.taxonomy.v1.reactivity_sterics_computed_v1 import compute_reactivity_sterics_poc_v1

    res = compute_reactivity_sterics_poc_v1("not-a-smiles")
    assert res.get("error")

