import pytest

from chemtools.util.rdkit_helpers import rdkit_available


pytestmark = pytest.mark.skipif(not rdkit_available(), reason="RDKit required for reactivity sterics v2 tests")


def test_sterics_v2_o_tolyl_bromide_has_positive_score():
    from chemtools.taxonomy.v1.ar_context_sterics_v2 import analyze_smiles_reactivity_sterics_v2

    res = analyze_smiles_reactivity_sterics_v2("Cc1ccccc1Br")
    assert res.get("error") is None
    assert res["aryl_anchor_site_count"] >= 1
    assert res["aryl_steric_score_0_10_max"] > 0


def test_sterics_v2_bromobenzene_score_zero():
    from chemtools.taxonomy.v1.ar_context_sterics_v2 import analyze_smiles_reactivity_sterics_v2

    res = analyze_smiles_reactivity_sterics_v2("Brc1ccccc1")
    assert res.get("error") is None
    assert res["aryl_anchor_site_count"] >= 1
    assert res["aryl_steric_score_0_10_max"] == 0


def test_sterics_v2_benzaldehyde_detects_ar_cho_and_score_zero():
    from chemtools.taxonomy.v1.ar_context_sterics_v2 import analyze_smiles_reactivity_sterics_v2

    res = analyze_smiles_reactivity_sterics_v2("O=Cc1ccccc1")
    assert res.get("error") is None
    assert res["aryl_anchor_site_count"] >= 1
    assert any(s.get("label") == "Ar-CHO" for s in res.get("sites", []))
    assert res["aryl_steric_score_0_10_max"] == 0


def test_sterics_v2_output_shapes_align_with_site_count():
    from chemtools.taxonomy.v1.ar_context_sterics_v2 import analyze_smiles_reactivity_sterics_v2

    res = analyze_smiles_reactivity_sterics_v2("Cc1ccccc1Br")
    n = res["aryl_anchor_site_count"]
    assert len(res["aryl_ortho_sub_count_list"]) == n
    assert len(res["aryl_ortho_bulk_heavy_list"]) == n
    assert len(res["aryl_ortho_bulk_score_list"]) == n
    assert len(res["aryl_steric_score_0_10_list"]) == n
    assert len(res["sites"]) == n
    assert all(isinstance(pair, list) and len(pair) == 2 for pair in res["aryl_ortho_bulk_heavy_list"])

