from chemtools.featurizers.analysis.reactions import canonical_family_label


def test_legacy_aliases_resolve_via_taxonomy() -> None:
    assert canonical_family_label("Suzuki_CC") == "suzuki_miyaura"
    assert canonical_family_label("Buchwald_CN") == "c_n_cross_coupling"
    assert canonical_family_label("SNAr-CO") == "c_o_coupling"
    assert canonical_family_label("SNAr-CS") == "c_s_coupling"
    assert canonical_family_label("SNAr-CN") == "c_n_cross_coupling"
    assert canonical_family_label("Fischer_esterification") == "esterification"
    assert canonical_family_label("RCM") == "ring_closing_metathesis"
