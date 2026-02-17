from chemtools.featurizers.formatters.utils import normalize_motif_id


def test_normalize_motif_id_maps_bpin_alias_to_boronate_ester() -> None:
    assert normalize_motif_id("Ar-Bpin") == "Ar-B(OR)2"
    assert normalize_motif_id("HeteroAr-Bpin") == "HeteroAr-B(OR)2"
    assert normalize_motif_id("Bpin") == "B(OR)2"
