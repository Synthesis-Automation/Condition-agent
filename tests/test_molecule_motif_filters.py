import pytest

from chemtools.featurizers import molecule as molecule_mod


def _motif(compound_id: str) -> dict:
    return {"compound_id": compound_id}


def test_motif_family_classification() -> None:
    assert molecule_mod._motif_family("Ar-Br") == "aryl"
    assert molecule_mod._motif_family("AromN-5") == "aryl"
    assert molecule_mod._motif_family("Pyridine-3") == "aryl"
    assert molecule_mod._motif_family("Alkyl-CH3") == "alkyl"
    assert molecule_mod._motif_family("Vinyl-Cl") == "alkyl"
    assert molecule_mod._motif_family("Unknown") is None


def test_background_motifs_filtered_when_non_background_present() -> None:
    motifs = [_motif("Ar-H"), _motif("R-H"), _motif("Ar-Br")]
    filtered = molecule_mod._filter_background_motifs(
        motifs,
        include_h_motifs=False,
        requested_ids=set(),
    )
    assert [m["compound_id"] for m in filtered] == ["Ar-Br"]


def test_background_motifs_kept_when_requested() -> None:
    motifs = [_motif("Ar-H"), _motif("R-H"), _motif("Ar-Br")]
    filtered, requested = molecule_mod._filter_motifs_by_targets(motifs, ["Ar-H"])
    filtered = molecule_mod._filter_background_motifs(
        filtered,
        include_h_motifs=False,
        requested_ids=requested,
    )
    assert [m["compound_id"] for m in filtered] == ["Ar-H"]


def test_target_group_suffix_matching() -> None:
    motifs = [_motif("Ar-Br"), _motif("Ar-Cl")]
    filtered, requested = molecule_mod._filter_motifs_by_targets(motifs, ["Br"])
    assert [m["compound_id"] for m in filtered] == ["Ar-Br"]
    assert requested == {"Ar-Br"}


def test_dedupe_background_motifs() -> None:
    motifs = [_motif("Ar-H"), _motif("Ar-H"), _motif("R-H"), _motif("Ar-Br")]
    deduped = molecule_mod._dedupe_background_motifs(motifs)
    assert [m["compound_id"] for m in deduped] == ["Ar-H", "R-H", "Ar-Br"]
