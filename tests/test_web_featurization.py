"""Tests for the Qt-free web featurization presentation service."""

from __future__ import annotations

from app.web_api.features import analyze_features, detect_input_kind


MAPPED_INTERCOMPONENT_ANNULATION = (
    "CS(=O)(=O)O[CH2:1][CH2:2][c:3]1[c:9](Cl)"
    "[cH:8][cH:7][cH:6][cH:4]1.[NH2:28][CH3:26]>>"
    "[CH3:26][N:28]1[CH2:1][CH2:2][c:3]2[c:9]1"
    "[cH:8][cH:7][cH:6][cH:4]2"
)
MAPPED_ACYL_TETRAZOLE_REARRANGEMENT = (
    "[CH3:1][CH:2]=[O:3]."
    "[CH3:4][CH2:5][O:6][C:7]([n:9]1[n:13][n:12][n:11][cH:10]1)"
    "=[O:8]>>"
    "[CH3:1][CH2:2][O:3][C:5]([O:6][CH:7]"
    "([n:9]1[n:13][n:12][n:11][cH:10]1)[CH3:4])=[O:8]"
)
MAPPED_CHANGED_CONTINUITY_MACROCYCLIZATION = (
    "CC1(C(C)(C)OB(B2OC(C)(C)C(C)(C)O2)O1)C."
    "[CH3:1][CH:2]([c:11]1[c:18]([C:19]([N:21]([CH2:23]"
    "[c:24]2[n:29]([CH2:30][CH2:31][O:32][Si:33]([C:36]"
    "([CH3:39])([CH3:38])[CH3:37])([CH3:35])[CH3:34])"
    "[n:28][c:26]([CH3:27])[c:25]2Br)[CH3:22])=[O:20])"
    "[cH:17][cH:16][c:14]([F:15])[c:12]1[Cl:13])[O:3]"
    "[c:4]3[c:9]([NH2:10])[n:8][cH:7][c:6](Br)[n:5]3>>"
    "[CH3:1][CH:2]1[c:11]([c:18]2[C:19](=[O:20])[N:21]"
    "([CH3:22])[CH2:23][c:24]([c:25]3[c:6]4[n:5][c:4]"
    "([c:9]([NH2:10])[n:8][cH:7]4)[O:3]1)[n:29]([CH2:30]"
    "[CH2:31][O:32][Si:33]([C:36]([CH3:39])([CH3:38])"
    "[CH3:37])([CH3:35])[CH3:34])[n:28][c:26]3[CH3:27])"
    "[c:12]([Cl:13])[c:14]([F:15])[cH:16][cH:17]2"
)


def test_detect_input_kind_supports_molecules_and_reactions() -> None:
    assert detect_input_kind("Brc1ccccc1") == "molecule"
    assert detect_input_kind("Brc1ccccc1.B(O)O>>c1ccccc1") == "reaction"


def test_molecule_feature_analysis_returns_compact_and_full_results() -> None:
    result = analyze_features("Brc1ccc(N)cc1C#N")

    assert result["input_kind"] == "molecule"
    assert result["valid"] is True
    assert result["overview"]["atom_count"] == 10
    assert result["overview"]["motif_count"] == len(result["motifs"])
    assert result["analysis"]["structure"]["canonical_smiles"]


def test_reaction_feature_analysis_includes_core_projection() -> None:
    result = analyze_features(
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )

    assert result["input_kind"] == "reaction"
    assert result["valid"] is True
    assert result["reaction_core"] is not None
    assert result["reaction_core"]["event_count"] >= 1
    assert result["core_projection"] is not None
    assert "<svg" in result["core_graphic_svg"]
    assert result["reactive_sites"]
    assert all(site["side"] == "reactant" for site in result["reactive_sites"])


def test_web_core_graphic_preserves_intercomponent_annulation_ring() -> None:
    result = analyze_features(MAPPED_INTERCOMPONENT_ANNULATION)

    assert result["valid"] is True
    assert result["core_graphic_error"] is None
    assert "<svg" in result["core_graphic_svg"]
    assert result["core_projection"]["reaction_scope"] == "intermolecular"
    assert result["core_projection"]["minimum_reaction_smiles"] == (
        "*COS(C)(=O)=O.*c1ccccc1Cl.*N>>*CN(*)c1ccccc1*"
    )
    assert result["core_graphic_svg"].count("stroke-dasharray") == 2


def test_web_core_graphic_completes_attached_aromatic_system() -> None:
    result = analyze_features(MAPPED_ACYL_TETRAZOLE_REARRANGEMENT)

    assert result["valid"] is True
    assert result["core_graphic_error"] is None
    assert "<svg" in result["core_graphic_svg"]
    assert result["core_projection"]["minimum_reaction_smiles"] == (
        "*C=O.CCOC(=O)n1cnnn1>>*COC(=O)OC(C)n1cnnn1"
    )


def test_web_core_graphic_connects_changed_macrocycle_tether() -> None:
    result = analyze_features(MAPPED_CHANGED_CONTINUITY_MACROCYCLIZATION)

    assert result["valid"] is True
    assert result["core_graphic_error"] is None
    assert result["core_projection"]["reaction_scope"] == "intramolecular"
    assert result["core_projection"]["annotation"] == (
        "Intramolecular; forms a 12-membered ring"
    )
    assert result["core_projection"]["minimum_reaction_smiles"] == (
        "*c1c(Br)cnn1*.*c1cncc(Br)n1>>*c1cncc(-c2cnn(*)c2*)n1"
    )
    assert result["core_graphic_svg"].count("stroke-dasharray") == 2
