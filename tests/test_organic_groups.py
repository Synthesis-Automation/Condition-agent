import json
from pathlib import Path

import pytest

from chemtools.taxonomy import loader as taxonomy_loader
from chemtools.core import rdkit as rdkit_helpers
from chemtools.core.smarts import compile_smarts


def _load_groups() -> dict:
    repo_root = Path(__file__).resolve().parents[1]
    groups_path = repo_root / "chemtools" / "taxonomy" / "data" / "organic_groups.v1.3.json"
    payload = json.loads(groups_path.read_text(encoding="utf-8"))
    return {
        group["id"]: group
        for group in payload.get("groups", []) or []
        if isinstance(group, dict) and group.get("id")
    }


def _load_compound_pairs() -> set[tuple[str, str]]:
    payload = taxonomy_loader.load_organic_compounds()
    return {
        (str(entry.get("A") or ""), str(entry.get("B") or ""))
        for entry in payload.get("compounds", []) or []
        if isinstance(entry, dict)
    }


@pytest.mark.skipif(not rdkit_helpers.rdkit_available(), reason="rdkit not available")
def test_isonitrile_group_matches_isocyanide() -> None:
    groups = _load_groups()

    assert "-NC" in groups, "Isonitrile group (-NC) missing from organic groups"
    smarts = groups["-NC"].get("smarts")
    assert smarts, "Isonitrile group missing SMARTS"

    pattern = compile_smarts(smarts, validate=True)
    assert pattern is not None, "Isonitrile SMARTS failed to compile"

    isonitrile = rdkit_helpers.parse_smiles("C[N+]#[C-]")
    nitrile = rdkit_helpers.parse_smiles("C#N")
    assert isonitrile is not None and nitrile is not None, "Failed to parse test SMILES"

    assert isonitrile.HasSubstructMatch(pattern)
    assert not nitrile.HasSubstructMatch(pattern)


@pytest.mark.skipif(not rdkit_helpers.rdkit_available(), reason="rdkit not available")
def test_iodonium_group_matches_diaryliodonium() -> None:
    groups = _load_groups()
    assert "-Iodonium" in groups, "Iodonium group (-Iodonium) missing from organic groups"
    smarts = groups["-Iodonium"].get("smarts")
    assert smarts, "Iodonium group missing SMARTS"

    pattern = compile_smarts(smarts, validate=True)
    assert pattern is not None, "Iodonium SMARTS failed to compile"

    diaryliodonium = rdkit_helpers.parse_smiles("[I+](c1ccccc1)(c1ccccc1)")
    aryl_iodide = rdkit_helpers.parse_smiles("Ic1ccccc1")
    assert diaryliodonium is not None and aryl_iodide is not None, "Failed to parse test SMILES"

    assert diaryliodonium.HasSubstructMatch(pattern)
    assert not aryl_iodide.HasSubstructMatch(pattern)


@pytest.mark.skipif(not rdkit_helpers.rdkit_available(), reason="rdkit not available")
def test_azide_group_matches_both_common_smiles_representations() -> None:
    groups = _load_groups()
    assert "-N3" in groups, "Azide group (-N3) missing from organic groups"
    smarts = groups["-N3"].get("smarts")
    assert smarts, "Azide group missing SMARTS"

    pattern = compile_smarts(smarts, validate=True)
    assert pattern is not None, "Azide SMARTS failed to compile"

    azide_triple_bond_form = rdkit_helpers.parse_smiles("c1ccc([N-][N+]#N)cc1")
    azide_double_bond_form = rdkit_helpers.parse_smiles("c1ccc(N=[N+]=[N-])cc1")
    aniline = rdkit_helpers.parse_smiles("c1ccccc1N")
    assert azide_triple_bond_form is not None and azide_double_bond_form is not None and aniline is not None

    assert azide_triple_bond_form.HasSubstructMatch(pattern)
    assert azide_double_bond_form.HasSubstructMatch(pattern)
    assert not aniline.HasSubstructMatch(pattern)


@pytest.mark.skipif(not rdkit_helpers.rdkit_available(), reason="rdkit not available")
def test_sulfonylhydrazone_group_matches_n_sulfonylhydrazone() -> None:
    groups = _load_groups()
    assert (
        "-Sulfonylhydrazone" in groups
    ), "Sulfonylhydrazone group (-Sulfonylhydrazone) missing from organic groups"
    smarts = groups["-Sulfonylhydrazone"].get("smarts")
    assert smarts, "Sulfonylhydrazone group missing SMARTS"

    pattern = compile_smarts(smarts, validate=True)
    assert pattern is not None, "Sulfonylhydrazone SMARTS failed to compile"

    n_sulfonylhydrazone = rdkit_helpers.parse_smiles("CC(=NNS(=O)(=O)c1ccc(C)cc1)C")
    plain_hydrazone = rdkit_helpers.parse_smiles("CC(=NN)C")
    assert n_sulfonylhydrazone is not None and plain_hydrazone is not None, "Failed to parse test SMILES"

    assert n_sulfonylhydrazone.HasSubstructMatch(pattern)
    assert not plain_hydrazone.HasSubstructMatch(pattern)


@pytest.mark.skipif(not rdkit_helpers.rdkit_available(), reason="rdkit not available")
def test_generic_heteroatom_wildcard_groups_exist_and_compile() -> None:
    groups = _load_groups()
    examples = {
        "-N*": "CN",
        "-O*": "CO",
        "-S*": "CS",
        "-P*": "CP",
    }

    for group_id, smiles in examples.items():
        assert group_id in groups, f"{group_id} missing from organic groups"
        smarts = str(groups[group_id].get("smarts") or "")
        assert smarts, f"{group_id} missing SMARTS"
        pattern = compile_smarts(smarts, validate=True)
        assert pattern is not None, f"{group_id} SMARTS failed to compile"
        mol = rdkit_helpers.parse_smiles(smiles)
        assert mol is not None
        assert mol.HasSubstructMatch(pattern), f"{group_id} did not match {smiles}"


def test_starter_generic_heteroatom_compound_motifs_exist() -> None:
    pairs = _load_compound_pairs()
    expected = {
        ("Ar", "-N*"),
        ("Ar", "-O*"),
        ("Ar", "-P*"),
        ("Ar", "-S*"),
        ("Alkyl", "-N*"),
        ("Alkyl", "-O*"),
        ("Alkyl", "-P*"),
        ("Alkyl", "-S*"),
        ("Alkenyl", "-N*"),
        ("Alkenyl", "-O*"),
        ("Alkenyl", "-P*"),
        ("Alkenyl", "-S*"),
        ("Alkynyl", "-N*"),
        ("Alkynyl", "-O*"),
        ("Alkynyl", "-P*"),
        ("Alkynyl", "-S*"),
    }
    missing = sorted(expected - pairs)
    assert not missing, f"Missing starter wildcard compound motifs: {missing}"
