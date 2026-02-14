import json
from pathlib import Path

import pytest

from chemtools.util import rdkit_helpers
from chemtools.util.smarts_cache import compile_smarts


def _load_groups() -> dict:
    repo_root = Path(__file__).resolve().parents[1]
    groups_path = repo_root / "chemtools" / "taxonomy" / "data" / "organic_groups.v1.3.json"
    payload = json.loads(groups_path.read_text(encoding="utf-8"))
    return {
        group["id"]: group
        for group in payload.get("groups", []) or []
        if isinstance(group, dict) and group.get("id")
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
