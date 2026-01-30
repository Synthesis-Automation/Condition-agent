import json
from pathlib import Path

import pytest

from chemtools.util import rdkit_helpers
from chemtools.util.smarts_cache import compile_smarts


@pytest.mark.skipif(not rdkit_helpers.rdkit_available(), reason="rdkit not available")
def test_isonitrile_group_matches_isocyanide() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    groups_path = repo_root / "chemtools" / "taxonomy" / "data" / "organic_groups.v1.3.json"

    payload = json.loads(groups_path.read_text(encoding="utf-8"))
    groups = {
        group["id"]: group
        for group in payload.get("groups", []) or []
        if isinstance(group, dict) and group.get("id")
    }

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
