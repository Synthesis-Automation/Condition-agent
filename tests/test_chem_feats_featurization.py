from __future__ import annotations

from typing import List

import pytest

pytest.importorskip("rdkit")

from chem_feats import featurize_mol, featurize_reaction
from chem_feats.registry import REGISTRY


def _expected_vector_length(roles: List[str]) -> int:
    length = len(REGISTRY.get("global", []))
    for role in roles:
        length += len(REGISTRY.get(role, []))
    fps = REGISTRY.get("fingerprints", {})
    for role in roles:
        cfg = fps.get(role, {})
        length += int(cfg.get("bits", 0))
    return length


def test_featurize_mol_vector_length_and_masks():
    roles = ["amine", "alcohol", "aryl_halide"]
    out = featurize_mol("Nc1ccccc1", roles=roles)

    vec = out.get("vector")
    fields = out.get("fields")
    masks = out.get("masks") or {}

    assert vec is not None
    assert len(vec) == _expected_vector_length(roles)
    assert fields is not None
    assert len(fields) == len(vec)

    assert masks.get("amine") == 1
    assert masks.get("alcohol") == 0
    assert masks.get("aryl_halide") == 0


def test_featurize_mol_detects_alcohol_role():
    roles = ["amine", "alcohol", "aryl_halide"]
    out = featurize_mol("CC(O)C", roles=roles)
    masks = out.get("masks") or {}
    assert masks.get("alcohol") == 1
    assert masks.get("amine") == 0
    assert masks.get("aryl_halide") == 0


def test_featurize_reaction_role_masks():
    rxn = "Brc1ccccc1.Nc1ccccc1>>"
    out = featurize_reaction(rxn)
    reactants = out.get("reactants") or []
    assert len(reactants) == 2

    aryl_entry = next((r for r in reactants if (r.get("masks") or {}).get("aryl_halide") == 1), None)
    amine_entry = next((r for r in reactants if (r.get("masks") or {}).get("amine") == 1), None)

    assert aryl_entry is not None
    assert amine_entry is not None

    for entry in reactants:
        vec = entry.get("vector")
        fields = entry.get("fields")
        assert vec is not None
        assert fields is not None
        assert len(vec) == len(fields) == _expected_vector_length(["amine", "alcohol", "aryl_halide"])
