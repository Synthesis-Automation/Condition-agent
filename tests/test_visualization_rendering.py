from __future__ import annotations

import pytest

from chemtools.util.rdkit_helpers import rdkit_available

pytestmark = pytest.mark.skipif(
    not rdkit_available(), reason="RDKit is required for visualization tests."
)

from chemtools.visualization import rendering  # noqa: E402
from rdkit import Chem  # noqa: E402


def test_promote_hetero_hydrogens_only_targets_hetero_atoms():
    mol = Chem.MolFromSmiles("CCO")
    promoted = rendering._promote_hetero_hydrogens(mol)

    carbons = [atom for atom in promoted.GetAtoms() if atom.GetAtomicNum() == 6]
    for atom in carbons:
        assert atom.GetNumExplicitHs() == 0

    oxygen = next(atom for atom in promoted.GetAtoms() if atom.GetAtomicNum() == 8)
    assert oxygen.GetNumExplicitHs() >= 1


def test_render_molecule_svg(tmp_path):
    output = tmp_path / "phenyl_boronic.svg"
    rendering.render_molecule_image(
        "c1ccc(B(O)O)cc1",
        output,
        image_format="svg",
    )
    data = output.read_text(encoding="utf-8")
    stripped = data.lstrip()
    assert stripped.startswith("<svg") or stripped.startswith("<?xml")
    assert len(data) > 200


def test_render_reaction_png(tmp_path):
    output = tmp_path / "suzuki.png"
    rendering.render_reaction_image(
        "BrC1(c2ccccc2)CC1.c1ccc(B(O)O)cc1>>c1ccc(C2(c3ccccc3)CC2)cc1",
        output,
    )
    assert output.exists()
    assert output.stat().st_size > 2000
