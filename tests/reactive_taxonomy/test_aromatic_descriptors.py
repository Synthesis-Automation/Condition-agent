from rdkit import Chem

from reactive_taxonomy.descriptors.aromatic import build_aromatic_context


def _context(smiles: str, center: int):
    molecule = Chem.MolFromSmiles(smiles)
    assert molecule is not None
    return build_aromatic_context(molecule, center)[0]


def test_aromatic_ring_family_matrix() -> None:
    examples = (
        ("c1ccccc1", 0, "benzene", (6,)),
        ("n1ccccc1", 1, "pyridine", (6,)),
        ("c1ncncc1", 0, "pyrimidine", (6,)),
        ("c1cc[nH]c1", 0, "pyrrole", (5,)),
        ("c1ccoc1", 0, "furan", (5,)),
        ("c1ccsc1", 0, "thiophene", (5,)),
        ("c1ncc[nH]1", 1, "imidazole", (5,)),
        ("c1ccc2ccccc2c1", 0, "naphthalene", (6, 6)),
        ("c1ccc2ncccc2c1", 0, "quinoline_like", (6, 6)),
        ("c1ccc2[nH]ccc2c1", 4, "indole", (5, 6)),
        ("c1ccc2[nH]c3ccccc3c2c1", 4, "carbazole", (5, 6, 6)),
    )
    for smiles, center, family, sizes in examples:
        context = _context(smiles, center)
        assert context.ring_family == family
        assert context.ring_sizes == sizes


def test_aromatic_heteroatoms_use_anchor_relative_positions_and_roles() -> None:
    pyrimidine = _context("c1ncncc1", 0)
    assert {
        (
            item.element,
            item.positional_relation,
            item.aromatic_role,
            item.same_anchor_ring,
        )
        for item in pyrimidine.heteroatoms
    } == {
        ("N", "ortho", "pyridine_like", True),
        ("N", "para", "pyridine_like", True),
    }

    quinoline = _context("c1ccc2ncccc2c1", 0)
    assert [
        item.positional_relation for item in quinoline.heteroatoms
    ] == ["fused_other_ring"]


def test_aromatic_ring_identity_is_smiles_serialization_invariant() -> None:
    first = _context("c1ncccc1", 0)
    second = _context("c1ccccn1", 0)
    assert first.ring_family == second.ring_family == "pyridine"
    assert first.ring_sizes == second.ring_sizes == (6,)
    assert [
        (item.element, item.positional_relation, item.aromatic_role)
        for item in first.heteroatoms
    ] == [
        (item.element, item.positional_relation, item.aromatic_role)
        for item in second.heteroatoms
    ]


def test_ortho_burden_distinguishes_occupancy_from_branch_size() -> None:
    one_methyl = _context("Brc1c(C)cccc1", 1)
    one_tert_butyl = _context("Brc1c(C(C)(C)C)cccc1", 1)

    assert one_methyl.ortho_occupancy_count == 1
    assert one_tert_butyl.ortho_occupancy_count == 1
    assert one_tert_butyl.ortho_burden_score > one_methyl.ortho_burden_score


def test_ring_fusion_is_not_misreported_as_an_ortho_substituent() -> None:
    naphthalene = _context("c1ccc2ccccc2c1", 0)
    indole_n = _context("c1ccc2[nH]ccc2c1", 4)
    assert naphthalene.fused and indole_n.fused
    assert naphthalene.ortho_occupancy_count == 0
    assert indole_n.ortho_occupancy_count == 0
