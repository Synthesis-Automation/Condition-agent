"""Shared deterministic molecular-graph editing primitives."""

from __future__ import annotations

from typing import Any, Callable, Sequence


def bond_type(name: str) -> Any:
    """Return an RDKit bond type for one supported localized bond state."""
    from rdkit import Chem

    values = {
        "SINGLE": Chem.BondType.SINGLE,
        "DOUBLE": Chem.BondType.DOUBLE,
        "TRIPLE": Chem.BondType.TRIPLE,
        "AROMATIC": Chem.BondType.AROMATIC,
    }
    try:
        return values[str(name).upper()]
    except KeyError as exc:
        raise ValueError(f"Unsupported bond order: {name}") from exc


def set_total_hydrogens(
    source_molecule: Any,
    target_molecule: Any,
    atom_index: int,
    delta: int,
) -> bool:
    """Apply an explicit total-H delta without implicit-valence guessing."""
    source_atom = source_molecule.GetAtomWithIdx(atom_index)
    target_count = int(source_atom.GetTotalNumHs(includeNeighbors=True)) + int(delta)
    if target_count < 0:
        return False
    target_atom = target_molecule.GetAtomWithIdx(atom_index)
    target_atom.SetNumExplicitHs(target_count)
    target_atom.SetNoImplicit(True)
    return True


def capture_join_stereochemistry(
    molecule: Any,
    *,
    removed_indices: set[int],
    join_indices: Sequence[int],
) -> tuple[tuple[int, int, int, int, Any], ...]:
    """Capture defined alkene stereo whose leaving substituent is replaced."""
    from rdkit import Chem

    if len(join_indices) != 2:
        return ()
    replacement_pairs = (
        (int(join_indices[0]), int(join_indices[1])),
        (int(join_indices[1]), int(join_indices[0])),
    )
    captured = []
    for bond in molecule.GetBonds():
        stereo = bond.GetStereo()
        stereo_atoms = tuple(int(index) for index in bond.GetStereoAtoms())
        if (
            bond.GetBondType() != Chem.BondType.DOUBLE
            or stereo == Chem.BondStereo.STEREONONE
            or len(stereo_atoms) != 2
            or not any(index in removed_indices for index in stereo_atoms)
        ):
            continue
        endpoints = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        remapped_stereo_atoms = []
        valid = True
        for stereo_atom in stereo_atoms:
            if stereo_atom not in removed_indices:
                remapped_stereo_atoms.append(stereo_atom)
                continue
            replacements = [
                replacement
                for endpoint, replacement in replacement_pairs
                if endpoint in endpoints
                and molecule.GetBondBetweenAtoms(endpoint, stereo_atom) is not None
                and replacement not in removed_indices
            ]
            if len(replacements) != 1:
                valid = False
                break
            remapped_stereo_atoms.append(replacements[0])
        if valid and not any(index in removed_indices for index in endpoints):
            captured.append(
                (
                    int(endpoints[0]),
                    int(endpoints[1]),
                    int(remapped_stereo_atoms[0]),
                    int(remapped_stereo_atoms[1]),
                    stereo,
                )
            )
    return tuple(captured)


def restore_join_stereochemistry(
    molecule: Any,
    *,
    captured: Sequence[tuple[int, int, int, int, Any]],
    shifted: Callable[[int], int],
) -> None:
    """Restore captured alkene stereo after atom deletion and bond formation."""
    from rdkit import Chem

    restored = False
    for endpoint_1, endpoint_2, stereo_atom_1, stereo_atom_2, stereo in captured:
        bond = molecule.GetBondBetweenAtoms(
            shifted(endpoint_1), shifted(endpoint_2)
        )
        if bond is None or bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        bond.SetStereoAtoms(
            shifted(stereo_atom_1),
            shifted(stereo_atom_2),
        )
        bond.SetStereo(stereo)
        restored = True
    if restored:
        Chem.SetDoubleBondNeighborDirections(molecule)
        Chem.AssignStereochemistry(molecule, cleanIt=False, force=True)


__all__ = [
    "bond_type",
    "capture_join_stereochemistry",
    "restore_join_stereochemistry",
    "set_total_hydrogens",
]
