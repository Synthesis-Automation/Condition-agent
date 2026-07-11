"""Local context classification shared by all reactive-site detectors."""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Set


CONTEXT_PRECEDENCE = (
    "SO2R", "P(O)R2", "C(O)NR", "C(O)OR", "C(O)R",
    "HeteroAr", "Ar", "Alkenyl", "Alkynyl", "Alkyl",
    "N", "O", "S", "Other",
)


def _bond_is_double(bond: Any) -> bool:
    return str(bond.GetBondType()) == "DOUBLE"


def _aromatic_ring_token(mol: Any, atom: Any) -> str:
    rings = [ring for ring in mol.GetRingInfo().AtomRings() if atom.GetIdx() in ring]
    if not rings:
        return "Ar"
    ring_atoms: Set[int] = set().union(*(set(ring) for ring in rings))
    return "HeteroAr" if any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in ring_atoms) else "Ar"


def classify_context(mol: Any, atom_index: int, excluded: Iterable[int] = ()) -> Dict[str, Any]:
    """Classify the retained fragment beginning at ``atom_index``."""
    atom = mol.GetAtomWithIdx(atom_index)
    excluded_set = set(excluded)
    symbol = atom.GetSymbol()
    token = "Other"

    if symbol == "S":
        oxo = sum(1 for bond in atom.GetBonds() if _bond_is_double(bond) and bond.GetOtherAtom(atom).GetSymbol() == "O")
        token = "SO2R" if oxo >= 2 else "S"
    elif symbol == "P":
        has_oxo = any(_bond_is_double(bond) and bond.GetOtherAtom(atom).GetSymbol() == "O" for bond in atom.GetBonds())
        token = "P(O)R2" if has_oxo else "Other"
    elif symbol == "C":
        if atom.GetIsAromatic():
            token = _aromatic_ring_token(mol, atom)
        else:
            double_o = any(_bond_is_double(bond) and bond.GetOtherAtom(atom).GetSymbol() == "O" for bond in atom.GetBonds())
            if double_o:
                single_hetero = [
                    bond.GetOtherAtom(atom).GetSymbol()
                    for bond in atom.GetBonds()
                    if not _bond_is_double(bond) and bond.GetOtherAtom(atom).GetIdx() not in excluded_set
                ]
                if "N" in single_hetero:
                    token = "C(O)NR"
                elif "O" in single_hetero:
                    token = "C(O)OR"
                else:
                    token = "C(O)R"
            elif str(atom.GetHybridization()) == "SP":
                token = "Alkynyl"
            elif str(atom.GetHybridization()) == "SP2":
                token = "Alkenyl"
            else:
                token = "Alkyl"
    elif symbol in {"N", "O", "S"}:
        token = symbol

    return {"token": token, "atom_indices": [atom_index], "features": {}}


def classify_neighbor_contexts(mol: Any, center_index: int, excluded: Iterable[int] = ()) -> List[Dict[str, Any]]:
    """Return retained heavy-atom contexts directly attached to a center."""
    excluded_set = set(excluded) | {center_index}
    center = mol.GetAtomWithIdx(center_index)
    contexts = [
        classify_context(mol, neighbor.GetIdx(), excluded_set)
        for neighbor in center.GetNeighbors()
        if neighbor.GetAtomicNum() > 1 and neighbor.GetIdx() not in excluded_set
    ]
    rank = {token: i for i, token in enumerate(CONTEXT_PRECEDENCE)}
    return sorted(contexts, key=lambda item: (rank.get(item["token"], 999), item["token"], item["atom_indices"]))


__all__ = ["CONTEXT_PRECEDENCE", "classify_context", "classify_neighbor_contexts"]
