"""Local context classification shared by all reactive-site detectors."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Set

from chemtools.core.smarts import compile_smarts


_CONTEXTS_PATH = Path(__file__).with_name("data") / "contexts.v1.json"


@lru_cache(maxsize=1)
def load_context_taxonomy() -> Dict[str, Any]:
    """Load context definitions owned by the reactive taxonomy."""
    with _CONTEXTS_PATH.open("r", encoding="utf-8") as handle:
        return json.load(handle)


CONTEXT_PRECEDENCE = tuple(load_context_taxonomy()["precedence"])


@lru_cache(maxsize=1)
def _context_definitions() -> Dict[str, Dict[str, Any]]:
    return {item["id"]: item for item in load_context_taxonomy()["contexts"]}


def _bond_is_double(bond: Any) -> bool:
    return str(bond.GetBondType()) == "DOUBLE"


def _aromatic_ring_token(mol: Any, atom: Any) -> str:
    rings = [ring for ring in mol.GetRingInfo().AtomRings() if atom.GetIdx() in ring]
    if not rings:
        return "Ar"
    ring_atoms: Set[int] = set().union(*(set(ring) for ring in rings))
    return "HeteroAr" if any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in ring_atoms) else "Ar"


def _query_map_position(query: Any, map_number: int) -> int | None:
    for atom in query.GetAtoms():
        if int(atom.GetAtomMapNum()) == int(map_number):
            return int(atom.GetIdx())
    return None


def _matches_mapped_context(
    mol: Any,
    atom_index: int,
    definition: Dict[str, Any],
    excluded: Set[int],
) -> bool:
    query = compile_smarts(str(definition.get("smarts") or ""), validate=False)
    if query is None:
        return False
    anchor_map = (definition.get("atom_roles") or {}).get("context_anchor")
    if anchor_map is None or isinstance(anchor_map, list):
        return False
    position = _query_map_position(query, int(anchor_map))
    if position is None:
        return False
    role_maps = definition.get("atom_roles") or {}
    substituent_maps = role_maps.get("substituent")
    if substituent_maps is None:
        substituent_maps = []
    elif not isinstance(substituent_maps, list):
        substituent_maps = [substituent_maps]
    substituent_positions = [
        mapped_position
        for map_number in substituent_maps
        if (mapped_position := _query_map_position(query, int(map_number))) is not None
    ]
    for match in mol.GetSubstructMatches(query, uniquify=True):
        if int(match[position]) != int(atom_index):
            continue
        if any(int(match[sub_position]) in excluded for sub_position in substituent_positions):
            continue
        return True
    return False


def classify_context(mol: Any, atom_index: int, excluded: Iterable[int] = ()) -> Dict[str, Any]:
    """Classify the retained fragment beginning at ``atom_index``."""
    atom = mol.GetAtomWithIdx(atom_index)
    excluded_set = set(excluded)
    symbol = atom.GetSymbol()
    token = "Other"

    definitions = _context_definitions()
    composite_tokens = [
        candidate for candidate in CONTEXT_PRECEDENCE
        if definitions[candidate].get("classification_method") == "mapped_smarts"
    ]
    for candidate in composite_tokens:
        if _matches_mapped_context(mol, atom_index, definitions[candidate], excluded_set):
            token = candidate
            break
    else:
        if symbol == "C":
            if atom.GetIsAromatic():
                token = _aromatic_ring_token(mol, atom)
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


__all__ = ["CONTEXT_PRECEDENCE", "classify_context", "classify_neighbor_contexts", "load_context_taxonomy"]
