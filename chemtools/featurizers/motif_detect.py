"""
Detect compound motifs in a molecule using compiled SMARTS queries.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from .motif_registry import CompoundPattern


def detect_motifs(
    mol: Any,
    compiled_compounds: Iterable[CompoundPattern],
    *,
    max_hits_per_compound: Optional[int] = None,
) -> List[Dict[str, Any]]:
    """
    Detect motif hits using compiled compound SMARTS.
    """
    hits: List[Dict[str, Any]] = []
    seen: Set[Tuple[str, int, int]] = set()

    for compound in compiled_compounds:
        query = compound.query
        if query is None:
            continue
        matches = mol.GetSubstructMatches(query, uniquify=True)
        if max_hits_per_compound:
            matches = matches[:max_hits_per_compound]
        for match in matches:
            match_atom_map, a_idx, b_idx = _extract_mapped_atoms(query, match)
            if a_idx is None:
                continue
            if b_idx is None:
                b_idx = _infer_b_idx(mol, query, match, a_idx)
                if b_idx is not None:
                    match_atom_map["2"] = b_idx
            if b_idx is None and _find_query_map_atom(query, map_num=2) is None:
                # Allow single-atom motifs (e.g., reagents) to reuse the A atom as B.
                b_idx = a_idx
                match_atom_map["2"] = a_idx
            if b_idx is None:
                continue
            key = (compound.compound_id, a_idx, b_idx)
            if key in seen:
                continue
            seen.add(key)
            hits.append(
                {
                    "compound_id": compound.compound_id,
                    "match_atom_map": match_atom_map,
                    "a_atom_idx": a_idx,
                    "b_atom_idx": b_idx,
                    "bond": [a_idx, b_idx],
                }
            )
    return hits


def _extract_mapped_atoms(query: Any, match: Tuple[int, ...]) -> Tuple[Dict[str, int], Optional[int], Optional[int]]:
    match_atom_map: Dict[str, int] = {}
    a_idx: Optional[int] = None
    b_idx: Optional[int] = None
    for atom in query.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if not map_num:
            continue
        mol_idx = match[atom.GetIdx()]
        match_atom_map[str(map_num)] = mol_idx
        if map_num == 1:
            a_idx = mol_idx
        elif map_num == 2:
            b_idx = mol_idx
    return match_atom_map, a_idx, b_idx


def _infer_b_idx(mol: Any, query: Any, match: Tuple[int, ...], a_idx: int) -> Optional[int]:
    q_a_idx = _find_query_map_atom(query, map_num=1)
    if q_a_idx is not None:
        q_a_atom = query.GetAtomWithIdx(q_a_idx)
        q_neighbors = list(q_a_atom.GetNeighbors())
        if len(q_neighbors) == 1:
            return match[q_neighbors[0].GetIdx()]

    matched_atoms = set(match)
    candidates = [
        nbr.GetIdx()
        for nbr in mol.GetAtomWithIdx(a_idx).GetNeighbors()
        if nbr.GetIdx() in matched_atoms
    ]
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]

    non_aromatic = [idx for idx in candidates if not mol.GetAtomWithIdx(idx).GetIsAromatic()]
    if non_aromatic:
        return min(non_aromatic)
    return min(candidates)


def _find_query_map_atom(query: Any, *, map_num: int) -> Optional[int]:
    for atom in query.GetAtoms():
        if atom.GetAtomMapNum() == map_num:
            return atom.GetIdx()
    return None
