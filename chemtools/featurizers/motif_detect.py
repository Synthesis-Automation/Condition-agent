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
    registry: Optional[Dict[str, Any]] = None,
    discovery_mode: bool = False,
) -> List[Dict[str, Any]]:
    """
    Detect motif hits using compiled compound SMARTS.
    Only returns the most specific matches (highest priority) per bond/site.
    """
    raw_hits: List[Dict[str, Any]] = []
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
            raw_hits.append(
                {
                    "compound_id": compound.compound_id,
                    "match_atom_map": match_atom_map,
                    "a_atom_idx": a_idx,
                    "b_atom_idx": b_idx,
                    "bond": tuple(sorted([a_idx, b_idx])),
                    "priority": compound.priority,
                    "undocumented": False,
                }
            )
    
    if discovery_mode and registry:
        undocumented = discover_undocumented_motifs(mol, registry)
        raw_hits.extend(undocumented)

    # Filter by bond: keep only highest priority per site
    sites: Dict[Tuple[int, ...], List[Dict[str, Any]]] = {}
    for hit in raw_hits:
        site = hit["bond"]
        if site not in sites:
            sites[site] = []
        sites[site].append(hit)
    
    final_hits = []
    for site_hits in sites.values():
        max_priority = max(h["priority"] for h in site_hits)
        # Keep all hits that share the maximum priority for this site
        final_hits.extend([h for h in site_hits if h["priority"] == max_priority])
        
    return final_hits


def discover_undocumented_motifs(
    mol: Any,
    registry: Dict[str, Any],
) -> List[Dict[str, Any]]:
    """
    Find bonds between any scaffold and any substituent that aren't in the registry.
    """
    compiled_groups = registry.get("compiled_groups", {})
    compound_map = registry.get("compound_map", {})
    
    scaffolds = [g for g in compiled_groups.values() if g["kind"] == "scaffold"]
    substituents = [g for g in compiled_groups.values() if g["kind"] == "substituent"]
    
    # 1. Detect all scaffold atoms
    scaffold_atoms: Dict[int, List[Dict[str, Any]]] = {}
    for s in scaffolds:
        matches = mol.GetSubstructMatches(s["query"])
        for match in matches:
            # Find the atom mapped as :1
            a_idx = None
            for atom in s["query"].GetAtoms():
                if atom.GetAtomMapNum() == 1:
                    a_idx = match[atom.GetIdx()]
                    break
            if a_idx is not None:
                scaffold_atoms.setdefault(a_idx, []).append(s)

    # 2. Detect all substituent atoms
    subst_atoms: Dict[int, List[Dict[str, Any]]] = {}
    for sub in substituents:
        matches = mol.GetSubstructMatches(sub["query"])
        for match in matches:
            # Find the atom mapped as :2
            b_idx = None
            for atom in sub["query"].GetAtoms():
                if atom.GetAtomMapNum() == 2:
                    b_idx = match[atom.GetIdx()]
                    break
            if b_idx is not None:
                subst_atoms.setdefault(b_idx, []).append(sub)

    hits = []
    # 3. Check for bonds between scaffold atoms and substituent atoms
    for a_idx, s_list in scaffold_atoms.items():
        mol_atom_a = mol.GetAtomWithIdx(a_idx)
        
        # Direct bonds: Scaffold-Substituent
        for neighbor in mol_atom_a.GetNeighbors():
            b_idx = neighbor.GetIdx()
            if b_idx in subst_atoms:
                for s in s_list:
                    for sub in subst_atoms[b_idx]:
                        _add_undocumented(hits, s, sub, a_idx, b_idx, compound_map)
        
        # Via Oxygen: Scaffold-O-Substituent
        for neighbor in mol_atom_a.GetNeighbors():
            if neighbor.GetSymbol() == "O" and neighbor.GetDegree() == 2:
                o_idx = neighbor.GetIdx()
                for o_neighbor in neighbor.GetNeighbors():
                    b_idx = o_neighbor.GetIdx()
                    if b_idx == a_idx:
                        continue
                    if b_idx in subst_atoms:
                        for s in s_list:
                            for sub in subst_atoms[b_idx]:
                                _add_undocumented(hits, s, sub, a_idx, b_idx, compound_map)

    return hits


def _add_undocumented(hits: List[Dict[str, Any]], s: Dict[str, Any], sub: Dict[str, Any], a_idx: int, b_idx: int, compound_map: Dict[str, Any]) -> None:
    compound_id = f"{s['id']}-{sub['id']}"
    if compound_id in compound_map:
        return
    
    # Avoid duplicates in hits
    for h in hits:
        if h["compound_id"] == compound_id and h["a_atom_idx"] == a_idx and h["b_atom_idx"] == b_idx:
            return

    priority = s["priority"] + sub["priority"]
    hits.append({
        "compound_id": compound_id,
        "match_atom_map": {"1": a_idx, "2": b_idx},
        "a_atom_idx": a_idx,
        "b_atom_idx": b_idx,
        "bond": tuple(sorted([a_idx, b_idx])),
        "priority": priority,
        "undocumented": True,
        "note": "Undocumented combination (Risk: High)",
    })


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
