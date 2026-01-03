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
            
            # Determine the "site bond": the bond from scaffold atom to its immediate neighbor in the motif
            site_bond = _get_site_bond(mol, query, match, a_idx)
            
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
                    "bond": site_bond,
                    "priority": compound.priority,
                    "undocumented": False,
                    "atoms": set(match),
                }
            )
    
    if discovery_mode and registry:
        undocumented = discover_undocumented_motifs(mol, registry)
        # Add undocumented motifs to raw_hits. 
        # We'll filter by priority later.
        raw_hits.extend(undocumented)

    # Filter by subsumption: if hit A's atoms are a subset of hit B's atoms, 
    # then B is more specific. If they match the same atoms, use priority.
    filtered_raw = []
    for i, h1 in enumerate(raw_hits):
        subsumed = False
        for j, h2 in enumerate(raw_hits):
            if i == j:
                continue
            
            h1_atoms = h1["atoms"]
            h2_atoms = h2["atoms"]
            
            if h1_atoms.issubset(h2_atoms):
                if h1_atoms == h2_atoms:
                    # Identical atom sets: use priority as tie-breaker
                    if h1["priority"] < h2["priority"]:
                        subsumed = True
                        break
                    elif h1["priority"] == h2["priority"]:
                        # If priorities are equal, prefer documented over undocumented
                        if h1.get("undocumented") and not h2.get("undocumented"):
                            subsumed = True
                            break
                else:
                    # h1 is a proper subset of h2. h2 is more specific.
                    subsumed = True
                    break
        if not subsumed:
            filtered_raw.append(h1)
    
    raw_hits = filtered_raw

    # Filter by substituent atom: keep only highest priority per functional group attachment
    sites: Dict[int, List[Dict[str, Any]]] = {}
    for hit in raw_hits:
        site = hit["b_atom_idx"]
        if site not in sites:
            sites[site] = []
        sites[site].append(hit)
    
    final_hits = []
    for site_hits in sites.values():
        max_priority = max(h["priority"] for h in site_hits)
        # Keep all hits that share the maximum priority for this site
        # If multiple hits have same priority, prefer documented ones
        best_hits = [h for h in site_hits if h["priority"] == max_priority]
        documented_best = [h for h in best_hits if not h.get("undocumented")]
        if documented_best:
            final_hits.extend(documented_best)
        else:
            final_hits.extend(best_hits)
        
    return final_hits


def discover_undocumented_motifs(
    mol: Any,
    registry: Dict[str, Any],
    skip_h: bool = True,
) -> List[Dict[str, Any]]:
    """
    Find bonds between any scaffold and any substituent that aren't in the registry.
    """
    compiled_groups = registry.get("compiled_groups", {})
    compound_map = registry.get("compound_map", {})
    combination_map = registry.get("combination_map", {})
    
    scaffolds = [g for g in compiled_groups.values() if g["kind"] == "scaffold" and g["priority"] > 0]
    substituents = [g for g in compiled_groups.values() if g["kind"] == "substituent" and g["priority"] > 0]
    
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
    subst_atoms: Dict[int, List[Tuple[Dict[str, Any], Set[int]]]] = {}
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
                # Skip Hydrogen if requested
                if skip_h and mol.GetAtomWithIdx(b_idx).GetSymbol() == "H":
                    continue
                subst_atoms.setdefault(b_idx, []).append((sub, set(match)))

    hits = []
    # 3. Check for bonds between scaffold atoms and substituent atoms
    for a_idx, s_list in scaffold_atoms.items():
        mol_atom_a = mol.GetAtomWithIdx(a_idx)
        
        # Direct bonds: Scaffold-Substituent
        for neighbor in mol_atom_a.GetNeighbors():
            b_idx = neighbor.GetIdx()
            if b_idx in subst_atoms:
                # Skip if the bond is in a ring (usually part of the scaffold itself)
                bond_obj = mol.GetBondBetweenAtoms(a_idx, b_idx)
                if bond_obj and bond_obj.IsInRing():
                    continue
                    
                for s in s_list:
                    for sub, sub_match_atoms in subst_atoms[b_idx]:
                        # Prevent self-intersection: scaffold atom cannot be part of substituent
                        if a_idx in sub_match_atoms:
                            continue
                        _add_undocumented(
                            hits, s, sub, a_idx, b_idx, compound_map, combination_map, 
                            site_bond=(a_idx, b_idx),
                            atoms={a_idx} | sub_match_atoms
                        )
        
        # Via Oxygen: Scaffold-O-Substituent
        for neighbor in mol_atom_a.GetNeighbors():
            if neighbor.GetSymbol() == "O" and neighbor.GetDegree() == 2:
                o_idx = neighbor.GetIdx()
                # Skip if the Scaffold-O bond is in a ring
                bond_ao = mol.GetBondBetweenAtoms(a_idx, o_idx)
                if bond_ao and bond_ao.IsInRing():
                    continue
                    
                for o_neighbor in neighbor.GetNeighbors():
                    b_idx = o_neighbor.GetIdx()
                    if b_idx == a_idx:
                        continue
                    if b_idx in subst_atoms:
                        for s in s_list:
                            for sub, sub_match_atoms in subst_atoms[b_idx]:
                                # Prevent self-intersection
                                if a_idx in sub_match_atoms:
                                    continue
                                _add_undocumented(
                                    hits, s, sub, a_idx, b_idx, compound_map, combination_map, 
                                    site_bond=(a_idx, o_idx),
                                    atoms={a_idx, o_idx} | sub_match_atoms
                                )

    return hits


def _add_undocumented(
    hits: List[Dict[str, Any]], 
    s: Dict[str, Any], 
    sub: Dict[str, Any], 
    a_idx: int, 
    b_idx: int, 
    compound_map: Dict[str, Any],
    combination_map: Dict[Tuple[str, str], str],
    site_bond: Tuple[int, int],
    atoms: Set[int]
) -> None:
    # Check if this combination is already documented
    combo_id = combination_map.get((s["id"], sub["id"]))
    if combo_id and combo_id in compound_map:
        return

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
        "bond": tuple(sorted(site_bond)),
        "priority": priority,
        "undocumented": True,
        "atoms": atoms,
        "note": "Undocumented combination (Risk: High)",
    })


def _get_site_bond(mol: Any, query: Any, match: Tuple[int, ...], a_idx: int) -> Tuple[int, int]:
    """
    Find the bond from the scaffold atom (a_idx) to its immediate neighbor in the motif.
    """
    q_a_idx = _find_query_map_atom(query, map_num=1)
    if q_a_idx is not None:
        q_a_atom = query.GetAtomWithIdx(q_a_idx)
        for q_neighbor in q_a_atom.GetNeighbors():
            # Find the corresponding atom in the molecule match
            m_neighbor_idx = match[q_neighbor.GetIdx()]
            # Return the bond as a sorted tuple
            return tuple(sorted([a_idx, m_neighbor_idx]))
    
    # Fallback to a_idx itself if no neighbor found (e.g. single atom motif)
    return (a_idx, a_idx)


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
