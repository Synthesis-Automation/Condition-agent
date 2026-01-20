"""
Detect compound motifs in a molecule using compiled SMARTS queries.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from .motif_registry import CompoundPattern


_DISCOVERY_SKIP_SUBSTITUENTS = {
    "R_Subst",
    "Alkyl_Subst",
    "Alkenyl_Subst",
}


def detect_motifs(
    mol: Any,
    compiled_compounds: Iterable[CompoundPattern],
    *,
    max_hits_per_compound: Optional[int] = None,
    registry: Optional[Dict[str, Any]] = None,
    discovery_mode: bool = False,
    site_filter: str = "bond",
) -> List[Dict[str, Any]]:
    """
    Detect motif hits using compiled compound SMARTS.
    Only returns the most specific matches (highest priority) per bond/site.
    
    Args:
        mol: RDKit molecule
        compiled_compounds: List of CompoundPattern
        max_hits_per_compound: Limit hits per compound type
        registry: Full registry (required for discovery_mode)
        discovery_mode: Whether to find undocumented motifs
        site_filter: "bond" (default) keeps best per (scaffold, substituent) bond.
                     "substituent" keeps best per substituent atom (collapses redundant scaffolds).
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
                    "complexity": compound.complexity,
                    "group_b": compound.group_b,
                    "reactivity_weight": compound.reactivity_weight,
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
                    # Identical atom sets.
                    # If they are the same compound type, deduplicate by priority/id.
                    if h1["compound_id"] == h2["compound_id"]:
                        if h1["priority"] < h2["priority"]:
                            subsumed = True
                            break
                        elif h1["priority"] == h2["priority"]:
                            if h1.get("undocumented") and not h2.get("undocumented"):
                                subsumed = True
                                break
                            elif h1.get("undocumented") == h2.get("undocumented"):
                                h1_key = (h1["a_atom_idx"], h1["b_atom_idx"])
                                h2_key = (h2["a_atom_idx"], h2["b_atom_idx"])
                                if h1_key > h2_key:
                                    subsumed = True
                                    h2.setdefault("alt_a_idxs", set()).add(h1["a_atom_idx"])
                                    break
                    else:
                        # Different compound types covering the same atoms (e.g., Ar-NH-R vs R-NH-Ar).
                        # We keep the higher priority one, but record the other's ID.
                        if h1["priority"] < h2["priority"]:
                            h2.setdefault("_alt_hits", []).append(h1)
                            h2.setdefault("alt_a_idxs", set()).add(h1["a_atom_idx"])
                            subsumed = True
                            break
                        elif h1["priority"] == h2["priority"]:
                            # Tie-break by complexity, then ID to be deterministic
                            if (h1["complexity"], h1["compound_id"]) < (h2["complexity"], h2["compound_id"]):
                                h2.setdefault("_alt_hits", []).append(h1)
                                h2.setdefault("alt_a_idxs", set()).add(h1["a_atom_idx"])
                                subsumed = True
                                break
                else:
                    # h1 is a proper subset of h2.
                    # We subsume h1 only if h2 has equal or higher priority.
                    # This prevents perspectival motifs (like H-SR) from eating core motifs (Alkyl-SH).
                    if (h1["priority"], h1["complexity"]) <= (h2["priority"], h2["complexity"]):
                        h2.setdefault("_alt_hits", []).append(h1)
                        subsumed = True
                        break
        if not subsumed:
            filtered_raw.append(h1)
    
    raw_hits = filtered_raw

    # Pass 2: Centric deduplication (for symmetric/multi-attach groups like Amines)
    # If multiple hits have the same compound_id and b_atom_idx, they are often 
    # different attachments/perspectives of the same functional group center.
    centric_filtered = []
    centric_map = {} # (compound_id, b_atom_idx) -> primary_hit
    for h in raw_hits:
        key = (h["compound_id"], h["b_atom_idx"])
        if key not in centric_map:
            centric_map[key] = h
            centric_filtered.append(h)
        else:
            primary = centric_map[key]
            if h["a_atom_idx"] != primary["a_atom_idx"]:
                primary.setdefault("alt_a_idxs", set()).add(h["a_atom_idx"])
            if "alt_a_idxs" in h:
                primary.setdefault("alt_a_idxs", set()).update(h["alt_a_idxs"])
            
            if h["bond"] != primary["bond"]:
                primary.setdefault("alt_bonds", set()).add(h["bond"])
            if h.get("alt_bonds"):
                primary.setdefault("alt_bonds", set()).update(h["alt_bonds"])
            
            # Merge atoms set so the hit covers the whole functional group unit
            primary["atoms"].update(h["atoms"])
            if h.get("_alt_hits"):
                primary.setdefault("_alt_hits", []).extend(h["_alt_hits"])

    raw_hits = centric_filtered

    # Filter by substituent atom and scaffold atom: 
    # keep only highest priority per (functional group, scaffold) pair.
    # This allows a single atom (e.g. Nitrogen) to be part of multiple motifs 
    # if it's attached to different scaffold atoms (e.g. an Ar and an R).
    sites: Dict[Any, List[Dict[str, Any]]] = {}
    for hit in raw_hits:
        if site_filter == "substituent":
            site_key = hit["b_atom_idx"]
        else:
            site_key = (hit["b_atom_idx"], hit["a_atom_idx"])
            
        if site_key not in sites:
            sites[site_key] = []
        sites[site_key].append(hit)
    
    final_hits = []
    for site_hits in sites.values():
        # Sort by specificity (narrowest first)
        site_hits.sort(key=lambda x: (-x["priority"], -x["complexity"], x["compound_id"]))
        
        # In most cases, we only want the absolute best (narrowest) winner.
        # But for symmetric groups, we might have multiple "best" (e.g. Ar-NH-Ar has 2 hits).
        # We group by compound_id to find the unique types.
        max_priority = site_hits[0]["priority"]
        max_complexity = site_hits[0]["complexity"]
        
        # Winners are those sharing the top rank
        winners = []
        others = []
        for h in site_hits:
            if h["priority"] == max_priority and h["complexity"] == max_complexity:
                winners.append(h)
            else:
                others.append(h)
        
        # If multiple winners (rare for one site unless same pattern), deduplicate by ID
        unique_winners = []
        seen_win_ids = set()
        for w in winners:
            if w["compound_id"] not in seen_win_ids:
                unique_winners.append(w)
                seen_win_ids.add(w["compound_id"])
            else:
                # Merge into existing winner
                for prev in unique_winners:
                    if prev["compound_id"] == w["compound_id"]:
                        prev["atoms"].update(w["atoms"])
                        break
        
        # Collect all alternatives for the winner(s)
        # We want a single list of IDs ordered specific -> general
        for best in unique_winners:
            # Gather all alternatives: recorded from Pass 1 + others at this site
            all_alts = best.get("_alt_hits", []) + others
            
            # Deduplicate alternatives and sort
            alt_map = {}
            for a in all_alts:
                cid = a["compound_id"]
                if cid == best["compound_id"]: continue
                if cid not in alt_map:
                    alt_map[cid] = a
                else:
                    if (a["priority"], a["complexity"]) > (alt_map[cid]["priority"], alt_map[cid]["complexity"]):
                        alt_map[cid] = a
            
            sorted_alts = sorted(alt_map.values(), key=lambda x: (-x["priority"], -x["complexity"], x["compound_id"]))
            best["alt_compound_ids"] = [a["compound_id"] for a in sorted_alts]
            final_hits.append(best)
    
    # Final pass: clean up internal tracking fields
    for h in final_hits:
        h.pop("_alt_hits", None)
        if "alt_a_idxs" in h:
            h["alt_a_idxs"] = sorted(list(h["alt_a_idxs"]))
        if "alt_bonds" in h:
            h["alt_bonds"] = sorted(list(h["alt_bonds"]))
        
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
    substituents = [
        g
        for g in compiled_groups.values()
        if g["kind"] == "substituent"
        and g["priority"] > 0
        and g["id"] not in _DISCOVERY_SKIP_SUBSTITUENTS
    ]
    
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
    
    # Avoid duplicates in hits (including symmetric counterparts)
    for h in hits:
        if h["compound_id"] == compound_id and h["atoms"] == atoms:
            # Already have a hit for this set of atoms and compound
            return
    
    # Check if this exact direction is already in hits (legacy check)
    for h in hits:
        if h["compound_id"] == compound_id and h["a_atom_idx"] == a_idx and h["b_atom_idx"] == b_idx:
            return

    priority = s["priority"] + sub["priority"]
    complexity = s.get("complexity", 0) + sub.get("complexity", 0)
    hits.append({
        "compound_id": compound_id,
        "match_atom_map": {"1": a_idx, "2": b_idx},
        "a_atom_idx": a_idx,
        "b_atom_idx": b_idx,
        "bond": tuple(sorted(site_bond)),
        "priority": priority,
        "complexity": complexity,
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
