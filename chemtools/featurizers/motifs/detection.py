"""Motif detection using compiled SMARTS patterns."""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from .models import CompoundPattern, _DISCOVERY_SKIP_SUBSTITUENTS


def compute_motif_fingerprint(mol: Any, matched_atoms: Set[int]) -> str:
    """
    Compute a fingerprint string that captures the state of matched atoms.
    
    This fingerprint includes:
    - For each matched heteroatom: element, H count, heavy-atom degree, hybridization
    - Sorted by atom index for consistency
    
    The fingerprint allows detecting when a functional group has changed
    (e.g., hydrazide NH2 -> hydrazone N=) even if the core SMARTS still matches.
    
    Works correctly with molecules that have explicit H's (via AddHs).
    
    Args:
        mol: RDKit molecule object
        matched_atoms: Set of atom indices that matched the SMARTS
        
    Returns:
        A fingerprint string like "N:H2:D1:SP3,N:H1:D2:SP2" (sorted by index)
    """
    if not matched_atoms:
        return ""
    
    try:
        parts = []
        for idx in sorted(matched_atoms):
            atom = mol.GetAtomWithIdx(idx)
            symbol = atom.GetSymbol()
            # Only fingerprint heteroatoms (N, O, S, P) - these are where H-count matters
            if symbol not in ("N", "O", "S", "P"):
                continue
            
            # Count H neighbors (works for both explicit and implicit H molecules)
            # GetTotalNumHs() returns 0 after AddHs(), so count explicit H neighbors
            h_count = 0
            heavy_neighbors = 0
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == "H":
                    h_count += 1
                else:
                    heavy_neighbors += 1
            
            # Also add implicit H if any (shouldn't be any after AddHs, but be safe)
            h_count += atom.GetNumImplicitHs()
            
            # Get hybridization
            hyb = str(atom.GetHybridization()).split(".")[-1]  # SP3, SP2, SP, etc
            
            parts.append(f"{symbol}:H{h_count}:D{heavy_neighbors}:{hyb}")
        return ",".join(sorted(parts))
    except Exception:
        return ""


def compute_group_fingerprint(
    mol: Any,
    matched_atoms: Set[int],
    *,
    a_idx: Optional[int],
    b_idx: Optional[int],
) -> str:
    """
    Compute a group-level fingerprint focused on the substituent side.

    Uses the matched atom subgraph connected to b_idx, excluding the scaffold
    attachment atom (a_idx) when possible. Falls back to full motif fingerprint.
    """
    if b_idx is None or a_idx is None or not matched_atoms:
        return compute_motif_fingerprint(mol, matched_atoms)
    try:
        visited: Set[int] = set()
        stack = [b_idx]
        exclude = {a_idx} if a_idx != b_idx else set()
        while stack:
            idx = stack.pop()
            if idx in visited or idx in exclude:
                continue
            visited.add(idx)
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                if n_idx in matched_atoms and n_idx not in visited:
                    stack.append(n_idx)
        if not visited:
            return compute_motif_fingerprint(mol, matched_atoms)
        return compute_motif_fingerprint(mol, visited)
    except Exception:
        return compute_motif_fingerprint(mol, matched_atoms)


def compute_extended_fingerprint(mol: Any, matched_atoms: Set[int], a_idx: int, b_idx: int) -> Dict[str, Any]:
    """
    Compute extended fingerprint with structured data for fine-grained comparison.
    
    Returns a dict with:
    - heteroatom_h_counts: Dict mapping atom index to H count for N/O/S/P
    - total_h_on_heteroatoms: Sum of H on all heteroatoms (quick comparison)
    - fingerprint_str: String representation for hashing/comparison
    
    Works correctly with molecules that have explicit H's (via AddHs).
    
    Args:
        mol: RDKit molecule
        matched_atoms: Set of matched atom indices
        a_idx: Scaffold attachment atom index
        b_idx: Substituent attachment atom index
    """
    result = {
        "heteroatom_h_counts": {},
        "total_h_on_heteroatoms": 0,
        "fingerprint_str": "",
    }
    
    if not matched_atoms:
        return result
    
    try:
        parts = []
        for idx in sorted(matched_atoms):
            atom = mol.GetAtomWithIdx(idx)
            symbol = atom.GetSymbol()
            if symbol in ("N", "O", "S", "P"):
                # Count H neighbors (works for molecules with explicit H's)
                h_count = 0
                heavy_neighbors = 0
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == "H":
                        h_count += 1
                    else:
                        heavy_neighbors += 1
                # Also add implicit H if any
                h_count += atom.GetNumImplicitHs()
                
                result["heteroatom_h_counts"][idx] = h_count
                result["total_h_on_heteroatoms"] += h_count
                parts.append(f"{symbol}:H{h_count}:D{heavy_neighbors}")
        
        result["fingerprint_str"] = ",".join(parts)
    except Exception:
        pass
    
    return result


def detect_motifs(
    mol: Any,
    compiled_compounds: Iterable[CompoundPattern],
    *,
    max_hits_per_compound: Optional[int] = None,
    registry: Optional[Dict[str, Any]] = None,
    discovery_mode: bool = False,
    site_filter: str = "bond",
) -> List[Dict[str, Any]]:
    """Detect compound motifs in a molecule using SMARTS patterns.
    
    Returns only the most specific matches (highest priority) per bond/site.
    Performs multi-stage filtering:
    1. Subsumption: Remove hits that are subsets of more specific hits
    2. Centric deduplication: Merge symmetric attachments to same functional group
    3. Site filtering: Keep best match per site (by priority/complexity)
    
    Args:
        mol: RDKit molecule
        compiled_compounds: List of CompoundPattern objects
        max_hits_per_compound: Limit hits per compound type
        registry: Full registry (required for discovery_mode)
        discovery_mode: Whether to find undocumented motif combinations
        site_filter: "bond" (default) keeps best per (scaffold, substituent) bond.
                     "substituent" keeps best per substituent atom (collapses scaffolds).
                     
    Returns:
        List of motif hit dictionaries with keys:
        - compound_id: Motif identifier (e.g., "Ar-Cl")
        - a_atom_idx: Scaffold attachment atom index
        - b_atom_idx: Substituent functional atom index
        - bond: Site bond tuple (sorted atom indices)
        - priority: Detection priority score
        - complexity: SMARTS complexity score
        - atoms: Set of all matched atom indices
        - alt_compound_ids: Alternative (less specific) compound IDs for same site
        - alt_a_idxs: Alternative scaffold attachment points (for symmetric groups)
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

            # Determine the "site bond": the bond from scaffold atom to its immediate neighbor
            site_bond = _get_site_bond(mol, query, match, a_idx)

            key = (compound.compound_id, a_idx, b_idx)
            if key in seen:
                continue
            seen.add(key)
            
            # Compute fingerprint capturing H-count and connectivity on heteroatoms
            matched_atom_set = set(match)
            fingerprint = compute_group_fingerprint(
                mol,
                matched_atom_set,
                a_idx=a_idx,
                b_idx=b_idx,
            )
            extended_fp = compute_extended_fingerprint(mol, matched_atom_set, a_idx, b_idx)
            
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
                    "atoms": matched_atom_set,
                    "fingerprint": fingerprint,
                    "h_count": extended_fp["total_h_on_heteroatoms"],
                }
            )

    if discovery_mode and registry:
        undocumented = discover_undocumented_motifs(mol, registry)
        raw_hits.extend(undocumented)

    # Pass 1: Filter by subsumption
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
                    # Identical atom sets
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
                                    h2.setdefault("alt_a_idxs", set()).add(
                                        h1["a_atom_idx"]
                                    )
                                    break
                    else:
                        # Different compound types covering same atoms
                        if h1["priority"] < h2["priority"]:
                            h2.setdefault("_alt_hits", []).append(h1)
                            h2.setdefault("alt_a_idxs", set()).add(h1["a_atom_idx"])
                            subsumed = True
                            break
                        elif h1["priority"] == h2["priority"]:
                            # Tie-break by complexity, then ID
                            if (h1["complexity"], h1["compound_id"]) < (
                                h2["complexity"],
                                h2["compound_id"],
                            ):
                                h2.setdefault("_alt_hits", []).append(h1)
                                h2.setdefault("alt_a_idxs", set()).add(h1["a_atom_idx"])
                                subsumed = True
                                break
                else:
                    # h1 is a proper subset of h2
                    if (h1["priority"], h1["complexity"]) <= (
                        h2["priority"],
                        h2["complexity"],
                    ):
                        h2.setdefault("_alt_hits", []).append(h1)
                        subsumed = True
                        break
        if not subsumed:
            filtered_raw.append(h1)

    raw_hits = filtered_raw

    # Pass 2: Centric deduplication (symmetric/multi-attach groups)
    centric_filtered = []
    centric_map = {}  # (compound_id, b_atom_idx) -> primary_hit
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

            # Merge atoms set
            primary["atoms"].update(h["atoms"])
            if h.get("_alt_hits"):
                primary.setdefault("_alt_hits", []).extend(h["_alt_hits"])

    raw_hits = centric_filtered

    # Pass 3: Filter by site - keep only highest priority per site
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
        site_hits.sort(
            key=lambda x: (-x["priority"], -x["complexity"], x["compound_id"])
        )

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

        # Deduplicate winners by compound ID
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

        # Collect all alternatives
        for best in unique_winners:
            all_alts = best.get("_alt_hits", []) + others

            # Deduplicate alternatives and sort
            alt_map = {}
            for a in all_alts:
                cid = a["compound_id"]
                if cid == best["compound_id"]:
                    continue
                if cid not in alt_map:
                    alt_map[cid] = a
                else:
                    if (a["priority"], a["complexity"]) > (
                        alt_map[cid]["priority"],
                        alt_map[cid]["complexity"],
                    ):
                        alt_map[cid] = a

            sorted_alts = sorted(
                alt_map.values(),
                key=lambda x: (-x["priority"], -x["complexity"], x["compound_id"]),
            )
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
    """Find bonds between scaffolds and substituents not in the registry.
    
    This "discovery mode" detects potential motifs that aren't explicitly
    defined in the compounds file by combining any scaffold with any
    substituent via direct bonds or via oxygen linkers.
    
    Args:
        mol: RDKit molecule
        registry: Full compiled registry
        skip_h: Whether to skip hydrogen substituents
        
    Returns:
        List of undocumented motif hits (marked with undocumented=True)
    """
    compiled_groups = registry.get("compiled_groups", {})
    compound_map = registry.get("compound_map", {})
    combination_map = registry.get("combination_map", {})

    scaffolds = [
        g
        for g in compiled_groups.values()
        if g["kind"] == "scaffold" and g["priority"] > 0
    ]
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
                # Skip if the bond is in a ring
                bond_obj = mol.GetBondBetweenAtoms(a_idx, b_idx)
                if bond_obj and bond_obj.IsInRing():
                    continue

                for s in s_list:
                    for sub, sub_match_atoms in subst_atoms[b_idx]:
                        # Prevent self-intersection
                        if a_idx in sub_match_atoms:
                            continue
                        _add_undocumented(
                            hits,
                            s,
                            sub,
                            a_idx,
                            b_idx,
                            compound_map,
                            combination_map,
                            site_bond=(a_idx, b_idx),
                            atoms={a_idx} | sub_match_atoms,
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
                                    hits,
                                    s,
                                    sub,
                                    a_idx,
                                    b_idx,
                                    compound_map,
                                    combination_map,
                                    site_bond=(a_idx, o_idx),
                                    atoms={a_idx, o_idx} | sub_match_atoms,
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
    atoms: Set[int],
) -> None:
    """Add an undocumented motif hit if not already documented."""
    # Check if this combination is already documented
    combo_id = combination_map.get((s["id"], sub["id"]))
    if combo_id and combo_id in compound_map:
        return

    compound_id = f"{s['id']}-{sub['id']}"
    if compound_id in compound_map:
        return

    # Avoid duplicates in hits
    for h in hits:
        if h["compound_id"] == compound_id and h["atoms"] == atoms:
            return

    for h in hits:
        if (
            h["compound_id"] == compound_id
            and h["a_atom_idx"] == a_idx
            and h["b_atom_idx"] == b_idx
        ):
            return

    priority = s["priority"] + sub["priority"]
    complexity = s.get("complexity", 0) + sub.get("complexity", 0)
    hits.append(
        {
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
        }
    )


def _get_site_bond(
    mol: Any, query: Any, match: Tuple[int, ...], a_idx: int
) -> Tuple[int, int]:
    """Find the bond from scaffold atom to its immediate neighbor in the motif."""
    q_a_idx = _find_query_map_atom(query, map_num=1)
    if q_a_idx is not None:
        q_a_atom = query.GetAtomWithIdx(q_a_idx)
        for q_neighbor in q_a_atom.GetNeighbors():
            # Find the corresponding atom in the molecule match
            m_neighbor_idx = match[q_neighbor.GetIdx()]
            # Return the bond as a sorted tuple
            return tuple(sorted([a_idx, m_neighbor_idx]))

    # Fallback to a_idx itself if no neighbor found
    return (a_idx, a_idx)


def _extract_mapped_atoms(
    query: Any, match: Tuple[int, ...]
) -> Tuple[Dict[str, int], Optional[int], Optional[int]]:
    """Extract atom map numbers from query and match."""
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


def _infer_b_idx(
    mol: Any, query: Any, match: Tuple[int, ...], a_idx: int
) -> Optional[int]:
    """Infer substituent atom index when not explicitly mapped."""
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

    non_aromatic = [
        idx for idx in candidates if not mol.GetAtomWithIdx(idx).GetIsAromatic()
    ]
    if non_aromatic:
        return min(non_aromatic)
    return min(candidates)


def _find_query_map_atom(query: Any, *, map_num: int) -> Optional[int]:
    """Find the query atom with the specified atom map number."""
    for atom in query.GetAtoms():
        if atom.GetAtomMapNum() == map_num:
            return atom.GetIdx()
    return None


__all__ = [
    "detect_motifs",
    "discover_undocumented_motifs",
]
