"""
Fix for atom index collision bug in aggregation.py

The bug: When combining motifs from multiple reactants before calling select_primary_motifs_by_atom(),
different reactants can have the same atom indices (e.g., piperazine atom 3 and aniline atom 3).
This causes select_primary_motifs_by_atom() to incorrectly choose between motifs from DIFFERENT molecules,
losing important nucleophile motifs like RCH2-NHR.

The fix: Apply select_primary_motifs_by_atom() SEPARATELY to each reactant's motifs, then combine the results.
This prevents atom index collisions between different molecules.

Replace lines ~508-555 in aggregation.py with this logic:
"""

# CORRECTED LOGIC (replace the existing loop):

# Phase 2: Collect motifs per reactant (keep separate to avoid atom index collisions)
reactant_motifs_per_mol: List[List[Dict[str, Any]]] = []

for reactant in reactant_list:
    # Steric features
    for entry in reactant.get("steric", {}).get("aryl", []):
        aryl_scores.extend(extract_scores(entry.get("result")))
    for entry in reactant.get("steric", {}).get("alkyl", []):
        alkyl_scores.extend(extract_scores(entry.get("result")))
    for entry in reactant.get("electronics", {}).get("aryl", []):
        electronic_scores.extend(extract_scores(entry.get("result")))
    
    motif_entries = reactant.get("motifs", [])
    context_entries = reactant.get("context_motifs", [])
    for motif in motif_entries:
        compound_id = normalize_motif_id(motif.get("compound_id") or motif.get("id") or "")
        if compound_id:
            motifs.add(str(compound_id))
    
    # Collect full motif dicts with bond info (per reactant)
    motifs_for_this_reactant: List[Dict[str, Any]] = []
    for motif in motif_entries:
        if isinstance(motif, dict):
            cid = motif.get("compound_id") or motif.get("id")
            if cid:
                reactant_motif_ids.append(normalize_motif_id(str(cid)))
                motif_info = extract_motif_with_bond_info(motif)
                if motif_info:
                    motifs_for_this_reactant.append(motif_info)
    if context_entries:
        for motif in context_entries:
            if isinstance(motif, dict):
                cid = motif.get("compound_id") or motif.get("id")
                if cid:
                    reactant_motif_ids.append(normalize_motif_id(str(cid)))
                    motif_info = extract_motif_with_bond_info(motif)
                    if motif_info:
                        motifs_for_this_reactant.append(motif_info)
    
    if motifs_for_this_reactant:
        reactant_motifs_per_mol.append(motifs_for_this_reactant)

# Phase 3: Select primary motif per attachment atom WITHIN EACH REACTANT
# This prevents atom index collisions between different reactants
reactant_motifs_full: List[Dict[str, Any]] = []
for motifs_for_mol in reactant_motifs_per_mol:
    primary_for_mol = select_primary_motifs_by_atom(motifs_for_mol)
    reactant_motifs_full.extend(primary_for_mol)

# Extract IDs from primary motifs for use in change analysis
primary_motif_ids_list = [m.get("id", "") for m in reactant_motifs_full if m.get("id")]
