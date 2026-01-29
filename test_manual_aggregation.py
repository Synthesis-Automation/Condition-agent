"""Test aggregation with manual tracing."""

from chemtools.featurizers.formatters.aggregation import (
    aggregate_reaction_features,
    normalize_motif_id,
    extract_motif_with_bond_info,
    select_primary_motifs_by_atom,
)

# Simulate reactant bundle
reactant_bundle = {
    "smiles": "CN1CCNCC1",
    "motifs": [
        {"id": "RCH2-NHR", "a_atom_idx": 3, "b_atom_idx": 4, "fingerprint": "N:H1:D2:SP3", "h_count": 1},
        {"id": "RCH2-NR2", "a_atom_idx": 2, "b_atom_idx": 3, "fingerprint": "N:H0:D3:SP3", "h_count": 0},
        {"id": "CH3-NR2", "a_atom_idx": 0, "b_atom_idx": 1, "fingerprint": "N:H0:D3:SP3", "h_count": 0},
    ],
    "context_motifs": [],
}

print("=== Step 1: Extract motifs with bond info ===")
reactant_motifs_full = []
for motif in reactant_bundle["motifs"]:
    if isinstance(motif, dict):
        cid = motif.get("compound_id") or motif.get("id")
        if cid:
            print(f"Processing: {cid}")
            motif_info = extract_motif_with_bond_info(motif)
            if motif_info:
                reactant_motifs_full.append(motif_info)
                print(f"  Added: {motif_info.get('id')}, a_idx={motif_info.get('a_idx')}")
            else:
                print(f"  FAILED to extract")

print(f"\n=== Step 2: Before select_primary_motifs_by_atom ===")
print(f"Count: {len(reactant_motifs_full)}")
print(f"IDs: {[m.get('id') for m in reactant_motifs_full]}")

print(f"\n=== Step 3: Call select_primary_motifs_by_atom ===")
reactant_motifs_full = select_primary_motifs_by_atom(reactant_motifs_full)
print(f"Count: {len(reactant_motifs_full)}")
print(f"IDs: {[m.get('id') for m in reactant_motifs_full]}")

print(f"\n=== Step 4: Extract IDs for primary_motif_ids ===")
primary_motif_ids = [m.get("id", "") for m in reactant_motifs_full if m.get("id")]
print(f"primary_motif_ids: {primary_motif_ids}")
