"""Debug primary motif selection for piperazine."""

from chemtools.featurizers.formatters.reaction import featurize_reaction
from chemtools.featurizers.formatters.aggregation import (
    extract_motif_with_bond_info,
    select_primary_motifs_by_atom
)

rxn = "CN1CCNCC1.Nc1cc(Cl)c(Cl)cc1[N+](=O)[O-]>>CN1CCN(c2cc(N)c([N+](=O)[O-])cc2Cl)CC1"
result = featurize_reaction(rxn)

# Manually trace the aggregation logic for reactant 1
r1 = result['reactants'][0]
print("=== Reactant 1 Motifs (Raw) ===")
for m in r1.get('motifs', []):
    print(f"  {m.get('id')}: a_idx={m.get('a_atom_idx')}, fp={m.get('fingerprint')}")

# Extract with bond info
motifs_with_bond = []
for m in r1.get('motifs', []):
    info = extract_motif_with_bond_info(m)
    if info:
        motifs_with_bond.append(info)
        print(f"\nExtracted: {info.get('id')}")
        print(f"  a_idx: {info.get('a_idx')}")
        print(f"  compound_id: {info.get('compound_id')}")
        print(f"  fingerprint: {info.get('fingerprint')}")

print(f"\n=== Before select_primary_motifs_by_atom ===")
print(f"Count: {len(motifs_with_bond)}")
print(f"IDs: {[m.get('id') for m in motifs_with_bond]}")

# Apply per-atom selection
primary = select_primary_motifs_by_atom(motifs_with_bond)

print(f"\n=== After select_primary_motifs_by_atom ===")
print(f"Count: {len(primary)}")
print(f"IDs: {[m.get('id') for m in primary]}")
for m in primary:
    print(f"  {m.get('id')}: a_idx={m.get('a_idx')}, fp={m.get('fingerprint')}")
