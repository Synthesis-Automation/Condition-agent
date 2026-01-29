"""Direct test of aggregation logic."""

from chemtools.featurizers.formatters.reaction import featurize_reaction
from chemtools.featurizers.formatters.aggregation import normalize_motif_id, extract_motif_with_bond_info

rxn = "CN1CCNCC1.Nc1cc(Cl)c(Cl)cc1[N+](=O)[O-]>>CN1CCN(c2cc(N)c([N+](=O)[O-])cc2Cl)CC1"
result = featurize_reaction(rxn)

# Manually trace first reactant
r1 = result['reactants'][0]
print("=== Reactant 1 Raw Motifs ===")
motif_entries = r1.get("motifs", [])
for motif in motif_entries:
    cid_field = motif.get("compound_id")
    id_field = motif.get("id")
    cid = cid_field or id_field
    print(f"  compound_id={cid_field}, id={id_field}, final_cid={cid}")
    if cid:
        norm = normalize_motif_id(str(cid))
        print(f"    Normalized: {norm}")
        motif_info = extract_motif_with_bond_info(motif)
        if motif_info:
            print(f"    Extracted: {motif_info.get('id')}")
        else:
            print(f"    Extraction FAILED")
