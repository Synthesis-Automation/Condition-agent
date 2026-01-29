"""Debug piperazine SNAr detection issue."""

from chemtools.featurizers.formatters.reaction import featurize_reaction

rxn = "CN1CCNCC1.Nc1cc(Cl)c(Cl)cc1[N+](=O)[O-]>>CN1CCN(c2cc(N)c([N+](=O)[O-])cc2Cl)CC1"
result = featurize_reaction(rxn)

print("=" * 60)
print("REACTION:", rxn)
print("=" * 60)
print(f"Detected as: {result['reaction_type']}")
print(f"Confidence: {result.get('confidence')}")
print()

# Check reactant 1 (piperazine)
r1 = result['reactants'][0]
print("=== Reactant 1 (Piperazine) ===")
print(f"SMILES: {r1.get('smiles')}")
print("All motifs:")
for m in r1.get('motifs', []):
    print(f"  {m.get('id')}: a_idx={m.get('a_atom_idx')}, weight={m.get('reactivity_weight')}, fp={m.get('fingerprint')}")

print()
print("=== Aggregation Results ===")
print(f"Primary motif IDs: {result['aggregates'].get('primary_motif_ids')}")
print(f"Reacted: {result['aggregates']['reacted_motifs']}")
print(f"Formed: {result['aggregates']['formed_motifs']}")
print(f"Spectators: {result['aggregates']['spectator_motifs']}")
