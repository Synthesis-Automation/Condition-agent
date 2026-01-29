"""Debug Suzuki reaction motif detection."""
from chemtools.featurizers.unified import featurize_reaction

rxn = 'Brc1cncnc1.OB(O)c1ccoc1>>C1(C2=COC=C2)=CN=CN=C1'
result = featurize_reaction(rxn, options={"detailed": True})

print('=== Reactant Motifs ===')
for i, r in enumerate(result.get('reactants', [])):
    smiles = r.get('smiles') or r.get('input') or f'reactant_{i}'
    print(f'Reactant {i+1} ({smiles}):')
    motifs = r.get('motifs', [])
    print(f'  Raw motifs type: {type(motifs)}, len: {len(motifs)}')
    for m in motifs:
        if isinstance(m, dict):
            # Try different key names
            cid = m.get("compound_id") or m.get("id") or "NO_ID"
            a_idx = m.get("a_atom_idx") or m.get("a_idx")
            rw = m.get("reactivity_weight", 0)
            print(f'  - {cid} (a_idx={a_idx}, rw={rw})')
        else:
            print(f'  - {m} (type: {type(m)})')
    if not motifs:
        print('  (no motifs detected)')

print()
print('=== Reaction Type ===')
rt = result.get('reaction_type', {})
if isinstance(rt, dict):
    print(f"  type: {rt.get('reaction_type')}")
    print(f"  confidence: {rt.get('confidence')}")
else:
    print(f"  type: {rt}")

print()
print('=== Aggregates ===')
agg = result.get('aggregates', {})
print(f'motif_ids: {agg.get("motif_ids")}')
print(f'primary_motif_ids: {agg.get("primary_motif_ids")}')
print(f'reacted_motifs: {agg.get("reacted_motifs")}')
print(f'formed_motifs: {agg.get("formed_motifs")}')
print(f'spectator_motifs: {agg.get("spectator_motifs")}')

print()
print('=== Reaction Key ===')
print(result.get('reaction_key'))
