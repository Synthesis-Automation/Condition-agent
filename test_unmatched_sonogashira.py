from chemtools.protocol.recommend import ProtocolRecommender

# Test the unmatched reactions
reactions = [
    ("Diphenylacetylene", "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1"),
    ("tert-butyl (Cl)", "Clc1ccc(C#N)cc1.C#CC(C)(C)C>>CC(C)(C)C#Cc1ccc(C#N)cc1"),
    ("TMS-acetylene (I)", "Ic1ccc(C(F)(F)F)cc1.C#C[Si](C)(C)C>>CC(C)(C)[Si]C#Cc1ccc(C(F)(F)F)cc1"),
    ("Ether-substituted (Cl)", "Clc1ccncc1.C#CCOC>>COC#Cc1ccncc1")
]

recommender = ProtocolRecommender()

for name, smiles in reactions:
    print(f"\n{'='*80}")
    print(f"Reaction: {name}")
    print(f"SMILES: {smiles}")
    print("="*80)
    
    result = recommender.recommend(smiles, k=3, use_smarts_filter=True)
    
    recommendations = result.get('recommended_conditions', [])
    print(f"Matches: {len(recommendations)}")
    
    if recommendations:
        for i, rec in enumerate(recommendations[:3], 1):
            metadata = rec.get('protocol_metadata', {})
            print(f"  {i}. {metadata.get('reaction_family', 'Unknown')} (confidence={rec.get('confidence', 0):.3f})")
    else:
        extras = result.get('extras', {})
        print(f"  No matches - {extras.get('num_after_smarts_filter', 0)} protocols remained after SMARTS filter")
