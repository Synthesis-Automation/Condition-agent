from chemtools.protocol.recommend import ProtocolRecommender

# Test Sonogashira reactions
test_reactions = [
    ("Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1", "Sonogashira - Diphenylacetylene"),
    ("Ic1ccncc1.C#CC>>C#Cc1ccncc1", "Sonogashira - Pyridine acetylene"),
    ("Clc1ccc(C#N)cc1.C#CC(C)(C)C>>CC(C)(C)C#Cc1ccc(C#N)cc1", "Sonogashira - tert-butyl"),
]

# Database protocol
protocol_smarts = '[c,n,o,s]I.C#C>>[c,n,o,s]C#C'
protocol_reaction = 'NC(=O)c1ccccc1I.C#CCCCC>>NC(=O)c1ccccc1C#CCCCC'

print("=" * 80)
print("SONOGASHIRA MATCHING ANALYSIS")
print("=" * 80)
print()
print(f"Protocol in database:")
print(f"  Reaction: {protocol_reaction}")
print(f"  SMARTS: {protocol_smarts}")
print()
print("=" * 80)
print()

# Test each reaction
recommender = ProtocolRecommender()

for smiles, desc in test_reactions:
    print(f"Testing: {desc}")
    print(f"  SMILES: {smiles}")
    
    # Test with SMARTS filtering
    result_with_filter = recommender.recommend(smiles, k=3, use_smarts_filter=True)
    num_matches = len(result_with_filter.get('protocols', []))
    
    # Test without SMARTS filtering
    result_no_filter = recommender.recommend(smiles, k=3, use_smarts_filter=False)
    num_no_filter = len(result_no_filter.get('protocols', []))
    
    print(f"  With SMARTS filter: {num_matches} match(es)")
    print(f"  Without SMARTS filter: {num_no_filter} match(es)")
    
    if num_matches == 0 and num_no_filter > 0:
        print(f"  --> SMARTS filtering removed all matches!")
        print(f"  --> Issue: Protocol SMARTS requires 'I' (iodide), but reaction has:")
        if 'Br' in smiles:
            print(f"      - 'Br' (bromide) instead")
        elif 'Cl' in smiles:
            print(f"      - 'Cl' (chloride) instead")
    
    if num_no_filter > 0:
        top = result_no_filter['protocols'][0]
        print(f"  Top match (by DRFP): {top.get('family', 'Unknown')} (sim={top.get('similarity', 0):.3f})")
    
    print()

print("=" * 80)
print("CONCLUSION:")
print("=" * 80)
print()
print("The Sonogashira protocol in the database has SMARTS pattern:")
print(f"  '{protocol_smarts}'")
print()
print("This pattern REQUIRES:")
print("  1. Aryl/heteroaryl IODIDE ([c,n,o,s]I)")
print("  2. Terminal alkyne (C#C)")
print()
print("But the test reactions use:")
print("  - Aryl BROMIDE (Brc1ccccc1) - doesn't match [c,n,o,s]I")
print("  - Aryl CHLORIDE (Clc1ccc...) - doesn't match [c,n,o,s]I")
print()
print("SOLUTIONS:")
print("  1. Update protocol SMARTS to accept Br/Cl: '[c,n,o,s][Br,Cl,I].C#C>>[c,n,o,s]C#C'")
print("  2. Or add separate protocols for Br/Cl variants")
print("  3. Or use --no-smarts-filter to see similarity matches")
