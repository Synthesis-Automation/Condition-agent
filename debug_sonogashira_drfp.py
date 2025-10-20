from chemtools.protocol.recommend import ProtocolRecommender
import logging

# Set up detailed logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(levelname)s: %(message)s'
)

# Test
reaction_smiles = "Ic1ccncc1.C#CC>>C#Cc1ccncc1"

print("=" * 80)
print("DETAILED DEBUG: IODIDE SONOGASHIRA WITH SMARTS FILTER")
print("=" * 80)
print()

recommender = ProtocolRecommender()

# Manually check the Sonogashira protocol
print("Checking Sonogashira-Coupling.json protocol:")
record = recommender.indexer.records.get('Sonogashira-Coupling.json')
if record:
    print(f"  Filename: {record.filename}")
    print(f"  Family: {record.reaction_family}")
    print(f"  SMARTS: {record.reaction_smarts}")
    print(f"  Has DRFP: {record.has_drfp}")
    
    # Try to get DRFP
    drfp = recommender.indexer.get_drfp_fingerprint(record.filename)
    print(f"  DRFP loaded: {drfp is not None}")
    if drfp is not None:
        print(f"  DRFP shape: {drfp.shape if hasattr(drfp, 'shape') else 'N/A'}")
        print(f"  DRFP type: {type(drfp)}")
    else:
        print("  ERROR: Could not load DRFP!")
else:
    print("  ERROR: Sonogashira-Coupling.json not found in records!")

print()
print("=" * 80)
print()

# Now test the recommendation
print("Running recommendation with SMARTS filter...")
result = recommender.recommend(reaction_smiles, k=3, use_smarts_filter=True)

print()
print(f"Protocols returned: {len(result.get('protocols', []))}")

if result.get('protocols'):
    for p in result['protocols']:
        print(f"  - {p.get('family')}: similarity={p.get('similarity'):.3f}")
else:
    print("  None!")
    print()
    print(f"Extras: {result.get('extras', {})}")
