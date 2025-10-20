from chemtools.protocol.recommend import ProtocolRecommender
from chemtools.protocol.indexer import ProtocolIndexer
import numpy as np

# Test
reaction_smiles = "Ic1ccncc1.C#CC>>C#Cc1ccncc1"

print("=" * 80)
print("MANUAL SIMILARITY COMPUTATION TEST")
print("=" * 80)
print()

recommender = ProtocolRecommender()

# Get the Sonogashira protocol
sonogashira_record = recommender.indexer.records.get('Sonogashira-Coupling.json')

if not sonogashira_record:
    print("ERROR: Sonogashira protocol not found!")
    exit(1)

print(f"Protocol: {sonogashira_record.filename}")
print(f"  Family: {sonogashira_record.reaction_family}")
print(f"  SMARTS: {sonogashira_record.reaction_smarts}")
print()

# Compute query DRFP
print("Computing query DRFP...")
query_drfp = recommender._compute_drfp(reaction_smiles)
print(f"  Query DRFP shape: {query_drfp.shape if hasattr(query_drfp, 'shape') else 'N/A'}")
print(f"  Query DRFP type: {type(query_drfp)}")
print(f"  Query DRFP non-zero: {np.count_nonzero(query_drfp) if hasattr(query_drfp, 'shape') else 'N/A'}")
print()

# Get protocol DRFP
print("Loading protocol DRFP...")
protocol_drfp = recommender.indexer.get_drfp_fingerprint(sonogashira_record.filename)
print(f"  Protocol DRFP shape: {protocol_drfp.shape if protocol_drfp is not None else 'N/A'}")
print(f"  Protocol DRFP type: {type(protocol_drfp)}")
print(f"  Protocol DRFP non-zero: {np.count_nonzero(protocol_drfp) if protocol_drfp is not None else 'N/A'}")
print()

# Compute similarity
if query_drfp is not None and protocol_drfp is not None:
    print("Computing cosine similarity...")
    
    # Check if they're the same length
    if hasattr(query_drfp, 'shape') and hasattr(protocol_drfp, 'shape'):
        print(f"  Query length: {query_drfp.shape[0]}")
        print(f"  Protocol length: {protocol_drfp.shape[0]}")
        
        if query_drfp.shape[0] != protocol_drfp.shape[0]:
            print(f"  ERROR: Shape mismatch!")
        else:
            try:
                similarity = recommender._cosine_similarity(query_drfp, protocol_drfp)
                print(f"  Similarity: {similarity:.4f}")
            except Exception as e:
                print(f"  ERROR computing similarity: {e}")
                import traceback
                traceback.print_exc()
else:
    print("ERROR: One or both DRFPs are None!")

print()
print("=" * 80)
