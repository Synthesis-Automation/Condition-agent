import numpy as np

# Test: Does _cosine_similarity work with list from query and ndarray from protocol?
print("Testing _cosine_similarity with mixed types...")
print()

# Simulate what happens in recommend.py:
# 1. Query DRFP is returned as a list
query_drfp_list = [0.1] * 2048  # Simulate list return from _compute_drfp

# 2. Protocol DRFP is loaded as ndarray from NPZ
protocol_drfp_ndarray = np.random.rand(2048)

print(f"Query type: {type(query_drfp_list)}, len: {len(query_drfp_list)}")
print(f"Protocol type: {type(protocol_drfp_ndarray)}, shape: {protocol_drfp_ndarray.shape}")
print()

# 3. Test cosine similarity
def cosine_similarity(vec1, vec2):
    v1 = np.array(vec1)
    v2 = np.array(vec2)
    
    print(f"After np.array conversion:")
    print(f"  v1 shape: {v1.shape}, dtype: {v1.dtype}")
    print(f"  v2 shape: {v2.shape}, dtype: {v2.dtype}")
    
    dot_product = np.dot(v1, v2)
    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)
    
    if norm1 == 0 or norm2 == 0:
        return 0.0
    
    similarity = dot_product / (norm1 * norm2)
    return float(similarity)

try:
    result = cosine_similarity(query_drfp_list, protocol_drfp_ndarray)
    print()
    print(f"Similarity: {result:.4f}")
    print("SUCCESS!")
except Exception as e:
    print()
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()
