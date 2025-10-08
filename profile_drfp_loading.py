"""
Profile exactly where the time is spent during precedent search with DRFP.
"""

import time
from chemtools import precedent, reaction_similarity as rs

print("=" * 70)
print("DRFP LOADING PERFORMANCE PROFILE")
print("=" * 70)

# Test reaction
test_reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
family = "C_N_Coupling_Pd"
k = 10

print(f"\nTest reaction: {test_reaction}")
print(f"Family: {family}")
print(f"K: {k}")
print()

# Step 1: Load dataset
print("[1] Loading dataset...")
start = time.perf_counter()

import json
import os

dataset_path = f"data/reaction_dataset/{family}.jsonl"
if not os.path.exists(dataset_path):
    dataset_path = f"C:/Git-softwares/Condition-agent/data/reaction_dataset/{family}.jsonl"

rows = []
with open(dataset_path, 'r', encoding='utf-8') as f:
    for line in f:
        if line.strip():
            rows.append(json.loads(line))

elapsed = time.perf_counter() - start
print(f"   ✓ Loaded {len(rows)} reactions in {elapsed:.3f}s ({elapsed/len(rows)*1000:.2f}ms per reaction)")
print()

# Step 2: Extract reaction SMILES
print("[2] Extracting reaction SMILES from dataset...")
start = time.perf_counter()

reaction_smiles_list = []
for row in rows[:100]:  # Test with first 100
    rsmi = row.get("precomputed", {}).get("reaction_smiles") or row.get("smiles", {})
    if isinstance(rsmi, dict):
        reactants = rsmi.get("reactants", "")
        products = rsmi.get("products", "")
        rsmi = f"{reactants}>>{products}"
    reaction_smiles_list.append(rsmi)

elapsed = time.perf_counter() - start
print(f"   ✓ Extracted {len(reaction_smiles_list)} reaction SMILES in {elapsed:.3f}s")
print()

# Step 3: Generate DRFP fingerprints (the suspected bottleneck)
print("[3] Generating DRFP fingerprints...")
print(f"   Testing on first 100 reactions...")
start = time.perf_counter()

fps = []
for i, rsmi in enumerate(reaction_smiles_list, 1):
    fp_start = time.perf_counter()
    fp = rs.encode_drfp(rsmi)
    fp_elapsed = time.perf_counter() - fp_start
    
    fps.append(fp)
    
    if i % 10 == 0:
        avg_time = (time.perf_counter() - start) / i
        print(f"   Progress: {i}/100 reactions, avg {avg_time*1000:.1f}ms per reaction")

total_elapsed = time.perf_counter() - start
avg_per_reaction = total_elapsed / len(fps)

print()
print(f"   ✓ Generated {len(fps)} fingerprints in {total_elapsed:.3f}s")
print(f"   ✓ Average: {avg_per_reaction*1000:.1f}ms per reaction")
print()

# Estimate for full dataset
print("[4] Estimates for full dataset:")
print(f"   Full dataset: {len(rows)} reactions")
print(f"   Estimated time: {avg_per_reaction * len(rows):.1f}s ({avg_per_reaction * len(rows) / 60:.1f} minutes)")
print()

# Step 4: Test cached access
print("[5] Testing cached DRFP access (second call)...")
start = time.perf_counter()

fps_cached = []
for rsmi in reaction_smiles_list[:20]:
    fp = rs.encode_drfp_cached(rsmi)  # Should be cached now
    fps_cached.append(fp)

elapsed = time.perf_counter() - start
print(f"   ✓ Retrieved {len(fps_cached)} cached fingerprints in {elapsed:.3f}s ({elapsed/len(fps_cached)*1000:.2f}ms per reaction)")
print()

print("=" * 70)
print("ANALYSIS:")
print("=" * 70)
print(f"First-time DRFP generation: {avg_per_reaction*1000:.1f}ms per reaction")
print(f"Cached DRFP retrieval: {elapsed/len(fps_cached)*1000:.2f}ms per reaction")
print(f"Speedup from caching: {(avg_per_reaction/(elapsed/len(fps_cached))):.1f}x")
print()
print(f"For {len(rows)} reactions (first time): ~{avg_per_reaction * len(rows):.0f}s")
print(f"For {len(rows)} reactions (cached): ~{(elapsed/len(fps_cached)) * len(rows):.0f}s")
